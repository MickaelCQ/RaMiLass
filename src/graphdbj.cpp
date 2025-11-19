#include "graphdbj.h"

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <algorithm>

// --- Helpers Binaires (Statics) ---

// Calcule le complément inverse d'un BitVector complet
static BitVector getBitVectorReverseComplement(const BitVector& seq) {
    BitVector rc;
    rc.reserve(seq.size());
    size_t len = seq.size() / 2;

    // Parcours inversé des nucléotides
    for (size_t i = 0; i < len; ++i) {
        // Index original : on part de la fin
        size_t original_idx = len - 1 - i;
        uint8_t nuc = seq.getNucleotideAt(original_idx);

        // Complément (XOR 3) : 00<->11, 01<->10
        uint8_t comp = nuc ^ 3;

        // Ajout au nouveau vecteur
        // Bit de poids fort (valeur 2)
        rc.push_back((comp & 2) != 0);
        // Bit de poids faible (valeur 1)
        rc.push_back((comp & 1) != 0);
    }
    return rc;
}

// Calcule l'overlap (chevauchement) en travaillant uniquement sur les bits
static int calculateBinaryOverlap(const BitVector& b1, const BitVector& b2, int min_len) {
    int max_ov = 0;
    size_t len1 = b1.size() / 2;
    size_t len2 = b2.size() / 2;
    size_t limit = std::min(len1, len2);

    for (size_t len = min_len; len <= limit; ++len) {
        int mismatches = 0;
        int allowed_errors = 2 + (len / 30);
        bool match = true;

        // Comparaison des 'len' derniers nucléotides de b1
        // avec les 'len' premiers de b2

        size_t start_idx_b1 = (len1 - len); // index en nucléotides

        for (size_t i = 0; i < len; ++i) {
            // Comparaison bit à bit (2 bits par nucléotide)
            uint8_t n1 = b1.getNucleotideAt(start_idx_b1 + i);
            uint8_t n2 = b2.getNucleotideAt(i);

            if (n1 != n2) {
                mismatches++;
                if (mismatches > allowed_errors) {
                    match = false;
                    break;
                }
            }
        }

        if (match) {
            max_ov = len;
        }
    }
    return max_ov;
}

// Vérifie si 'small' est inclus dans 'large' (recherche de sous-séquence binaire exacte)
static bool isContained(const BitVector& large, const BitVector& small) {
    size_t l_len = large.size() / 2;
    size_t s_len = small.size() / 2;

    // Impossible que le petit contienne le grand
    if (s_len > l_len) return false;

    // Tolérance d'erreur :
    // On autorise 1 erreur de base + 1 erreur tous les 50 nucléotides (~2%)
    // Cela permet de gérer les petites variations ou les bouts mal coupés.
    size_t max_mismatches = 1 + (s_len / 50);

    // Recherche naïve avec tolérance (Sliding Window)
    for (size_t i = 0; i <= l_len - s_len; ++i) {
        size_t mismatches = 0;
        bool possible = true;

        for (size_t j = 0; j < s_len; ++j) {
            // Comparaison des nucléotides
            if (large.getNucleotideAt(i + j) != small.getNucleotideAt(j)) {
                mismatches++;
                // Si on dépasse le seuil d'erreur, on arrête cette fenêtre
                if (mismatches > max_mismatches) {
                    possible = false;
                    break;
                }
            }
        }

        // Si 'possible' est vrai ici, c'est qu'on a trouvé une correspondance
        // avec un nombre d'erreurs acceptable. Le petit est donc "contenu".
        if (possible) return true;
    }

    return false;
}


std::vector<BitVector> GraphDBJ::mergeContigs(std::vector<BitVector> contigs, int min_overlap) {
    std::cout << "--- Post-traitement : Fusion des Contigs ---" << std::endl;

    // Tri par longueur (en bits) décroissante
    std::sort(contigs.begin(), contigs.end(), [](const BitVector& a, const BitVector& b) {
        return a.size() > b.size();
    });

    bool global_change = true;
    while (global_change) {
        global_change = false;
        std::vector<bool> absorbed(contigs.size(), false);

        for (size_t i = 0; i < contigs.size(); ++i) {
            if (absorbed[i]) continue;

            bool local_change = true;
            while (local_change) {
                local_change = false;

                for (size_t j = 0; j < contigs.size(); ++j) {
                    if (i == j || absorbed[j]) continue;

                    BitVector& master = contigs[i];
                    const BitVector& candidate = contigs[j];
                    bool merged = false;

                    // --- TENTATIVE 1 : Orientation Standard ---
                    if (isContained(master, candidate)) {
                        merged = true;
                    }
                    else if (int ov = calculateBinaryOverlap(master, candidate, min_overlap); ov > 0) {
                        // Master -> Candidate (Overlap à la fin de Master)
                        // On ajoute Candidate en sautant les 'ov' premiers nucléotides
                        master.append(candidate, ov);
                        merged = true;
                    }
                    else if (int ov = calculateBinaryOverlap(candidate, master, min_overlap); ov > 0) {
                        // Candidate -> Master (Overlap à la fin de Candidate)
                        // On crée un nouveau master : Candidate + (Master sans le début)
                        BitVector newMaster = candidate;
                        newMaster.append(master, ov);
                        master = newMaster; // Move assignment
                        merged = true;
                    }

                    // --- TENTATIVE 2 : Orientation Inverse ---
                    if (!merged) {
                        BitVector candRC = getBitVectorReverseComplement(candidate);

                        if (isContained(master, candRC)) {
                            merged = true;
                        }
                        else if (int ov = calculateBinaryOverlap(master, candRC, min_overlap); ov > 0) {
                            master.append(candRC, ov);
                            merged = true;
                        }
                        else if (int ov = calculateBinaryOverlap(candRC, master, min_overlap); ov > 0) {
                            BitVector newMaster = candRC;
                            newMaster.append(master, ov);
                            master = newMaster;
                            merged = true;
                        }
                    }

                    if (merged) {
                        absorbed[j] = true;
                        local_change = true;
                        global_change = true;
                        break;
                    }
                }
            }
        }

        if (global_change) {
            std::vector<BitVector> next_pass;
            for (size_t k = 0; k < contigs.size(); ++k) {
                if (!absorbed[k]) next_pass.push_back(contigs[k]);
            }
            contigs = next_pass;
            std::cout << "Cycle termine. Contigs restants : " << contigs.size() << std::endl;
        }
    }
    return contigs;
}

void GraphDBJ::addKmerToBitVector(BitVector& bv, uint64_t val, int length) const {
    // Reconstitution binaire (poids fort = premier nuc)
    for (int i = length - 1; i >= 0; --i) {
        uint64_t mask = 3ULL << (i * 2);
        uint64_t two_bits = (val & mask) >> (i * 2);

        // two_bits est 0(A), 2(C), 1(G), 3(T) selon l'encodage uint64
        // Mapping vers BitVector::push_back :
        // 0 (00) -> false, false
        // 2 (10) -> true, false
        // 1 (01) -> false, true
        // 3 (11) -> true, true

        bv.push_back((two_bits & 2) != 0); // Bit de poids fort
        bv.push_back((two_bits & 1) != 0); // Bit de poids faible
    }
}

std::vector<BitVector> GraphDBJ::generateContigs() const {
    std::vector<BitVector> contigs;
    std::unordered_map<Noeud*, bool> visited;
    const double COVERAGE_RATIO = 1.0;
    const size_t MAX_CONTIG_LEN = 1000000;

    for (const auto& pair : nodes_map) {
        Noeud* startNode = pair.second;

        if (startNode->removed || visited[startNode]) continue;

        if (startNode->parents.empty() || startNode->parents.size() > 1) {

            // INITIALISATION DU CONTIG EN BINAIRE
            BitVector currentContig;
            // On ajoute la séquence complète du premier k-1 mer
            addKmerToBitVector(currentContig, startNode->p, k - 1);

            visited[startNode] = true;
            Noeud* curr = startNode;
            size_t sanity_check = 0;

            while (sanity_check < MAX_CONTIG_LEN) {
                sanity_check++;
                Noeud* next = nullptr;

                // Choix du prochain noeud (inchangé)
                if (curr->c.size() == 1) {
                    next = curr->c[0];
                }
                else if (curr->c.size() > 1) {
                    Noeud* best_candidate = nullptr;
                    uint32_t max_cov = 0;
                    uint32_t second_max_cov = 0;

                    for (auto* child : curr->c) {
                        if (child->removed) continue;
                        if (child->coverage > max_cov) {
                            second_max_cov = max_cov;
                            max_cov = child->coverage;
                            best_candidate = child;
                        } else if (child->coverage > second_max_cov) {
                            second_max_cov = child->coverage;
                        }
                    }
                    if (best_candidate != nullptr && max_cov > (second_max_cov * COVERAGE_RATIO)) {
                        next = best_candidate;
                    } else { break; }
                } else { break; }

                if (next == nullptr || next->removed) break;

                // EXTRACTION ET AJOUT BINAIRE
                uint64_t val = next->p;
                uint64_t last_nuc_val = val & 3ULL; // Les 2 derniers bits

                // Ajout direct au BitVector
                currentContig.push_back((last_nuc_val & 2) != 0);
                currentContig.push_back((last_nuc_val & 1) != 0);

                if (visited[next]) break;

                visited[next] = true;
                curr = next;
            }
            contigs.push_back(currentContig);
        }
    }
    return contigs;
}

// ... (Le reste des méthodes : constructeur, resolveBubbles, clipTips restent inchangées
// car elles utilisent déjà uint64_t ou la structure de graphe) ...
// Assurez-vous d'inclure le reste du fichier original ici (constructeur, exportToGFA, etc.)

// (Je remets les méthodes non modifiées nécessaires à la compilation pour que le fichier soit complet)

GraphDBJ::GraphDBJ(const Convert& converter, int kmer_size) : k(kmer_size) {
    if (k < 2 || k > 32) throw std::invalid_argument("K invalide (2-32)");
    const BitVector& bv = converter.get_bitVector();
    const std::vector<size_t>& read_ends = converter.get_read_end_positions();
    size_t current_read_start = 0;

    for (size_t end_pos : read_ends) {
        size_t read_len_nuc = (end_pos - current_read_start) / 2;
        if (read_len_nuc >= (size_t)k) {
            for (size_t i = 0; i <= read_len_nuc - k; ++i) {
                uint64_t u_fwd = extractKmerValue(bv, current_read_start + i * 2, k - 1);
                uint64_t v_fwd = extractKmerValue(bv, current_read_start + i * 2 + 2, k - 1);
                uint64_t v_rev = getReverseComplement(v_fwd, k - 1);
                uint64_t u_rev = getReverseComplement(u_fwd, k - 1);

                auto add_edge = [&](uint64_t source_val, uint64_t target_val) {
                    Noeud*& sourceNode = nodes_map[source_val];
                    if (!sourceNode) sourceNode = new Noeud(source_val);
                    Noeud*& targetNode = nodes_map[target_val];
                    if (!targetNode) targetNode = new Noeud(target_val);
                    sourceNode->coverage++; targetNode->coverage++;
                    bool edgeExists = false;
                    for (auto* child : sourceNode->c) { if (child == targetNode) { edgeExists = true; break; } }
                    if (!edgeExists) { sourceNode->c.push_back(targetNode); targetNode->parents.push_back(sourceNode); }
                };
                add_edge(u_fwd, v_fwd);
                add_edge(v_rev, u_rev);
            }
        }
        current_read_start = end_pos;
    }
}

GraphDBJ::~GraphDBJ() {
    for (auto& pair : nodes_map) delete pair.second;
    nodes_map.clear();
}

std::vector<Noeud*> GraphDBJ::getNodes() const {
    std::vector<Noeud*> result;
    for (const auto& pair : nodes_map) result.push_back(pair.second);
    return result;
}

// Note: extractKmerValue reste identique car elle convertit BitVector -> uint64_t
uint64_t GraphDBJ::extractKmerValue(const BitVector& bv, size_t start_bit_idx, int len_nucleotides) const {
    uint64_t val = 0;
    for (int i = 0; i < len_nucleotides; ++i) {
        size_t pos = start_bit_idx + (i * 2);
        bool b1 = bv.test(pos);
        bool b2 = bv.test(pos + 1);
        val = val << 2;
        if (b1) val |= 2ULL;
        if (b2) val |= 1ULL;
    }
    return val;
}

// Note: getReverseComplement (uint64) reste identique
uint64_t GraphDBJ::getReverseComplement(uint64_t val, int length) const {
    uint64_t res = 0;
    for (int i = 0; i < length; ++i) {
        uint64_t nuc = val & 3ULL;
        uint64_t comp = nuc ^ 3ULL;
        res |= (comp << (2 * (length - 1 - i)));
        val >>= 2;
    }
    return res;
}

// Les méthodes clipTips et resolveBubbles sont identiques à votre version,
// simplement s'assurer d'inclure disconnectNodes et findConvergence.
void GraphDBJ::disconnectNodes(Noeud* parent, Noeud* child) {
    if (!parent || !child) return;
    parent->c.erase(std::remove(parent->c.begin(), parent->c.end(), child), parent->c.end());
    child->parents.erase(std::remove(child->parents.begin(), child->parents.end(), parent), child->parents.end());
}

int GraphDBJ::clipTips() {
    // Inspiré par la méthode utilisé par minia
    int tips_removed_count = 0;
    bool changed = true;
    const int TOPO_MAX_LEN = (int)(k * 2.5);
    const int RCTC_MAX_LEN = k * 10;
    const double RCTC_RATIO = 2.0;

    while (changed) {
        changed = false;
        std::vector<uint64_t> keys;
        keys.reserve(nodes_map.size());
        for(auto& p : nodes_map) keys.push_back(p.first);

        for (uint64_t key : keys) {
            if (nodes_map.find(key) == nodes_map.end()) continue;
            Noeud* node = nodes_map[key];
            if (node->removed) continue;
            if (node->c.empty()) {
                std::vector<Noeud*> tip_path;
                Noeud* curr = node;
                bool is_valid_tip = true;
                uint64_t sum_coverage = 0;
                while (tip_path.size() <= RCTC_MAX_LEN) {
                    tip_path.push_back(curr);
                    sum_coverage += curr->coverage;
                    if (curr->parents.empty()) { is_valid_tip = false; break; }
                    Noeud* parent = curr->parents[0];
                    if (parent->c.size() > 1) break;
                    if (parent->parents.size() != 1) { is_valid_tip = false; break; }
                    curr = parent;
                }
                if (!is_valid_tip || tip_path.empty()) continue;
                Noeud* last_tip_node = tip_path.back();
                if (last_tip_node->parents.empty()) continue;
                Noeud* anchor = last_tip_node->parents[0];
                int len = tip_path.size();
                double tip_avg_cov = (double)sum_coverage / len;
                double anchor_cov = (double)anchor->coverage;
                bool should_remove = false;
                if (len <= TOPO_MAX_LEN) should_remove = true;
                else if (len <= RCTC_MAX_LEN && anchor_cov > (tip_avg_cov * RCTC_RATIO)) should_remove = true;

                if (should_remove) {
                    disconnectNodes(anchor, last_tip_node);
                    for (Noeud* n : tip_path) n->removed = true;
                    tips_removed_count++;
                    changed = true;
                }
            }
        }
    }
    return tips_removed_count;
}

int GraphDBJ::resolveBubbles() {
    // ... (Code identique à votre version fournie)
    // Reprenez le corps de resolveBubbles du fichier original.
    int search_depth = k * 20;
    int bubbles_collapsed = 0;
    bool changed = true;
    int max_passes = 50;
    int pass = 0;
    while (changed && pass < max_passes) {
        changed = false; pass++;
        std::vector<uint64_t> keys;
        keys.reserve(nodes_map.size());
        for(auto& p : nodes_map) keys.push_back(p.first);
        for (uint64_t key : keys) {
            if (nodes_map.find(key) == nodes_map.end()) continue;
            Noeud* s = nodes_map[key];
            if (s->removed || s->c.size() != 2) continue;
            Noeud* branch1 = s->c[0]; Noeud* branch2 = s->c[1];
            std::vector<Noeud*> path1, path2;
            path1.push_back(branch1); path2.push_back(branch2);
            Noeud* convergence = nullptr;
            for (int step = 0; step < search_depth; ++step) {
                Noeud* last1 = path1.back();
                if (last1->c.size() == 1 && !last1->c[0]->removed) path1.push_back(last1->c[0]);
                Noeud* last2 = path2.back();
                if (last2->c.size() == 1 && !last2->c[0]->removed) path2.push_back(last2->c[0]);
                for(Noeud* n2 : path2) { if (path1.back() == n2) { convergence = n2; break; } }
                if(convergence) break;
                for(Noeud* n1 : path1) { if (path2.back() == n1) { convergence = n1; break; } }
                if(convergence) break;
                if (path1.back() == last1 && path2.back() == last2) break;
            }
            if (convergence && convergence != branch1 && convergence != branch2) {
                uint32_t cov1 = 0; for(auto* n : path1) { if(n==convergence) break; cov1 += n->coverage; }
                uint32_t cov2 = 0; for(auto* n : path2) { if(n==convergence) break; cov2 += n->coverage; }
                std::vector<Noeud*>& to_remove = (cov1 >= cov2) ? path2 : path1;
                Noeud* branch_to_keep = (cov1 >= cov2) ? branch1 : branch2;
                for(auto* n : to_remove) { if (n != convergence) n->removed = true; }
                Noeud* branch_to_remove = (cov1 >= cov2) ? branch2 : branch1;
                disconnectNodes(s, branch_to_remove);
                bool connected = false;
                for(auto* child : s->c) if(child == branch_to_keep) connected = true;
                if(!connected) s->c.push_back(branch_to_keep);
                bubbles_collapsed++; changed = true;
            }
        }
    }
    return bubbles_collapsed;
}

std::string GraphDBJ::kmerToString(uint64_t val, int length) const {
    std::string seq = "";
    for (int i = length - 1; i >= 0; --i) {
        uint64_t mask = 3ULL << (i * 2);
        uint64_t two_bits = (val & mask) >> (i * 2);
        if (two_bits == 0) seq += 'A';
        else if (two_bits == 2) seq += 'C';
        else if (two_bits == 1) seq += 'G';
        else seq += 'T';
    }
    return seq;
}

void GraphDBJ::exportToGFA(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) return;
    out << "H\tVN:Z:1.0\n";
    std::unordered_map<uint64_t, int> id_map;
    int current_id = 1;
    for (const auto& pair : nodes_map) {
        if (!pair.second->removed) id_map[pair.first] = current_id++;
    }
    for (const auto& pair : nodes_map) {
        Noeud* node = pair.second;
        if (node->removed) continue;
        std::string seq = kmerToString(node->p, k - 1);
        out << "S\t" << id_map[node->p] << "\t" << seq << "\tDP:i:" << node->coverage << "\n";
    }
    for (const auto& pair : nodes_map) {
        Noeud* src = pair.second;
        if (src->removed) continue;
        for (Noeud* dest : src->c) {
            if (dest->removed) continue;
            if (id_map.find(src->p) != id_map.end() && id_map.find(dest->p) != id_map.end()) {
                out << "L\t" << id_map[src->p] << "\t+\t" << id_map[dest->p] << "\t+\t" << (k - 2) << "M\n";
            }
        }
    }
    out.close();
}