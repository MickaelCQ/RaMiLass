#include "graphdbj.h"
#include <iostream>
#include <stdexcept>

GraphDBJ::GraphDBJ(const Convert& converter, int kmer_size) : k(kmer_size) {
    if (k <= 1) {
        throw std::invalid_argument("K doit etre superieur a 1 pour construire un graphe.");
    }
    // Un uint64_t peut stocker jusqu'à 32 nucléotides (64 bits). 
    // Donc k-1 doit être <= 32 => k <= 33.
    if (k > 33) {
        throw std::invalid_argument("Taille de K trop grande pour etre stockee sur 64 bits (max 33).");
    }

    const BitVector& bv = converter.get_bitVector();
    const std::vector<size_t>& read_ends = converter.get_read_end_positions();

    size_t current_read_start = 0;

    // 1. Parcourir chaque lecture
    for (size_t end_pos : read_ends) {
        // Calculer la longueur de la lecture en bits
        size_t read_len_bits = end_pos - current_read_start;
        size_t read_len_nuc = read_len_bits / 2;

        // Si la lecture est assez longue pour contenir au moins un k-mer
        if (read_len_nuc >= (size_t)k) {
            // Parcourir tous les k-mers de cette lecture
            // On s'arrête quand il ne reste plus assez de place pour un k-mer
            for (size_t i = 0; i <= read_len_nuc - k; ++i) {
                
                // Position en bits du début du k-mer courant
                size_t bit_idx = current_read_start + (i * 2);

                // --- CONSTRUCTION NOEUD PREFIXE (k-1) ---
                // Prend k-1 nucléotides à partir de bit_idx
                uint64_t prefix_val = extractKmerValue(bv, bit_idx, k - 1);

                // --- CONSTRUCTION NOEUD SUFFIXE (k-1) ---
                // Le suffixe commence 1 nucléotide (2 bits) après le préfixe
                // et a aussi une longueur de k-1
                uint64_t suffix_val = extractKmerValue(bv, bit_idx + 2, k - 1);

                // 1. Récupérer ou créer le noeud Prefix
                Noeud* prefixNode;
                if (nodes_map.find(prefix_val) == nodes_map.end()) {
                    prefixNode = new Noeud(prefix_val);
                    nodes_map[prefix_val] = prefixNode;
                } else {
                    prefixNode = nodes_map[prefix_val];
                }

                // 2. Récupérer ou créer le noeud Suffix
                Noeud* suffixNode;
                if (nodes_map.find(suffix_val) == nodes_map.end()) {
                    suffixNode = new Noeud(suffix_val);
                    nodes_map[suffix_val] = suffixNode;
                } else {
                    suffixNode = nodes_map[suffix_val];
                }

                // 3. Créer l'arête (ajouter suffixe dans les enfants du préfixe)
                // Note : On peut ajouter une vérification ici pour éviter les doublons d'arêtes
                // si on veut un multigraphe ou non. Ici on ajoute simplement.
                prefixNode->c.push_back(suffixNode);
            }
        }

        // Mise à jour pour la prochaine lecture
        current_read_start = end_pos;
    }
}

GraphDBJ::~GraphDBJ() {
    // Nettoyage de la mémoire : delete de tous les pointeurs stockés dans la map
    for (auto& pair : nodes_map) {
        delete pair.second;
    }
    nodes_map.clear();
}

uint64_t GraphDBJ::extractKmerValue(const BitVector& bv, size_t start_bit_idx, int len_nucleotides) const {
    uint64_t val = 0;

    // On parcourt nucléotide par nucléotide
    for (int i = 0; i < len_nucleotides; ++i) {
        size_t pos = start_bit_idx + (i * 2);
        
        // Lecture des 2 bits
        bool b1 = bv.test(pos);
        bool b2 = bv.test(pos + 1);

        // On décale la valeur actuelle de 2 bits vers la gauche pour faire de la place
        val = val << 2;

        // Construction de la valeur 2 bits (b1b2)
        // Si b1 est vrai, on ajoute 2 (10 en binaire)
        // Si b2 est vrai, on ajoute 1 (01 en binaire)
        if (b1) val |= 2ULL;
        if (b2) val |= 1ULL;
    }

    return val;
}

std::vector<Noeud*> GraphDBJ::getNodes() const {
    std::vector<Noeud*> result;
    result.reserve(nodes_map.size());
    
    for (const auto& pair : nodes_map) {
        result.push_back(pair.second);
    }
    return result;
}