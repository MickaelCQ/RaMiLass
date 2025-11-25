#ifndef ASSEMBLEUR_GRAPHDBJ_H
#define ASSEMBLEUR_GRAPHDBJ_H

#include <vector>
#include <cstdint>
#include <unordered_map>
#include <string>
#include <algorithm>
#include "convert.h"

struct Noeud {
    uint64_t p;
    std::vector<Noeud*> c;
    std::vector<Noeud*> parents;
    bool removed = false;
    uint32_t coverage = 0;
    Noeud(uint64_t val) : p(val) {}
};

struct GraphDBJConfig {
    // --- Paramètres de Fusion (Post-traitement) ---
    double ERROR_PERCENT_OVERLAP = 0.05;
    double ERROR_PERCENT_CONTAINED = 0.02;

    // --- NOUVEAUX PARAMETRES ---
    size_t MAX_SCAN_DEPTH = 5000; ///< Profondeur max de scan pour l'extension (en pb)
    size_t MAX_SEED_DEPTH = 1500; ///< Profondeur max pour chercher une graine (seed)

    // --- Paramètres de Nettoyage du Graphe ---
    size_t TOPO_MAX_LEN;
    size_t RCTC_MAX_LEN;

    double TOPO_MAX_RATIO = 2.5;
    double RCTC_MAX_RATIO = 5;

    uint32_t MIN_DEPTH = 1;
    double COVERAGE_RATIO = 1.0;

    size_t MAX_CONTIG_LEN = 1000000;
    double SEARCH_DEPTH_FACTOR = 20.0;
    int MAX_PASSES = 50;

    GraphDBJConfig(int k_size) {
        TOPO_MAX_LEN = (size_t)(k_size * TOPO_MAX_RATIO);
        RCTC_MAX_LEN = (size_t)(k_size * RCTC_MAX_RATIO);
    }
};

class GraphDBJ {
private:
    // ... (Membres privés inchangés) ...
    int k;
    std::unordered_map<uint64_t, Noeud*> nodes_map;
    GraphDBJConfig config;

    /**
     * @brief Retourne complément inverse k-mer codé sur 2*length bits.
     * @complexity : time : O(k), space : O(1).
     * k = taille du k-mer.
     */
    uint64_t getReverseComplement(uint64_t val, int k) const;

    /**
     * @brief Extrait valeur d'un k-mer depuis le BitVector.
     * @complexity : time : O(ln), space : O(1).
     * ln = nb de nucléotides à extraire (tend vers k).
     */
    uint64_t extractKmerValue(const BitVector& bv, size_t start_bit_idx, int len_nucleotides) const;

    /**
     * @brief Convertit un entier k-mer en chaîne de caractères.
     * @complexity : time :  O(k), space : O(k).
     * k = taille du k-mer.
     */
    std::string kmerToString(uint64_t val, int length) const;

    /**
     * @brief Déconnecte des nœuds.
     * @complexity : time : O(dp + dc), space : O(1).
     * dp = degré parent, dc = degré enfant.
     */
    void disconnectNodes(Noeud* parent, Noeud* child);

    /**
     * @brief Trouve une convergence.
     * @complexity : time : O(d^dc), Space : O(d^dc).
     * d = depth_limit (profondeur max de recherche), dc = degrè noeud (pire des cas sinon d).
     */
    static Noeud* findConvergence(Noeud* branch1, Noeud* branch2, int depth_limit);

    /**
     * @brief Ajoute au BitVector un kmer.
     * @complexity : time : O(k), space : O(1).
     * k = taille du k-mer. espace constant car bitvector gère d'avance l'allocation.
     */
    void addKmerToBitVector(BitVector& bv, uint64_t val, int length) const;

    /**
     * @brief Calcule la longueur de chevauchement entre deux BitVectors).
     * @complexity : time : O(lbmin + lbmax), space : O(1).
     * lbmin = longueur du plus petit bitvector. lbmax = longueur du plus grand bitvector (généralement lb * 2).
     */
    static int calculateBinaryOverlap(const BitVector& b1, const BitVector& b2, int min_len, double error_percent);

    /**
     * @brief observe si un bitvector et contenu dans un autre.
     * @complexity :time : O(lbmin + lb max), space : O(1).
     * lbmin = longueur du plus petit bitvector. lbmax = longueur du plus grand bitvector. (pire cas peut aller a lbmin * lbmax).
     */
    static bool isContained(const BitVector& large, const BitVector& small, double error_percent);

public:
    /**
     * @brief Construit le graphe de De Bruijn.
     * @complexity : time : O(N * k), space : O(V + E).
     * N = nombre total de nucléotides, k = taille du k-mer (simplifiable à O(N)).
     * V = nombre de nœuds, E = nombre d'arêtes.
     */
    GraphDBJ(const Convert& converter, int kmer_size, const GraphDBJConfig& config);

    /**
     * @brief Destructeur.
     * @complexity : time ; O(V), space : O(1).
     * V = nombre de nœuds.
     */
    ~GraphDBJ();

    /**
     * @brief récupére liste des Noeud*.
     * @complexity : time : O(V), space : O(V).
     * v = nombre de nœuds.
     */
    std::vector<Noeud*> getNodes() const;

    /**
     * @brief Coupe les bras terminaux courts du graphe.
     * @complexity : time : O(p * V * t), space : O(R).
     * V = nombre de nœuds, t = longueur max d'un bout terminale. p = nb de passages. R = seuil de taille max (RCTC_MAX_LEN).
     */
    int clipTips();

    /**
     * @brief Résout bulles ( = bifurcations convergentes) du graphe.
     * @complexity : time : O(p * V * D^2), space : O(D).
     * p = nb de passages, V = nombre de nœuds, D = profondeur max de recherche (pire cas en temps v^D).
     */
    int resolveBubbles();

    /**
     * @brief Génère les contigs en parcourant les chemins linéaires du graphe.
     * @complexity : time : O(V + E), space : O(M)
     * V = nombre de nœuds, E = nombre d'arêtes, M = taille totale des contigs.
     */
    std::vector<BitVector> generateContigs() const;

    /**
     * @brief Exporte le graphe au format GFA.
     * @complexity : time : O(V * k + E), space : O(V).
     * v = nombre de nœuds, E = nombre d'arêtes, k = taille du k-mer.
     */
    void exportToGFA(const std::string& filename) const;

    /**
     * @brief Fusionne les contiges (Mise à jour pour utiliser Config).
     * @complexity : time : O(M * m + S * D * sd), space : O(M / ss)
     * M = taille totale des contigs, m = taille moyenne d'un contig,
     * S = nombre de seeds, D = profondeur max de recherche,
     * sd = profondeur max de scan pour l'extension.
     * ss = intervalle entre les seeds.
     */
    static void mergeContigs(std::vector<BitVector>& contigs, int k, int min_overlap, const GraphDBJConfig& config);

    /**
     * @brief Supprime les k-mers sous un seuil de couverture.
     * @complexity : time : O(V + E), space : O(1).
     * V = nombre de nœuds, E = nombre d'arêtes.
     */
    int removeLowDepthKmers(uint32_t threshold);
};

#endif //ASSEMBLEUR_GRAPHDBJ_H
