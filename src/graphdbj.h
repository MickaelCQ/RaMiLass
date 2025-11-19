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

/// Boilerplate de configuration pour GraphDBJ
struct GraphDBJConfig {
    // Paramètres d'erreur pour la fusion des contigs (pour calculateBinaryOverlap / isContained)
    double ERROR_PERCENT_OVERLAP = 0.05; // 5% de tolérance pour le chevauchement (min_mismatches / len)
    double ERROR_PERCENT_CONTAINED = 0.02; // 2% de tolérance pour l'inclusion

    // Paramètres de nettoyage du graphe
    size_t TOPO_MAX_LEN;    ///< Longueur max (en noeuds) pour suppression topologique (tip)
    size_t RCTC_MAX_LEN;    ///< Longueur max (en noeuds) pour suppression par couverture (tip)
    double RCTC_MAX_RATIO = 10;    ///< Ratio de couverture ancrage/bout
    double TOPO_MAX_RATIO = 2.5;    ///< Ratio de couverture ancrage/bout pour suppression topologique
    double RCTC_RATIO = 2.0;    ///< Ratio de couverture ancrage/bout
    double COVERAGE_RATIO = 1.0; ///< Ratio de couverture pour la résolution des bifurcations (contigs)
    size_t MAX_CONTIG_LEN = 1000000; ///< Longueur max de contig à générer (garde-fou)
    double SEARCH_DEPTH_FACTOR = 20.0; ///< Profondeur de recherche = k * facteur
    int MAX_PASSES = 50;               ///< Nombre maximum de cycles de simplification

    // Le constructeur est utile pour initialiser les valeurs dépendantes de k
    GraphDBJConfig(int k_size) {
        // Initialisation des valeurs en fonction de k (k-mer size)
        TOPO_MAX_LEN = (size_t)(k_size * TOPO_MAX_RATIO);
        RCTC_MAX_LEN = (size_t)(k_size * RCTC_MAX_RATIO);
    }
};


class GraphDBJ {
private:
    int k;
    std::unordered_map<uint64_t, Noeud*> nodes_map;
    GraphDBJConfig config;

    uint64_t getReverseComplement(uint64_t val, int k) const;
    uint64_t extractKmerValue(const BitVector& bv, size_t start_bit_idx, int len_nucleotides) const;
    std::string kmerToString(uint64_t val, int length) const;
    void disconnectNodes(Noeud* parent, Noeud* child);
    static Noeud* findConvergence(Noeud* branch1, Noeud* branch2, int depth_limit);
    void addKmerToBitVector(BitVector& bv, uint64_t val, int length) const;
    static int calculateBinaryOverlap(const BitVector& b1, const BitVector& b2, int min_len, double error_percent);
    static bool isContained(const BitVector& large, const BitVector& small, double error_percent);

public:
    GraphDBJ(const Convert& converter, int kmer_size, const GraphDBJConfig& config);
    ~GraphDBJ();

    std::vector<Noeud*> getNodes() const;

    // --- Algorithmes d'Assemblage ---
    int clipTips();
    int resolveBubbles();
    std::vector<BitVector> generateContigs() const;
    void exportToGFA(const std::string& filename) const;
    static std::vector<BitVector> mergeContigs(std::vector<BitVector> contigs, int min_overlap, double overlap_error_percent, double contained_error_percent);
};

#endif //ASSEMBLEUR_GRAPHDBJ_H