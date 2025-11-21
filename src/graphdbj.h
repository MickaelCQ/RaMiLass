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

    int clipTips();
    int resolveBubbles();
    std::vector<BitVector> generateContigs() const;
    void exportToGFA(const std::string& filename) const;

    /**
     * @brief Fusionne les contigs (Mise à jour pour utiliser Config).
     */
    static void mergeContigs(std::vector<BitVector>& contigs, int k, int min_overlap, const GraphDBJConfig& config);

    int removeLowDepthKmers(uint32_t threshold);
};

#endif //ASSEMBLEUR_GRAPHDBJ_H