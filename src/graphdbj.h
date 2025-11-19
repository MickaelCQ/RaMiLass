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

class GraphDBJ {
private:
    int k;
    std::unordered_map<uint64_t, Noeud*> nodes_map;

    uint64_t getReverseComplement(uint64_t val, int k) const;
    uint64_t extractKmerValue(const BitVector& bv, size_t start_bit_idx, int len_nucleotides) const;
    std::string kmerToString(uint64_t val, int length) const;
    void disconnectNodes(Noeud* parent, Noeud* child);
    static Noeud* findConvergence(Noeud* branch1, Noeud* branch2, int depth_limit);

    // Helper interne pour reconstruire un noeud vers un bitvector
    void addKmerToBitVector(BitVector& bv, uint64_t val, int length) const;

public:
    GraphDBJ(const Convert& converter, int kmer_size);
    ~GraphDBJ();

    std::vector<Noeud*> getNodes() const;

    // --- Algorithmes d'Assemblage ---
    int clipTips();
    int resolveBubbles();

    // RETOURNE MAINTENANT DES BITVECTORS
    std::vector<BitVector> generateContigs() const;

    void exportToGFA(const std::string& filename) const;

    // FUSIONNE VIA BITVECTOR
    static std::vector<BitVector> mergeContigs(std::vector<BitVector> contigs, int min_overlap);
};

#endif //ASSEMBLEUR_GRAPHDBJ_H