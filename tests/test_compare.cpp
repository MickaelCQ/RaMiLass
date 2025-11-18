#include "gtest/gtest.h"
#include "compare.h" // Inclut notre classe modifiée
#include <vector>
#include <stdexcept>

// Helper: convert std::vector<bool> literal to BitVector for easy test migration
static BitVector bv_from_vector(const std::vector<bool>& v) {
    BitVector bv;
    bv.reserve(v.size());
    for (bool b : v) bv.push_back(b);
    return bv;
}

// A=00, C=10, G=01, T=11

TEST(CompareKMersTest, GettersAndSetters) {
    // Read 1: "AC" -> {0,0, 1,0, 0,1}
    // Read 2: "GT" -> {0,1, 1,1}
    std::vector<bool> raw_bv = {false, false, true, false, false, true, true, true};
    std::vector<size_t> re = {4, 8}; // Fin de R1 à 4 bits, fin de R2 à 8 bits
    BitVector bb = bv_from_vector(raw_bv);
    CompareKMers comp(bb, re, 2); // k=2

    EXPECT_EQ(comp.get_kmersize(), 2);

    // CORRECTION: get_nReads
    EXPECT_EQ(comp.get_nReads(), 2);

    // CORRECTION: get_bitVector
    EXPECT_EQ(comp.get_bitVector().to_vector(), raw_bv);
    EXPECT_EQ(comp.get_reads(), re);

    comp.set_kmersize(3);
    EXPECT_EQ(comp.get_kmersize(), 3);
}

TEST(CompareKMersTest, ReadAndKmerCounting) {
    // Read 1: "ACG" -> {0,0, 1,0, 0,1}
    // Read 2: "T"   -> {1,1}
    // Read 3: "CGTA"-> {1,0, 0,1, 1,1, 0,0}
    std::vector<bool> raw_bv = {
        false, false, true, false, false, true, // ACG (6)
        true, true,                             // T (2)
        true, false, false, true, true, true, false, false // CGTA (8)
    };
    std::vector<size_t> re = {6, 8, 16}; // Positions de fin: 6, 8, 16
    BitVector bb = bv_from_vector(raw_bv);
    CompareKMers comp(bb, re, 3); // k=3

    // Test get_read_end_pos
    EXPECT_EQ(comp.get_read_end_pos(0), 6);
    EXPECT_EQ(comp.get_read_end_pos(1), 8);
    EXPECT_EQ(comp.get_read_end_pos(2), 16);
    EXPECT_THROW(comp.get_read_end_pos(3), std::out_of_range);

    // Test get_nKmers (k=3)
    // Read 1 ("ACG"): 3 nucs. (3 - 3 + 1) = 1 k-mer
    // CORRECTION: get_nKmers
    EXPECT_EQ(comp.get_nKmers(0), 1);
    // Read 2 ("T"): 1 nuc. (1 < 3) = 0 k-mers
    EXPECT_EQ(comp.get_nKmers(1), 0);
    // Read 3 ("CGTA"): 4 nucs. (4 - 3 + 1) = 2 k-mers
    EXPECT_EQ(comp.get_nKmers(2), 2);

    // Test get_all_nKmers
    // CORRECTION: get_all_nKmers
    EXPECT_EQ(comp.get_all_nKmers(), 1 + 0 + 2); // 3 k-mers au total
}

TEST(CompareKMersTest, CompareLine) {
    // k=3, donc chevauchement k-1 = 2
    // K-mer 1: "ACG" -> {0,0, 1,0, 0,1} (idx 0)
    // K-mer 2: "CGT" -> {1,0, 0,1, 1,1} (idx 6)
    // K-mer 3: "CGA" -> {1,0, 0,1, 0,0} (idx 12)
    std::vector<bool> raw_bv = {
        false, false, true, false, false, true, // ACG
        true, false, false, true, true, true,  // CGT
        true, false, false, true, false, false  // CGA
    };
    std::vector<size_t> re = {18}; // Une seule longue lecture
    BitVector bb = bv_from_vector(raw_bv);
    CompareKMers comp(bb, re, 3);

    // Compare "ACG" (idx 0) et "CGT" (idx 6)
    // Chevauchement: "CG" vs "CG". Doit matcher (1).
    EXPECT_EQ(comp.compare_line(0, 6), 1);

    // Compare "ACG" (idx 0) et "CGA" (idx 12)
    // Chevauchement: "CG" vs "CG". Doit matcher (1).
    EXPECT_EQ(comp.compare_line(0, 12), 1);

    // Compare "CGT" (idx 6) et "CGA" (idx 12)
    // Chevauchement: "GT" vs "CG". Ne doit pas matcher (0).
    EXPECT_EQ(comp.compare_line(6, 12), 0);
}

TEST(CompareKMersTest, CompareLinesAndAll) {
    // k=3
    // Read 0: "ACG" -> {0,0, 1,0, 0,1} (Début: 0)
    // Read 1: "CGT" -> {1,0, 0,1, 1,1} (Début: 6)
    // Read 2: "CGA" -> {1,0, 0,1, 0,0} (Début: 12)
    std::vector<bool> raw_bv = {
        false, false, true, false, false, true, // ACG
        true, false, false, true, true, true,  // CGT
        true, false, false, true, false, false  // CGA
    };
    std::vector<size_t> re = {6, 12, 18};
    BitVector bb = bv_from_vector(raw_bv);
    CompareKMers comp(bb, re, 3);

    // Test compare_lines(0) (Compare "ACG" à tous)
    std::vector<size_t> expected0 = {1, 1, 1};
    EXPECT_EQ(comp.compare_lines(0), expected0);

    // Test compare_lines(1) (Compare "CGT" à tous)
    std::vector<size_t> expected1 = {0, 1, 0};
    EXPECT_EQ(comp.compare_lines(1), expected1);

    // Test compare_all()
    std::vector<std::vector<size_t>> expected_all_final = {
        {1, 1, 1},
        {0, 1, 0},
        {0, 0, 1}
    };
    EXPECT_EQ(comp.compare_all(), expected_all_final);
}