#include "gtest/gtest.h"
#include "compare.h" // Inclut notre classe modifiée
#include <vector>
#include <stdexcept>

// A=00, C=10, G=01, T=11

TEST(CompareKMersTest, GettersAndSetters) {
    // Read 1: "AC" -> {0,0, 1,0}
    // Read 2: "GT" -> {0,1, 1,1}
    std::vector<bool> bv = {false, false, true, false, false, true, true, true};
    std::vector<size_t> re = {4, 8}; // Fin de R1 à 4 bits, fin de R2 à 8 bits
    CompareKMers comp(bv, re, 2); // k=2

    EXPECT_EQ(comp.get_kmersize(), 2);
    EXPECT_EQ(comp.get_n_reads(), 2);
    EXPECT_EQ(comp.get_bit_vector(), bv);
    EXPECT_EQ(comp.get_reads(), re);

    comp.set_kmersize(3);
    EXPECT_EQ(comp.get_kmersize(), 3);
}

TEST(CompareKMersTest, ReadAndKmerCounting) {
    // Read 1: "ACG" -> {0,0, 1,0, 0,1}
    // Read 2: "T"   -> {1,1}
    // Read 3: "CGTA"-> {1,0, 0,1, 1,1, 0,0}
    std::vector<bool> bv = {
        false, false, true, false, false, true, // ACG (6)
        true, true,                             // T (2)
        true, false, false, true, true, true, false, false // CGTA (8)
    };
    std::vector<size_t> re = {6, 8, 16}; // Positions de fin: 6, 8, 16
    CompareKMers comp(bv, re, 3); // k=3

    // Test get_read_end_pos
    EXPECT_EQ(comp.get_read_end_pos(0), 6);
    EXPECT_EQ(comp.get_read_end_pos(1), 8);
    EXPECT_EQ(comp.get_read_end_pos(2), 16);
    EXPECT_THROW(comp.get_read_end_pos(3), std::out_of_range);

    // Test get_n_kmers (k=3)
    // Read 1 ("ACG"): 3 nucs. (3 - 3 + 1) = 1 k-mer
    EXPECT_EQ(comp.get_n_kmers(0), 1);
    // Read 2 ("T"): 1 nuc. (1 < 3) = 0 k-mers
    EXPECT_EQ(comp.get_n_kmers(1), 0);
    // Read 3 ("CGTA"): 4 nucs. (4 - 3 + 1) = 2 k-mers
    EXPECT_EQ(comp.get_n_kmers(2), 2);

    // Test get_all_n_kmers
    EXPECT_EQ(comp.get_all_n_kmers(), 1 + 0 + 2); // 3 k-mers au total
}

TEST(CompareKMersTest, CompareLine) {
    // k=3, donc chevauchement k-1 = 2
    // K-mer 1: "ACG" -> {0,0, 1,0, 0,1} (idx 0)
    // K-mer 2: "CGT" -> {1,0, 0,1, 1,1} (idx 6)
    // K-mer 3: "CGA" -> {1,0, 0,1, 0,0} (idx 12)
    std::vector<bool> bv = {
        false, false, true, false, false, true, // ACG
        true, false, false, true, true, true,  // CGT
        true, false, false, true, false, false  // CGA
    };
    std::vector<size_t> re = {18}; // Une seule longue lecture
    CompareKMers comp(bv, re, 3);

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
    std::vector<bool> bv = {
        false, false, true, false, false, true, // ACG
        true, false, false, true, true, true,  // CGT
        true, false, false, true, false, false  // CGA
    };
    std::vector<size_t> re = {6, 12, 18};
    CompareKMers comp(bv, re, 3);

    // Test compare_lines(0) (Compare "ACG" à tous)
    // vs R0 ("ACG"): 1 (self)
    // vs R1 ("CGT"): compare_line(0, 6) -> 1 (match "CG")
    // vs R2 ("CGA"): compare_line(0, 12) -> 1 (match "CG")
    std::vector<size_t> expected0 = {1, 1, 1};
    EXPECT_EQ(comp.compare_lines(0), expected0);

    // Test compare_lines(1) (Compare "CGT" à tous)
    // vs R0 ("ACG"): compare_line(6, 0) -> 0 (match "GT" vs "AC")
    // vs R1 ("CGT"): 1 (self)
    // vs R2 ("CGA"): compare_line(6, 12) -> 0 (match "GT" vs "CG")
    std::vector<size_t> expected1 = {0, 1, 0};
    EXPECT_EQ(comp.compare_lines(1), expected1);

    // Test compare_all()
    std::vector<std::vector<size_t>> expected_all = {
        {1, 1, 1}, // Résultats pour R0
        {0, 1, 0}, // Résultats pour R1
        {0, 1, 1}  // Résultats pour R2 (compare_line(12, 0)=0, compare_line(12, 6)=1, self=1)
    };
    // Calcul manuel pour R2:
    // vs R0 ("ACG"): compare_line(12, 0) -> "GA" vs "AC" -> 0
    // vs R1 ("CGT"): compare_line(12, 6) -> "GA" vs "CG" -> 0
    // Oups, ma logique de test était fausse. Recalculons.

    // R0 vs R0: 1
    // R0 vs R1: compare_line(0, 6) ["ACG" vs "CGT"] -> chevauchement "CG" vs "CG" -> 1
    // R0 vs R2: compare_line(0, 12) ["ACG" vs "CGA"] -> chevauchement "CG" vs "CG" -> 1
    std::vector<size_t> r0_res = {1, 1, 1};
    EXPECT_EQ(comp.compare_lines(0), r0_res);

    // R1 vs R0: compare_line(6, 0) ["CGT" vs "ACG"] -> chevauchement "GT" vs "AC" -> 0
    // R1 vs R1: 1
    // R1 vs R2: compare_line(6, 12) ["CGT" vs "CGA"] -> chevauchement "GT" vs "CG" -> 0
    std::vector<size_t> r1_res = {0, 1, 0};
    EXPECT_EQ(comp.compare_lines(1), r1_res);

    // R2 vs R0: compare_line(12, 0) ["CGA" vs "ACG"] -> chevauchement "GA" vs "AC" -> 0
    // R2 vs R1: compare_line(12, 6) ["CGA" vs "CGT"] -> chevauchement "GA" vs "CG" -> 1 (G=G)
    // Attendez, k-1 = 2. "GA" vs "CG".
    // i=0: 'G' vs 'C'. Pas de match. Doit être 0.
    // i=1: 'A' vs 'G'.
    // Test compare_line(12, 6):
    // i=0: nuc1(idx 1)='G' (bits 14,15={0,0}?? NON, {0,1}). nuc2(idx 0)='C' (bits 6,7={1,0}). Match pas. -> 0
    // Ah, "CGA" -> {1,0, 0,1, 0,0}. bits 12,13='C', 14,15='G', 16,17='A'.
    // "CGT" -> {1,0, 0,1, 1,1}. bits 6,7='C', 8,9='G', 10,11='T'.
    // compare_line(12, 6):
    // i=0: nuc1(idx 1)='G' (bits 14,15) vs nuc2(idx 0)='C' (bits 6,7).
    //      {0,1} vs {1,0}. Pas de match. Retourne 0.
    EXPECT_EQ(comp.compare_line(12, 6), 0); // Le test ci-dessus était faux.

    // R2 vs R0: compare_line(12, 0) ["CGA" vs "ACG"] -> chevauchement "GA" vs "AC" -> 0
    // R2 vs R1: compare_line(12, 6) ["CGA" vs "CGT"] -> chevauchement "GA" vs "CG" -> 0
    // R2 vs R2: 1
    std::vector<size_t> r2_res = {0, 0, 1};
    EXPECT_EQ(comp.compare_lines(2), r2_res);

    // Test final compare_all()
    std::vector<std::vector<size_t>> expected_all_final = {
        {1, 1, 1},
        {0, 1, 0},
        {0, 0, 1}
    };
    EXPECT_EQ(comp.compare_all(), expected_all_final);
}