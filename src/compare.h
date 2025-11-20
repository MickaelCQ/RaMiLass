#ifndef ASSEMBLEUR_COMPARE_H
#define ASSEMBLEUR_COMPARE_H
#include <string>
#include <vector>
#include <cstddef> // Pour size_t
#include "bitvector.h"

/**
 * @struct CompareKMers
 * @brief Compare les k-mers à partir d'un vecteur binaire de nucléotides.
 *
 * Cette classe permet d'analyser les similarités entre lectures sans construire 
 * explicitement le graphe, en travaillant directement sur le BitVector compressé.
 */
struct CompareKMers {
    BitVector bitVector;            // Copie ou référence au vecteur binaire
    std::vector<size_t> reads;      // Vecteur des positions de fin de lecture (en bits)
    size_t kmersize = 31;

    CompareKMers(BitVector bitVector, std::vector<size_t> reads, size_t kmersize);
    CompareKMers(BitVector bitVector, std::vector<size_t> reads);

    // Setters & Getters
    /** @brief
     *  @complexity : space : O(1), time : O(1). pour tous les getters et setters.
     * */
    void set_kmersize(size_t k);
    size_t get_kmersize() const;
    BitVector& get_bitVector();
    std::vector<size_t>& get_reads();
    
    /**
     * @brief Récupère la position de fin (en bits) d'une lecture donnée.
     * @complexity : space : O(1), time : O(1).
     */
    size_t get_read_end_pos(size_t read_idx) const;
    
    /**
     * @brief Calcule le nombre théorique de k-mers dans une lecture donnée.
     * @complexity : space : O(1), time : O(1).
     */
    size_t get_nKmers(size_t ref_read_idx) const;
    /**
     * @brief Calcule le nombre total de k-mers dans toutes les lectures.
     * @complexity : space : O(1), time : O(r)
     * r  = nb de lectures.
     */
    size_t get_all_nKmers() const;
    /**
     * @brief Récupère le nombre total de lectures.
     * @complexity : space : O(1), time : O(1).
     */
    size_t get_nReads() const;

    // Méthodes de comparaison
    
    /**
     * @brief Compare le chevauchement (k-1) de deux k-mers bit à bit.
     * @param bit_idx1 Index binaire de début du premier k-mer.
     * @param bit_idx2 Index binaire de début du second k-mer.
     * @return 1 si le suffixe (k-1) de k-mer1 correspond au préfixe (k-1) de k-mer2, sinon 0.*
     * @complexity : space : O(1), time : O(k)
     * k = taille du kmer.
     */
    size_t compare_line(size_t bit_idx1, size_t bit_idx2) const;

    /**
     * @brief Compare le premier k-mer d'une lecture de référence à tous les autres.
     * @param ref_read_idx Index de la lecture de référence.
     * @return Vecteur de résultats (1 pour match, 0 sinon).
     * @complexity : space : O(r), time : O(r*k)
     * r = nb de lectures, k = taille du kmer.
     */
    std::vector<size_t> compare_lines(size_t ref_read_idx) const;

    /**
     * @brief Génère une matrice complète de comparaison (toutes les lectures contre toutes).
     * @return Matrice 2D des résultats.
     * @complexity : space : O(r^2), time : O(r^2*k)
     * r = nb de lectures, k = taille du kmer.
     */
    std::vector<std::vector<size_t>> compare_all() const;
};

#endif //ASSEMBLEUR_COMPARE_H