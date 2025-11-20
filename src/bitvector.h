#pragma once

#include <vector>
#include <bitset>
#include <string>
#include <utility>
#include <cstddef>
#include <cstdint>

class BitVector
{
public:
    static constexpr size_t BLOCK_BITS = 64;

protected:
    std::vector<std::bitset<BLOCK_BITS>> _bitVector;

    size_t _sizeElement = 2;
    size_t _bitCount = 0;
    std::vector<std::pair<char, int>> _listBit;

public :
    BitVector();
    BitVector(int sizeElement);
    BitVector(int sizeElement, size_t dataSize);

    void reserve(size_t bits);
    void push_back(bool bit);
    size_t size() const;
    void clear();
    std::vector<bool> to_vector() const;

    void setSizeElement(int sizeElement);
    int getSizeElement() const;

    void createListBit(const std::string& list);
    const std::vector<std::pair<char, int>>& getListBit() const;

    void addCha(char cha);
    std::string readBitVector() const;

    bool test(size_t idx) const;
    bool operator[](size_t idx) const { return test(idx); }

    /// TRUCS QUE J'AI RAJOUTÉ
    /**
     * @brief Ajoute le contenu d'un autre BitVector à la fin de celui-ci.
     * @param other Le vecteur à ajouter.
     * @param skip_nucleotides Nombre de nucléotides à ignorer au début de 'other' (pour la fusion).
     */
    void append(const BitVector& other, size_t skip_nucleotides = 0);

    /**
     * @brief Récupère la valeur (0-3) du nucléotide à l'index donné.
     */
    uint8_t getNucleotideAt(size_t index) const;

    void resize(size_t num_bits);
};