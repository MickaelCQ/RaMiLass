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

    /**
     *  @brief  save memory for n bits.
     *  @complexity : space : O(n), time : O(b)
     *  n = nb de bits. b =  nb de block.
     * **/
    void reserve(size_t bits);

    /**
     *  @brief  add bit at the end.
     *  @complexity : space : O(1), time : O(b).
     *  b = nb de block. (pire cas, sinon O(1))
     * **/
    void push_back(bool bit);

    /**
     *  @brief  return quantity of bits.
     *  @complexity : space : O(1), time : O(1).
     * **/
    size_t size() const;


    /**
     *  @brief  clear bitvector.
     *  @complexity : space : O(1), time : O(b+l)
     *  b = nb de block. l = taille de alphabet bit.
     * **/
    void clear();

    /**
     *  @brief  convert bitvector to bool vector.
     *  @complexity : space : O(n), time : O(n)
     *  n = nb de bits.
     * **/
    std::vector<bool> to_vector() const;

    /**
     *  @brief  set nb of bits per caracter.
     *  @complexity : space : O(1), time : O(1).
     * **/
    void setSizeElement(int sizeElement);

    /**
     *  @brief  get nb of bits per caracter.
     *  @complexity : space : O(1), time : O(1).
     * **/
    int getSizeElement() const;

    /**
     *  @brief  create list bit from string.
     *  @complexity : space : O(m), time : O(m)
     *  m = taille de la chaine de caractére.
     * **/
    void createListBit(const std::string& list);

    /**
     *  @brief  get list bit.
     *  @complexity : space : O(1), time : O(1).
     * **/
    const std::vector<std::pair<char, int>>& getListBit() const;


   /**
    *  @brief  add caracter to bitvector.
    *  @complexity : space : O(b), time : O(b).
    *  b = nb de block. (pire cas, sinon O(1) en temps et espace)
    * **/
    void addCha(char cha);

    /**
     * @brief from the bits, show the sequence on a string.
     * @complexity : space : O(n/2), time : O(n/2)
     * n = nb de bits.
     * **/
    std::string readBitVector() const;

    /** @brief Test bit at given index.
     * @complexity : space : O(1), time : O(1).
     * **/
    bool test(size_t idx) const;

    /** @brief Operator to access bit at given index.
     * @complexity : space : O(1), time : O(1).
     * **/
    bool operator[](size_t idx) const { return test(idx); }

    /// TRUCS QUE J'AI RAJOUTÉ
    /**
     * @brief Ajoute le contenu d'un autre BitVector à la fin de celui-ci.
     * @param other Le vecteur à ajouter.
     * @param skip_nucleotides Nombre de nucléotides à ignorer au début de 'other' (pour la fusion).
     * @complexity : space : O(m + b), time : O(m)
     * m = nombre de nouveau bit. b = nb de block (pire cas, sinon O(1) en espace).
     * **/
    void append(const BitVector& other, size_t skip_nucleotides = 0);

    /**
     * @brief Récupère la valeur (0-3) du nucléotide à l'index donné.
     * @complexity : space : O(1), time : O(1).
     * **/
    uint8_t getNucleotideAt(size_t index) const;

    /**
     * @brief Redimensionne bitVector.
     * @complexity : space : O(1), time : O(1)
     * **/
    void resize(size_t num_bits);
};