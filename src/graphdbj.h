#ifndef ASSEMBLEUR_GRAPHDBJ_H
#define ASSEMBLEUR_GRAPHDBJ_H

#include <vector>
#include <cstdint> // Pour uint64_t
#include <unordered_map>
#include "convert.h"

/**
 * @struct Noeud
 * @brief Représente un (k-1)-mer dans le graphe de De Bruijn.
 */
struct Noeud {
	uint64_t p;             // L'entier représentant la séquence du (k-1)-mer
	std::vector<Noeud*> c;  // Vecteur de pointeurs vers les enfants (arêtes)

	// Constructeur pour faciliter la création
	Noeud(uint64_t val) : p(val) {}
};

/**
 * @class GraphDBJ
 * @brief Construit le graphe à partir des données compressées dans Convert.
 */
class GraphDBJ {
private:
	int k; // Taille du k-mer (les noeuds seront de taille k-1)

	// Hash map pour stocker l'unicité des noeuds :
	// Clé = valeur entière du (k-1)-mer, Valeur = Pointeur vers le Noeud unique
	std::unordered_map<uint64_t, Noeud*> nodes_map;

	/**
	 * @brief Extrait une valeur entière sur 64 bits représentant (len) nucléotides
	 * à partir d'une position donnée dans le BitVector.
	 */
	uint64_t extractKmerValue(const BitVector& bv, size_t start_bit_idx, int len_nucleotides) const;

public:
	/**
	 * @brief Constructeur. Construit immédiatement le graphe.
	 * @param converter L'objet Convert contenant les lectures chargées.
	 * @param kmer_size La taille k (les noeuds représentent k-1).
	 */
	GraphDBJ(const Convert& converter, int kmer_size);

	/**
	 * @brief Destructeur pour nettoyer la mémoire allouée aux noeuds.
	 */
	~GraphDBJ();

	/**
	 * @brief Récupère tous les noeuds du graphe.
	 * @return Un vecteur contenant les pointeurs vers tous les noeuds créés.
	 */
	std::vector<Noeud*> getNodes() const;
};

#endif //ASSEMBLEUR_GRAPHDBJ_H