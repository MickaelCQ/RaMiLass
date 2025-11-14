#ifndef GRAPH_DBJ_H
#define GRAPH_DBJ_H

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/* Début Mickael 10/11/2025 */

/**
 * @brief Classe responsable de la construction d'un graphe de k-mers, 
 * i.e. une version spécialisée du graphe de De Bruijn de taille k.
 *
 * Cette classe s’inscrit dans la continuité de la chaîne de traitement que nous avons 
 * conçue ensemble :
 * - Convert produit un vecteur de bits compact représentant les séquences.
 * - CompareKmers structure la comparaison entre lectures.
 * - GraphDBJ extrait explicitement les k-mers et construit une structure exploitable pour l’assemblage et l’analyse topologique.
 *
 * Nos choix, discutés en équipe et débattus:
 * - une structure simple (liste d’adjacence),
 * - un tri final des nœuds pour garantir la reproductibilité,
 * - un contrôle raisonnable des coûts mémoire/temps,
 * - une implémentation fidèle à l’algorithme abordé en cours.
 *
 * @complexity
 * Construction en O(N * L) où :
 * - N = nombre de séquences,
 * - L = longueur moyenne des séquences,
 * car chaque position j génère un k-mer et une insertion prefix→suffix.
 */
 
class GraphDBJ 
{
public:

    /**
     * @brief Constructeur principal.
     * @param k Taille des k-mers (typiquement 31).
     */
    GraphDBJ(size_t k);

    /**
     * @brief Construit le graphe à partir d’une liste de séquences déjà décodées.
     * @param sequences Liste de séquences brutes (ADN).
     */
    void build_from_sequences(const std::vector<std::string>& sequences);

    /**
     * @brief Retourne la liste triée des nœuds ((k-1)-mers) du graphe.
     * @return Un vecteur trié contenant l’ensemble des nœuds.
     */
    std::vector<std::string> get_sorted_nodes() const;

    /**
     * @brief Accès direct à la structure d’adjacence.
     * @return Une table associant chaque prefixe (k-1)-mer aux suffixes possibles.
     */
    const std::unordered_map<std::string, std::unordered_set<std::string>>& 
    get_graph() const;

private:

    size_t k; ///< Taille des k-mers.

    /** 
     * @brief Structure principale du graphe : 
     * chaque (k-1)-mer → ensemble des suffixes possibles.
     */
    std::unordered_map<std::string, std::unordered_set<std::string>> adjacency;

    /**
     * @brief Insère un k-mer en décomposant prefixe/suffixe ((k-1)-mers).
     * @param kmer Le k-mer complet.
     */
    void insert_kmer(const std::string& kmer);

    /**
     * @brief Extrait tous les k-mers d’une séquence donnée et les insère dans le graphe.
     * @param seq Séquence ADN brute.
     */
    void process_sequence(const std::string& seq);
};

#endif // GRAPH_DBJ_H
