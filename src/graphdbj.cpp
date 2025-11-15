#include "graphdbj.h"

#include <algorithm>  // pour std::sort
#include <stdexcept> // pour std::runtime_error
#include <iostream>	// logs optionnels
using namespace std;
/* Début Mickael – Implémentation 10/11/2025
 *
 * Dans ce fichier, comme  discuté en équipe : la construction du graphe de De Bruijn..
 * La philosophie de cette implémentation est de rester fidèle à l’algorithme du cours (double boucle i/j sur les séquences). Découpler autant que possible les responsabilités.
 */

 // CONSTRUCTEUR

GraphDBJ::GraphDBJ(size_t k): k(k)
{
    if (k < 2)
        throw runtime_error("GraphDBJ: k doit être >= 2 pour construire un graphe de De Bruijn.");
}

// CONSTRUCTION À PARTIR D’UNE LISTE DE SEQUENCES

/**
 * @brief
 * On applique directement le pseudocode vu en cours :
 *
 * Début
 *    pour i = 0 .. n-1
 *       pour j = 0 .. |Fi| - k
 *          inserer( L , sous_chaine(Fi, j, k) )
 *       fin pour
 *    fin pour
 * Fin
 *
 * Où Fi est la i-ème séquence.
 *
 * Ici, nous déléguons la génération des k-mers à process_sequence(),
 * ce qui nous donne une meilleure lisibilité et une cohérence interne.
 *
 * @complexity
 * O(N * L) k-mers générés, chaque insertion coût amorti O(1).
 */
void GraphDBJ::build_from_sequences(const vector<string>& sequences)
{
    adjacency.clear();  // Par sécurité : permet de reconstruire proprement

    for (const std::string& seq : sequences)
    {
        if (seq.size() >= k)
            process_sequence(seq);
        // en-dessous, la séquence est trop courte : nous avons décidé  simplement de l'ignorer (pas de k-mer possible)
    }
}


// EXTRACTION DES KMERS

/**
 * @brief
 * Parcourt la séquence brute et génère tous les k-mers possibles.
 * Chaque k-mer est ensuite passé à insert_kmer(), qui gère le découpage
 * prefixe/suffixe ((k-1)-mers).
 *
 * Cette fonction traduit directement la logique du j = 0 .. |Fi| - k.
 */
void GraphDBJ::process_sequence(const string& seq)
{
    const size_t limit = seq.size() - k + 1;

    for (size_t j = 0; j < limit; ++j)
    {
        string kmer = seq.substr(j, k);
        insert_kmer(kmer);
    }
}


//  INSERTION PREFIX/SUFFIXE
/**
 * @brief
 * Concrétise la structure d’un graphe de De Bruijn :
 * 
 *    k-mer = BANANE  (k=6)
 *    prefix = BANAN  (k-1)
 *    suffix = ANANE  (k-1)
 *
 * Nous stockons prefix -> {suffix1, suffix2, ...}
 * Pourquoi "unordered_set" ? Nous avons discuté de l’usage d’une liste versus set. Le set garantit l’unicité des edges sans, à priori, de tri trop coûteux. 
 * Cela limite la duplication inutile et nous simplifie grandement la phase de tri final.
 */
void GraphDBJ::insert_kmer(const string& kmer)
{
    if (kmer.size() != k)
        return; // garde fou

    const string prefix = kmer.substr(0, k - 1);
    const string suffix = kmer.substr(1, k - 1);

    adjacency[prefix].insert(suffix);
}

// RECUPERATION DES NOEUDS TRIES

/**
 * @brief
 * Nous avons souvent eu besoin dans certaines étapes en aval d’obtenir une liste
 * canonique et déterministe des (k-1)-mers. Le tri permet de la  reproductibilité, iteration cohérente entre runs,  alignement avec CompareKMers si des analyses croisées sont réalisées.
 * Le tri est effectué sur une copie (pas de mutation interne).
 * @complexity
 * O(V log V) où V = nombre de noeuds (k-1)-mers distincts.
 */
vector<std::string> GraphDBJ::get_sorted_nodes() const
{
    vector<string> nodes;
    nodes.reserve(adjacency.size());

    for (const auto& entry : adjacency)
    {
        nodes.push_back(entry.first);
    }

    sort(nodes.begin(), nodes.end());
    return nodes;
}

// ACCÈS DIRECT À LA STRUCTURE D’ADJACENCE

const unordered_map<string, unordered_set<string>>& 
GraphDBJ::get_graph() const
{
    return adjacency;
}
