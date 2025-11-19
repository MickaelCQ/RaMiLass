#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <bitset>

#include "graphdbj.h"
#include "convert.h"
#include "compare.h"

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <fichier_fasta.fa> <fichier_sortie_prefixe>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    std::string output_filename = argv[2];
    const size_t KMER_SIZE = 31;

    Convert converter;
    try {
        std::cout << "Traitement du fichier FASTA: " << filename << "..." << std::endl;
        converter.processFile(filename);
    } catch (const std::runtime_error& e) {
        std::cerr << "Erreur: " << e.what() << std::endl;
        return 1;
    }

    const BitVector& bitVector_ref = converter.get_bitVector();
    std::vector<size_t> read_ends = converter.get_read_end_positions();

    if (read_ends.empty()) {
        std::cout << "Aucune lecture trouvee." << std::endl;
        return 0;
    }

    // Stats (inchangé)
    CompareKMers comparator(bitVector_ref, read_ends, KMER_SIZE);
    std::cout << "Lectures: " << comparator.get_nReads() << ", K-mers: " << comparator.get_all_nKmers() << std::endl;

    std::cout << "Construction du graphe de Bruijn..." << std::endl;
    // Rechargement nécessaire dû au std::move dans CompareKMers précédemment
    converter.processFile(filename);
    GraphDBJ graph(converter, KMER_SIZE);

    auto nodes = graph.getNodes();
    std::cout << "Noeuds initiaux : " << nodes.size() << std::endl;

    std::cout << "--- Simplification ---" << std::endl;
    bool changed = true;
    int iter = 0;
    while (changed && iter < 10) {
        changed = false; iter++;
        if (graph.resolveBubbles() > 0) changed = true;
        if (graph.clipTips() > 0) changed = true;
    }

    std::cout << "Export GFA..." << std::endl;
    graph.exportToGFA(output_filename + ".gfa");

    std::cout << "Generation des Contigs (Format BitVector)..." << std::endl;

    // CHANGEMENT : Le type reçu est maintenant vector<BitVector>
    std::vector<BitVector> contigs = graph.generateContigs();
    std::cout << "Nombre de contigs bruts : " << contigs.size() << std::endl;

    int min_overlap = 15;
    std::cout << "Fusion (BitVector logic) overlap min : " << min_overlap << std::endl;
    contigs = GraphDBJ::mergeContigs(contigs, min_overlap);

    std::cout << "Nombre de contigs finaux : " << contigs.size() << std::endl;

    // ÉCRITURE FINALE : Conversion BitVector -> String uniquement ici
    std::string contig_filename = output_filename + ".contigs.fasta";
    std::ofstream out_contigs(contig_filename);

    int contigs_exported = 0;
    for (size_t i = 0; i < contigs.size(); ++i) {
        size_t len = contigs[i].size() / 2; // Taille en nucléotides
        if (len >= (size_t)KMER_SIZE * 2) {
            out_contigs << ">contig_" << i << "_len_" << len << "\n";
            // Utilise la méthode existante readBitVector() pour obtenir la string
            out_contigs << contigs[i].readBitVector() << "\n";
            contigs_exported++;
        }
    }
    out_contigs.close();
    std::cout << "Contigs exportes : " << contigs_exported << std::endl;

    return 0;
}