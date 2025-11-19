#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <chrono>

#include "graphdbj.h"
#include "convert.h"
#include "compare.h"

// --- Fonction utilitaire pour mesurer le temps ---
// Utilisation de high_resolution_clock pour la précision
template<typename TimeT = std::chrono::milliseconds>
struct Timer
{
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::string name;

    Timer(const std::string& taskName) : name(taskName) {
        start = std::chrono::high_resolution_clock::now();
        std::cout << "\n--- DEBUT : " << name << " ---" << std::endl;
    }

    ~Timer() {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<TimeT>(end - start);
        std::cout << "--- FIN : " << name << " en " << duration.count() << "ms ---" << std::endl;
    }
};
// -----------------------------------------------------------

int main(int argc, char* argv[]) {
    Timer<std::chrono::milliseconds> total_timer("EXECUTION TOTALE");

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <fichier_fasta.fa> <fichier_sortie_prefixe>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    std::string output_filename = argv[2];
    const size_t KMER_SIZE = 31;

    Convert converter;

    // 1. ÉTAPE : Traitement du fichier FASTA
    try {
        Timer<std::chrono::milliseconds> convert_timer("Conversion FASTA (2 Passes)");
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

    // Stats
    CompareKMers comparator(bitVector_ref, read_ends, KMER_SIZE);
    std::cout << "Lectures: " << comparator.get_nReads() << ", K-mers: " << comparator.get_all_nKmers() << std::endl;

    GraphDBJConfig config((int)KMER_SIZE);
    GraphDBJ graph(converter, (int)KMER_SIZE, config);

    auto nodes = graph.getNodes();
    std::cout << "Noeuds initiaux : " << nodes.size() << std::endl;

    // 2. ÉTAPE : Simplification du Graphe (Bubbles & Tips)
    {
        Timer<std::chrono::milliseconds> simplification_timer("Simplification du Graphe (Bubbles & Tips)");
        std::cout << "--- Simplification ---" << std::endl;
        bool changed = true;
        int iter = 0;
        while (changed && iter < config.MAX_PASSES) {
            changed = false; iter++;
            std::cout << "Passage " << iter << "..." << std::endl;
            if (graph.resolveBubbles() > 0) changed = true;
            if (graph.clipTips() > 0) changed = true;
            if (!changed) std::cout << "Simplification completee en " << iter << " passages." << std::endl;
        }
    }

    // 3. ÉTAPE : Génération des Contigs
    std::vector<BitVector> contigs;
    {
        Timer<std::chrono::milliseconds> generation_timer("Generation des Contigs (Chaines)");
        std::cout << "Generation des Contigs (Format BitVector)..." << std::endl;
        contigs = graph.generateContigs();
        std::cout << "Nombre de contigs bruts : " << contigs.size() << std::endl;
    }

    // 4. ÉTAPE : Fusion des Contigs
    {
        Timer<std::chrono::milliseconds> merge_timer("Fusion des Contigs (Overlap/Contained)");
        int min_overlap = (int)KMER_SIZE / 2;
        std::cout << "Fusion (BitVector logic) overlap min : " << min_overlap << std::endl;

        contigs = GraphDBJ::mergeContigs(
            contigs,
            min_overlap,
            config.ERROR_PERCENT_OVERLAP,
            config.ERROR_PERCENT_CONTAINED
        );

        std::cout << "Nombre de contigs finaux : " << contigs.size() << std::endl;
    }

    // 5. ÉTAPE : Export des fichiers
    {
        Timer<std::chrono::milliseconds> export_timer("Export GFA et FASTA");

        std::cout << "Export GFA..." << std::endl;
        graph.exportToGFA(output_filename + ".gfa");

        std::string contig_filename = output_filename + ".contigs.fasta";
        std::ofstream out_contigs(contig_filename);

        int contigs_exported = 0;
        for (size_t i = 0; i < contigs.size(); ++i) {
            size_t len = contigs[i].size() / 2;
            if (len >= (size_t)KMER_SIZE * 2) {
                out_contigs << ">contig_" << i << "_len_" << len << "\n";
                out_contigs << contigs[i].readBitVector() << "\n";
                contigs_exported++;
            }
        }
        out_contigs.close();
        std::cout << "Contigs exportes : " << contigs_exported << std::endl;
    }

    return 0;
}