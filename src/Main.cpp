// src/Main.cpp
#include <iostream>
#include <string>
#include <chrono>
#include <filesystem>

#include "MultiTFMarkovChain.h"
#include "SingleStrandTopo.h"
#include "ProteinA.h"
#include "ProteinB.h"

int main(int argc, char* argv[]) {
    try {
        // Defaults
        int length = 50;
        double runTime = 50.0;

        //std::string outFile = "data/events.csv";
        //std::string outFile = "data/mini_events.csv";
        //std::string outFile = "data/reduced_events.csv";
        std::string outFile = "data/test_width_one_protein.csv";

        // Optional CLI args: ./sim [length] [runTime] [outFile]
        if (argc >= 2) length = std::stoi(argv[1]);
        if (argc >= 3) runTime = std::stod(argv[2]);
        if (argc >= 4) outFile = argv[3];

        // Build topology + chain
        std::vector<std::unique_ptr<Protein>> proteins;
        proteins.push_back(std::make_unique<ProteinA>());
        proteins.push_back(std::make_unique<ProteinB>());

        SingleStrandTopo topo(std::move(proteins));

        MultiTFMarkovChain chain(length, topo);

        // Run simulation
        std::cout << "Running sim: length=" << length
                  << ", runTime=" << runTime
                  << ", outFile=" << outFile << "\n";

        auto start = std::chrono::high_resolution_clock::now();

        chain.runSimulation(runTime);

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = end - start;

        std::cout << "Simulation runtime: "
                  << elapsed.count()
                  << " seconds\n";
        
        std::cout << "Current working directory: "
                  << std::filesystem::current_path() << "\n";

        std::cout << "Absolute CSV path: "
                  << std::filesystem::absolute(outFile) << "\n";

        // Write results
        chain.writeToCSV(outFile, runTime);

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        std::cerr << "Usage: ./sim [length] [runTime] [outFile]\n";
        return 1;
    }
}