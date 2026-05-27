// src/Main.cpp
#include <iostream>
#include <string>

#include "MultiTFMarkovChain.h"
#include "SingleStrandTopo.h"

int main(int argc, char* argv[]) {
    try {
        // Defaults
        //int length = 100;
        int length = 5;
        double runTime = 1000000.0;
        //double runTime = 50.0;

        //std::string outFile = "data/events.csv";
        //std::string outFile = "data/mini_events.csv";
        std::string outFile = "data/reduced_events.csv";



        // Optional CLI args: ./sim [length] [runTime] [outFile]
        if (argc >= 2) length = std::stoi(argv[1]);
        if (argc >= 3) runTime = std::stod(argv[2]);
        if (argc >= 4) outFile = argv[3];

        // Build topology + chain
        SingleStrandTopo topo;                 // adjust if your topo ctor needs args
        MultiTFMarkovChain chain(length, topo);

        // Run simulation
        std::cout << "Running sim: length=" << length
                  << ", runTime=" << runTime
                  << ", outFile=" << outFile << "\n";

        chain.runSimulation(runTime);

        // Write results
        chain.writeToCSV(outFile, runTime);

        std::cout << "Done. Wrote CSV to: " << outFile << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        std::cerr << "Usage: ./sim [length] [runTime] [outFile]\n";
        return 1;
    }
}
