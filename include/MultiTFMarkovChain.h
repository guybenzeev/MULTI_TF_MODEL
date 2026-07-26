#pragma once

#include <vector>
#include <memory>
#include <random>
#include <string>

#include "Edit.h"
#include "SubChain.h"
#include "SubChainTopo.h"
#include "FenwickTree.h"

/**
 * @file MultiTFMarkovChain.h
 * @brief Declares the multi-protein continuous-time Markov-chain simulator.
 */

/**
 * @class MultiTFMarkovChain
 * @brief Multi-protein continuous-time Markov chain on a DNA strand.
 *
 * This class models a strand of DNA with multiple proteins by running
 * trajectories through a continuous-time Markov chain (CTMC). It stores the
 * current configuration, the applied edits, and their event times. A
 * Gillespie-style algorithm selects events, while local rate updates and a
 * Fenwick tree avoid rescanning the full rate vector after every edit.
 *
 * The referenced topology must outlive the chain.
 */
class MultiTFMarkovChain {
public:
    /**
     * @brief Constructs a chain whose sites are initially free.
     * @param length Number of sites in the modeled DNA strand.
     * @param topo Topology defining candidate edits and transition rates.
     */
    MultiTFMarkovChain(int length, const SubChainTopo& topo);

    /** @brief Releases all owned sub-chains. */
    ~MultiTFMarkovChain();

    /**
     * @brief Runs a new stochastic trajectory up to a time limit.
     *
     * Existing event history and simulation time are reset before the run.
     *
     * @param runTime Simulation time limit.
     * @throws std::runtime_error If no transition can be selected.
     */
    void runSimulation(double runTime);

    /**
     * @brief Writes the current trajectory history to CSV.
     * @param filename Output path.
     * @param runTime Requested simulation time, recorded as CSV metadata.
     * @throws std::runtime_error If the output file cannot be opened.
     */
    void writeToCSV(const std::string& filename, double runTime) const;

private:
    /**
     * @struct PastEvent
     * @brief Associates an applied edit with its zero-based strand position.
     */
    struct PastEvent {
        int index; /**< Head-site index at which the edit fired. */
        Edit edit; /**< Applied transition. */
    };

    /** Non-owning topology reference. */
    const SubChainTopo& topo_;

    /** Number of sites in the modeled strand. */
    int length_;

    /** Owned local states in strand order. */
    std::vector<std::unique_ptr<SubChain>> currentChain_;

    /** Current trajectory time. */
    double currentTime_;

    /** Applied events in chronological order. */
    std::vector<PastEvent> pastEvents_;

    /** Event times aligned by index with pastEvents_. */
    std::vector<double> pastTimes_;

    /** Total exit rate of each sub-chain. */
    std::vector<double> subChainExitRates_;

    /** Sum of all values in subChainExitRates_. */
    double totalExitRate_;

    /** Random-number engine used for event and holding-time sampling. */
    std::mt19937 rng_;

    /** Zero-based sub-chain selected for the current event. */
    int selectedSubChainIndex_;

    /** Cumulative-rate index used for logarithmic event selection. */
    FenwickTree ratesTree_ = FenwickTree(0);

    /**
     * @brief Samples an exponentially distributed holding time.
     * @return Time until the next transition.
     * @throws std::runtime_error If the total exit rate is not positive.
     */
    double findHoldingTime();

    /**
     * @brief Selects a sub-chain from the cumulative rate tree.
     * @param threshold Global cumulative threshold; converted to a local threshold.
     * @return Zero-based selected index, or -1 when the threshold is out of range.
     */
    int selectSubChainFromTree(double& threshold);

    /**
     * @brief Initializes the per-site rates and cumulative-rate tree.
     */
    void initializeExitRates();

    /**
     * @brief Recomputes rates in the neighborhood affected by an edit.
     * @param editHeadIndex Zero-based head index after the edit.
     * @param proteinWidth Footprint of the edited protein.
     */
    void updateExitRates(int editHeadIndex, int proteinWidth);

    /**
     * @brief Applies an edit to every affected footprint site.
     * @param nextEdit Selected edit.
     */
    void applyEdit(Edit& nextEdit);
};
