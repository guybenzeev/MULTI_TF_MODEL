#pragma once

#include <vector>
#include <memory>
#include <functional>
#include <random>
#include <string>


#include "Edit.h"
#include "State.h"
#include "SubChain.h"
#include "SubChainTopo.h"

class SubChain;  ///< Forward declaration of SubChain.

/**
 * @class MultiTFMarkovChain
 * @brief Multi–transcription-factor continuous-time Markov chain on a DNA strand.
 *
 * This class models a Strand of DNA with multiple transcription factors (TFs).
 * At any time there is a current configuration (`currentState`), a history of 
 * previously visited configurations (`pastStates`)
 * and their visit times (`pastTimes`). Transitions between configurations are
 * driven by exit rates associated with each sub-chain and simulated using a
 * Gillespie-style algorithm.
 */
class MultiTFMarkovChain {

public:

    struct PastEvent {
        int index;
        Edit edit;
    };

    /**
    * @brief Update exit rates for all sub-chains after a state change.
    * This refreshes \c subChainExitRates and \c totalExitRate
    * to reflect the current configuration.
    */  
    void updateExitRates(int editHeadIndex, int proteinWidth);

    /**
     * @brief Select which sub-chain will fire next.
     *
     * Given a threshold in the range \f$[0, \text{totalExitRate})\f$,
     * this method walks through \c subChainExitRates and returns the
     * first sub-chain whose cumulative rate exceeds the threshold.
     * The returned pointer is non-owning; ownership is still managed
     * by the \c MultiTFMarkovChain.
     *
     * @param threshold A uniform random value in \f$[0,\text{totalExitRate})\f$
     *                  used to select the next sub-chain.
     * @return Index of the selected \c SubChain.
     */
    int selectSubChain(double& threshold);

    /**
    * @brief Construct a new MultiTFMarkovChain with a default initial state [all states free].
    *
    * @param length      Length of the DNA strand being modeled.
    * @param topo        Topology defining sub-chain structure and rates.
    */
    MultiTFMarkovChain(int length, const SubChainTopo& topo);

    /**
    * @brief Destroy the MultiTFMarkovChain.
    *
    * Frees all owned \c SubChain instances and any associated resources.
    */
    ~MultiTFMarkovChain();

    /**
    * @brief Run the CTMC simulation up to a runTime.
    *
    * Repeatedly calls \c updateChain() to generate events until the
    * internal time reaches or exceeds \p maxTime, or until no further
    * transitions are possible (e.g., \c totalExitRate becomes zero).
    *
    * @param  runTime simulation time to run to.
    */
    void runSimulation(double runTime);

    // DEBUGING METHODS 
    void printCurrentStates() const;
    void printExitRates();
    double getTotalExitRate() const;
    const std::vector<double>& getSubChainExitRates() const;
    void debugSetSelectedSubChainIndex(int idx) { selectedSubChainIndex_ = idx; }
    void debugUpdateExitRates(int centerIndex, int proteinWidth) { updateExitRates(centerIndex, proteinWidth); }
    std::vector<std::unique_ptr<SubChain>>& debugChain() { return currentChain_; }

    double getSpecificExitRate(int index){ 
        return currentChain_[index]->getTotalExitRate(currentChain_);
    }

    void debugApplyEditToSubChain(int index, Edit& edit){
        currentChain_[index]->applyEdit(edit);
    }

        /**
    * @brief Write the history of past events to a CSV file.
    *
    * Each row in the CSV corresponds to a past event, with columns for
    * the time of the event, the index of the sub-chain, and details of the edit.
    *
    * @param filename The path to the output CSV file.
    * @param runTime The total simulation time.
    */
    void writeToCSV(const std::string& filename, double runTime) const;

private:

    /// The topology defining sub-chain structure and rates.
    const SubChainTopo& topo_;

    /// Length of the DNA strand being modeled (in base pairs or sites).
    int length_;

    /// The current chain (ownership held by the MultiTFMarkovChain).
    std::vector<std::unique_ptr<SubChain>> currentChain_;

    /// The current simulation time.
    double currentTime_;

    /// The collection of past events visited, held as tuples of index and edit.
    std::vector<PastEvent> pastEvents_;

    /// The times at which past events were visited (aligned with \c pastEvents).
    std::vector<double> pastTimes_;

    /// Exit rates for each sub-chain (per-configuration total rates).
    std::vector<double> subChainExitRates_;

    /// Total exit rate across all sub-chains (sum of \c subChainExitRates).
    double totalExitRate_;

    /// Random number generator for stochastic sampling.
    std::mt19937 rng_;

    int selectedSubChainIndex_;


    // ===================== PRIVATE METHODS =====================

    /**
     * @brief Sample the holding time until the next event.
     *
     * Uses the current \c totalExitRate to draw the time increment to the
     * next state change, via an exponential distribution. The
     * returned value is in the same time units as the rates.
     *
     * @return The holding time until the next transition.
     */
    double findHoldingTime();



    /**
     * @brief Recompute possible edits and exit rates for all sub-chains.
     *
     * This updates the list of possible edits on each sub-chain, refreshes
     * \c subChainExitRates, and recomputes \c totalExitRate accordingly.
     *
     * @return An integer status code (e.g., number of active sub-chains or 0 on success).
     */
    int updatePossibleEdits();

    /**
     * @brief Apply the next event and update the global chain state.
     *
     * This function performs a single simulation step:
     *   - samples a holding time,
     *   - selects which sub-chain fires,
     *   - applies the corresponding edit,
     *   - records the new state in \c pastStates and \c pastTimes,
     *   - and updates all dependent rates.
     */
    void updateChain();

    /**
    * @brief Initialize exit rates for all sub-chains at the start of the simulation.
    *
    * This sets up \c subChainExitRates and \c totalExitRate based on the
    * initial configuration of the chain.
    */
    void initializeExitRates();

    void applyEdit(Edit& nextEdit);


};
