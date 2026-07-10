#pragma once
#include <vector>
#include <memory>
#include <stdexcept>
#include <string>
#include "Edit.h"
#include "State.h"
#include "SubChainTopo.h"


/**
 * @brief Class representing a single sub-chain representing an offset in the larger DNA strand.
 *
 * A SubChain models an offset in the larger continuous-time Markov process. Each
 * sub-chain tracks its current state and maintains internal/external possible
 * edits (transitions), each with an associated transition rate.
 *
 * Derived classes must implement:
 *  - findNextEdit() for Gillespie event selection
 *  - updatePossibleEdits() to refresh transition rates based on the global chain
 */


class SubChain {

    private:

    const SubChainTopo& topo_;

    /// Number of possible states in this sub-chain.
    const int numStates_;

    /// Unique identifier for this node/sub-chain within the global model.
    const int nodeId_;

    /// The current state of this sub-chain.
    State currentState;

    /// Current state ID of this sub-chain.
    int currentStateID;

    int proteinHead_;

    /// Possible transitions within the sub-chain.
    std::vector<Edit> possibleInternalEdits_;

    /// Possible transitions to slide out of the sub-chain.
    std::vector<Edit> possibleExternalEdits_;

    /**
     * @brief Computes the rate.
     *
     * This method calculates the transition rate from the current state to the "toState"
     * based on the provided topo and the states of other sub-chains in the global model.
     *
     * @param toState The state we may transition to.
     * @param Chain The larger DNA chain of SubChains.
     * @return The rate to toState.
     */
    virtual double computeInternalRate(int toState, const std::vector<std::unique_ptr<SubChain>>& Chain);
    
    /**
     * @brief Computes the rate to slide.
     *
     * This method calculates the transition rate from the current state to the "toState"
     * based on the provided topo and the states of other sub-chains in the global model.
     *
     * @param toState The state we may transition to.
     * @param Chain The larger DNA chain of SubChains.
     * @return The rate to toState.
     */
    virtual double computeSlidingRate(int toState, const std::vector<std::unique_ptr<SubChain>>& Chain);
    
    public:

    /**
    * @brief Constructor used by derived classes to initialize immutable members.
    *
    * @param effectRadius Radius of effect for TF changes.
    * @param numStates   Number of possible states in this sub-chain.
    * @param nodeId      Unique identifier for this sub-chain in the global model.
    * @param initialState Initial state of this sub-chain.
    */
    SubChain(const SubChainTopo& topo,
            int nodeId,
            const int initialState)
        : topo_(topo),
        numStates_(topo.getNumStates()),
        nodeId_(nodeId),
        proteinHead_(nodeId),
        currentStateID(initialState)
    {currentState = State();}


    /**
     * @brief Virtual destructor for polymorphic cleanup.
     */
    virtual ~SubChain() = default;

    /**
     * @brief Selects the next edit based on a random threshold.
     *
     * This method implements a subsbtep in the Gillespie "event selection" step, 
     * where the caller already selected this sub-chain to execute an edit.
     * The caller supplies a random threshold in \f$[0, R_{\text{total}})\f$,
     * and the edit whose cumulative rate interval contains the threshold
     * is chosen.
     *
     * @param threshold A random value in \f$[0, R_{\text{total}})\f$.
     * @return The State reached by executing the selected edit.
     */
    virtual Edit findNextEdit(double threshold, const std::vector<std::unique_ptr<SubChain>>& chain);

    /**
     * @brief Returns the total exit rate from the current state.
     *
     * This function gets the transition rates for every possible edit as listed in 
     * the edit table and sums them up.
     *
     * @param chain Vector of all sub-chain states in the global model.
     * @return The total exit rate \f$R_{\text{total}}\f$ from the current state.
     */
    virtual double getTotalExitRate(
        const std::vector<std::unique_ptr<SubChain>>& chain
    );

    /**
    * @brief Apply the given edit to this sub-chain.
    *
    * This method updates the internal state of the sub-chain
    * according to the specified edit.
    *
    * @param edit The edit to apply.
    */
    virtual void applyEdit(Edit& edit, int proteinHeadIndex);

    void applyEdit(Edit& edit) {
        applyEdit(edit, nodeId_);
    }
    /**
     * @brief Get the effect radius for this sub-chain.
     *
     * @return The effect radius.
     */
    int getRadius() const {
        if (currentState.protein == 0) {
            return 1;
        }
        const Protein* protein = topo_.getProteinForSubChain(currentState.protein-1);
        return protein->getEffectRadius();
    }

    const SubChainTopo::EditDict& getEditTable() const {
        return topo_.getPossibleEditsByState();
    }

    // Access edits for current state
    const SubChainTopo::EditList& getPossibleEditsForCurrentState() const {
        return topo_.getPossibleEditsForState(currentState);
    }

    int getCurrentStateID() const { return currentStateID; }
    const State& getCurrentState() const { return currentState; }
};
