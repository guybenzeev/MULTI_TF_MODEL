#pragma once

#include <vector>
#include <memory>

#include "Edit.h"
#include "State.h"
#include "SubChainTopo.h"

/**
 * @file SubChain.h
 * @brief Declares the local state associated with one strand site.
 */

/**
 * @class SubChain
 * @brief Represents one site in the larger DNA-strand Markov chain.
 *
 * Each sub-chain stores its local binding state. Candidate edits and their
 * rates are supplied by the shared SubChainTopo.
 */
class SubChain {
private:
    /** Topology shared with the owning chain. */
    const SubChainTopo& topo_;

    /** Zero-based position within the strand. */
    const int nodeId_;

    /** Current local binding state. */
    State currentState;

public:
    /**
     * @brief Constructs a free site.
     * @param topo Topology that defines edits and rates.
     * @param nodeId Zero-based position in the strand.
     */
    SubChain(const SubChainTopo& topo, int nodeId)
        : topo_(topo),
          nodeId_(nodeId),
          currentState()
    {}

    /** @brief Destroys the site. */
    virtual ~SubChain() = default;

    /**
     * @brief Selects an edit from this site's cumulative rate intervals.
     *
     * @param threshold Local threshold in `[0, getTotalExitRate(chain))`.
     * @param chain Complete strand configuration.
     * @return Selected edit.
     * @throws std::runtime_error When no rate interval contains @p threshold.
     */
    virtual Edit findNextEdit(double threshold, const std::vector<std::unique_ptr<SubChain>>& chain);

    /**
     * @brief Sums the currently available edit rates.
     * @param chain Complete strand configuration.
     * @return Total exit rate from this site.
     */
    virtual double getTotalExitRate(
        const std::vector<std::unique_ptr<SubChain>>& chain
    );

    /**
     * @brief Applies an edit and updates footprint occupancy.
     * @param edit Edit to apply.
     * @param proteinHeadIndex Zero-based head site of the affected protein.
     */
    virtual void applyEdit(Edit& edit, int proteinHeadIndex);

    /**
     * @brief Returns candidate edits for the current binding state.
     * @return Immutable topology-owned edit list.
     */
    const SubChainTopo::EditList& getPossibleEditsForCurrentState() const {
        return topo_.getPossibleEditsForState(currentState);
    }

    /**
     * @brief Returns the local binding state.
     * @return Immutable reference to the current state.
     */
    const State& getCurrentState() const { return currentState; }
};
