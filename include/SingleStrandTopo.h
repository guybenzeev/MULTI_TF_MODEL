#pragma once
#include <unordered_map>
#include <memory>
#include <utility>
#include <vector>

#include "SubChainTopo.h"

/**
 * @file SingleStrandTopo.h
 * @brief Declares the one-sided DNA strand topology.
 */

/**
 * @class SingleStrandTopo
 * @brief Topology for a strand with one binding side.
 *
 * The topology builds candidate edits for free, non-specifically bound, and
 * specifically bound states. Protein models determine the corresponding rates.
 */
class SingleStrandTopo : public SubChainTopo {
private: 
    /** @brief Builds the state-to-candidate-edit dictionary. */
    void buildPossibleEditsDict();

public:
    /**
     * @brief Constructs a single-sided topology.
     * @param proteins Protein models transferred to the topology.
     */
    explicit SingleStrandTopo(std::vector<std::unique_ptr<Protein>> proteins = {})
        : SubChainTopo(1, std::move(proteins))
    {
        buildPossibleEditsDict();
    }

protected:
    /** @copydoc SubChainTopo::computeRate */
    double computeRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId

    ) const override;

    /** @return Number of local states supported by this topology. */
    int getNumStates() const override {
        return (2 * num_sides_ * getNumProteins() + 1);
    }
};
