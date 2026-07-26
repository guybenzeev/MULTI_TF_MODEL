#pragma once
#include <vector>
#include <memory>
#include <unordered_map>
#include <utility>

#include "State.h"
#include "Edit.h"
#include "Protein.h"

class SubChain;

/**
 * @file SubChainTopo.h
 * @brief Declares the topology interface shared by strand models.
 */

/**
 * @class SubChainTopo
 * @brief Defines valid edits, protein types, and rates for a strand topology.
 *
 * The topology owns its protein models. A MultiTFMarkovChain and its SubChain
 * objects retain a non-owning reference to the topology, so the topology must
 * outlive them.
 */
class SubChainTopo {
public:
    /** @brief Ordered list of edits available from one state. */
    using EditList = std::vector<Edit>;

    /** @brief Mapping from a binding state to its candidate edits. */
    using EditDict = std::unordered_map<State, EditList, StateHash>;

protected:
    int num_sides_; /**< Number of independently modeled strand sides. */
    std::vector<std::unique_ptr<Protein>> proteins_; /**< Owned protein rate models. */
    EditDict possibleEditsByState_; /**< Candidate edits indexed by binding state. */
    int largestProteinWidth_ = 1; /**< Maximum footprint across all protein types. */
    int largestRadius_ = 1; /**< Maximum dependency radius across all protein types. */

public:
    /**
     * @brief Constructs a topology and caches its maximum interaction extents.
     * @param num_sides_val Number of strand sides.
     * @param proteins Protein models transferred to this topology.
     */
    SubChainTopo(
        int num_sides_val,
        std::vector<std::unique_ptr<Protein>> proteins = {}
    )
        : num_sides_(num_sides_val),
          proteins_(std::move(proteins))
    {
        for (const auto& protein : proteins_) {
            if (protein->getWidth() > largestProteinWidth_) {
                largestProteinWidth_ = protein->getWidth();
            }
            if (protein->getEffectRadius() > largestRadius_) {
                largestRadius_ = protein->getEffectRadius();
            }
        }
    }

    /** @brief Enables safe destruction through a topology base pointer. */
    virtual ~SubChainTopo() = default;

    /** @return Number of modeled strand sides. */
    virtual int getNumSides() const {
        return num_sides_;
    }

    /** @return Number of available protein types. */
    virtual int getNumProteins() const {
        return static_cast<int>(proteins_.size());
    }

    /** @return Largest protein footprint in sites. */
    virtual int getLargestProteinWidth() const {
        return largestProteinWidth_;
    }

    /** @return Largest transition dependency radius in sites. */
    virtual int getLargestRadius() const {
        return largestRadius_;
    }

    /** @return Number of binding states supported at one site. */
    virtual int getNumStates() const = 0;

    /** @return Candidate-edit dictionary for the topology. */
    virtual const EditDict& getPossibleEditsByState() const {
        return possibleEditsByState_;
    }

    /**
     * @brief Looks up candidate edits for a state.
     * @param state Binding state to query.
     * @return Stored edit list, or an immutable empty list when absent.
     */
    const EditList& getPossibleEditsForState(const State& state) const {
        const auto& dict = getPossibleEditsByState();
        auto it = dict.find(state);
        static const EditList empty{};
        return (it == dict.end()) ? empty : it->second;
    }

    /**
     * @brief Computes the rate of an edit in the current chain.
     * @param currentState State at the edit's head site.
     * @param edit Candidate transition.
     * @param chain Complete strand configuration.
     * @param nodeId Zero-based head-site index.
     * @return Transition rate, or zero when the edit is unavailable.
     */
    virtual double computeRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const = 0;

    /**
     * @brief Returns a protein model by zero-based identifier.
     * @param proteinId Zero-based protein identifier.
     * @return Non-owning protein pointer, or `nullptr` when out of range.
     */
    virtual const Protein* getProteinForSubChain(int proteinId) const {
        if (proteinId < 0 || static_cast<std::size_t>(proteinId) >= proteins_.size()) {
            return nullptr;
        }
        return proteins_[static_cast<std::size_t>(proteinId)].get();
    }

    /**
     * @brief Returns a protein footprint by zero-based identifier.
     * @param proteinId Zero-based protein identifier.
     * @return Protein footprint, or one when the topology contains no proteins.
     * @pre @p proteinId is valid when the topology contains proteins.
     */
    virtual int getProteinWidth(int proteinId) const {
        if (proteins_.empty()) {
            return 1;
        }
        return proteins_[proteinId]->getWidth();
    }
};
