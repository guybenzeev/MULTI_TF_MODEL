#pragma once

#include <vector>
#include <memory>

#include "State.h"
#include "Edit.h"

class SubChain;

/**
 * @file Protein.h
 * @brief Declares the interface implemented by protein rate models.
 */

/**
 * @class Protein
 * @brief Base class for protein-specific geometry and transition rates.
 *
 * A protein model supplies its footprint, dependency radius, and the rate of
 * every edit that the topology may offer.
 */
class Protein {
protected:
    int transition_dependency_radius_; /**< Sites on either side that can affect a rate. */
    int width_; /**< Number of consecutive sites occupied by the protein. */

    /**
     * @brief Tests whether the protein footprint beginning at a site is free.
     * @param chain Complete strand configuration.
     * @param nodeId Zero-based head site of the proposed binding.
     * @return `true` when every required site exists and is free.
     */
    virtual bool checkForSpaceToBind(
        const std::vector<std::unique_ptr<SubChain>>& chain,
        int nodeId
    ) const;

public:
    /**
     * @brief Constructs a protein description.
     * @param transition_dependency_radius Number of neighboring sites that can affect a rate.
     * @param width Number of consecutive sites occupied when bound.
     */
    Protein(int transition_dependency_radius = 0, int width = 1)
        : transition_dependency_radius_(transition_dependency_radius),
          width_(width)
    {}

    /** @brief Enables safe destruction through a base pointer. */
    virtual ~Protein() = default;

    /**
     * @brief Computes the rate of a proposed edit.
     * @param currentState State at the proposed edit's head site.
     * @param edit Proposed transition.
     * @param chain Complete strand configuration.
     * @param nodeId Zero-based head-site index.
     * @return Transition rate; zero when the edit is not currently possible.
     */
    virtual double computeRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const = 0;

    /**
     * @brief Returns the transition dependency radius.
     * @return Number of neighboring sites that can affect this protein's rates.
     */
    virtual int getEffectRadius() const {
        return transition_dependency_radius_;
    }

    /**
     * @brief Returns the protein footprint.
     * @return Number of consecutive occupied sites.
     */
    virtual int getWidth() const {
        return width_;
    }
};
