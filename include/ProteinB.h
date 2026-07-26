#pragma once

#include <memory>
#include <vector>

#include "Protein.h"

/**
 * @file ProteinB.h
 * @brief Declares the second built-in protein kinetic model.
 */

/**
 * @class ProteinB
 * @brief Five-site protein with the Protein B transition-rate constants.
 */
class ProteinB : public Protein {
private: 
    const double sliding_rate_ = 1.0;    /**< Left/right sliding rate. */
    const double switch_side_rate_ = 0.1; /**< Strand-side switching rate. */
    const double bind_ns_rate_ = 1.0;    /**< Non-specific binding rate. */
    const double bind_s_rate_ = 0.1;     /**< Specific binding rate. */
    const double unbind_ns_rate_ = 0.1;  /**< Non-specific unbinding rate. */
    const double unbind_s_rate_ = 0.05;  /**< Specific unbinding rate. */

protected:
    /**
     * @brief Computes a sliding transition rate.
     * @param currentState State at the protein head.
     * @param edit Proposed left or right slide.
     * @param chain Complete strand configuration.
     * @param nodeId Zero-based protein-head index.
     * @return Sliding rate, or zero when movement is blocked.
     */
    double slideRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        int nodeId
    ) const;


public:
    /** @brief Constructs Protein B with radius 1 and width 5. */
    ProteinB()
        : Protein(1, 5)
    {}

    /** @copydoc Protein::computeRate */
    double computeRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        int nodeId
    ) const override;
};
