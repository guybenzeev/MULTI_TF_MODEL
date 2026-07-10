#pragma once

#include <memory>
#include <vector>

#include "Protein.h"

class ProteinA : public Protein {
private: 
    const double sliding_rate_ = 2.0;
    const double switch_side_rate_ = 0.1;
    const double bind_ns_rate_ = 1.0;
    const double bind_s_rate_ = 0.5;
    const double unbind_ns_rate_ = 0.1;
    const double unbind_s_rate_ = 0.05;

protected:
    double slideRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        int nodeId
    ) const;


public:

    ProteinA()
        : Protein(1, 1)
    {}

    double computeRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        int nodeId
    ) const override;
};