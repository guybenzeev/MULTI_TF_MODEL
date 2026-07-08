#pragma once
#include <unordered_map>
#include <memory>
#include <utility>
#include <vector>

#include "SubChainTopo.h"

class SingleStrandTopo : public SubChainTopo {
private:
 //events and mini_events
 
    const double sliding_rate_ = 2.0;
    const double bind_ns_rate_ = 1.0;
    const double bind_s_rate_ = 0.5;
    const double unbind_ns_rate_ = 0.1;
    const double unbind_s_rate_ = 0.05;


 //reduced onbinding
 /*
    const double sliding_rate_ = 2.0;
    const double bind_ns_rate_ = 0.001;
    const double bind_s_rate_ = 0.5;
    const double unbind_ns_rate_ = 0.1;
    const double unbind_s_rate_ = 0.05;
*/

/* //Josh
    const double sliding_rate_ = 2.0;
    const double bind_ns_rate_ = 0.50;
    const double bind_s_rate_ = 0.0;
    const double unbind_ns_rate_ = 0.7;
    const double unbind_s_rate_ = 0.00;
*/

    EditDict possibleEditsByState_;

    void buildPossibleEditsDict();


public:

    explicit SingleStrandTopo(std::vector<std::unique_ptr<Protein>> proteins = {})
        : SubChainTopo(1, 1, std::move(proteins))   // num_sides_ = 1, radius = 1
    {
        buildPossibleEditsDict();
    }

    const EditDict& getPossibleEditsByState() const override {
        return possibleEditsByState_;
    }

protected:

    double slideRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const override;

    double bindNsRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const override;

    double bindSRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const override;

    double unbindNsRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const override;

    double unbindSRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const override;

    double switchSideRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const override;
};
