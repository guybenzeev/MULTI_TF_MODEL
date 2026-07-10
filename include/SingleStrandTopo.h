#pragma once
#include <unordered_map>
#include <memory>
#include <utility>
#include <vector>

#include "SubChainTopo.h"

class SingleStrandTopo : public SubChainTopo {
private: 
    void buildPossibleEditsDict();

public:

    explicit SingleStrandTopo(std::vector<std::unique_ptr<Protein>> proteins = {})
        : SubChainTopo(1, std::move(proteins))   // num_sides_ = 1
    {
        buildPossibleEditsDict();
    }

protected:

    double computeRate(
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
    ); 

    int getNumStates() const override {
        return (2 * num_sides_ * getNumProteins() + 1);
    }
};
