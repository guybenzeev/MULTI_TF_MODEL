#pragma once
#include <vector>
#include <memory>
#include <unordered_map>
#include <utility>

#include "State.h"
#include "Edit.h"
#include "Protein.h"

class SubChain;

class SubChainTopo {

public:
    using EditList = std::vector<Edit>;
    using EditDict = std::unordered_map<State, EditList, StateHash>;

protected:
    int num_sides_;
    std::vector<std::unique_ptr<Protein>> proteins_;
    EditDict possibleEditsByState_;


public:
    SubChainTopo(
        int num_sides_val,
        std::vector<std::unique_ptr<Protein>> proteins = {}
    )
        : num_sides_(num_sides_val),
          proteins_(std::move(proteins))
    {}


    virtual ~SubChainTopo() = default;

    virtual int getNumSides() const {
        return num_sides_;
    }
    virtual int getNumProteins() const {
        return proteins_.size();
    }
    virtual int getNumStates() const = 0;

    virtual const EditDict& getPossibleEditsByState() const {
        return possibleEditsByState_;
    }

    const EditList& getPossibleEditsForState(const State& state) const {
            const auto& dict = getPossibleEditsByState();
            auto it = dict.find(state);
            static const EditList empty{};
            return (it == dict.end()) ? empty : it->second;
        }


    
    virtual double computeRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const = 0;

    virtual const Protein* getProteinForSubChain(int nodeId) const {
        if (nodeId < 0 || static_cast<size_t>(nodeId) >= proteins_.size()) {
            return nullptr;
        }
        return proteins_[static_cast<size_t>(nodeId)].get();
    }

};
