#pragma once
#include <vector>
#include <memory>
#include <unordered_map>

#include "State.h"
#include "Edit.h"

class SubChainTopo {

public:
    using EditList = std::vector<Edit>;
    using EditDict = std::unordered_map<int, EditList>; // key = stateId

protected:
    int num_sides_;
    int transition_dependancy_radius_;

    virtual double slideRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<State>>& chain,
        const int nodeId
    ) const = 0;
    
    virtual double bindNsRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<State>>& chain,
        const int nodeId
    ) const = 0;

    virtual double bindSRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<State>>& chain,
        const int nodeId
    ) const = 0;

    virtual double unbindNsRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<State>>& chain,
        const int nodeId
    ) const = 0;

    virtual double unbindSRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<State>>& chain,
        const int nodeId

    ) const = 0;

    virtual double switchSideRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<State>>& chain,
        const int nodeId

    ) const = 0;

public:
    SubChainTopo(int num_sides_val, int transition_dependancy_radius_val)
        : num_sides_(num_sides_val),
          transition_dependancy_radius_(transition_dependancy_radius_val)
    {}


    virtual ~SubChainTopo() = default;

    virtual int getNumSides() const {
        return num_sides_;
    }

    virtual const EditDict& getPossibleEditsByState() const = 0;

    const EditList& getPossibleEditsForState(int stateId) const {
            const auto& dict = getPossibleEditsByState();
            auto it = dict.find(stateId);
            static const EditList empty{};
            return (it == dict.end()) ? empty : it->second;
        }


    virtual double computeRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<State>>& chain,
        const int nodeId

    ) const {
        double rate = 0.0;

        switch(edit.type) {
            case EditType::SLIDE_RIGHT:
            case EditType::SLIDE_LEFT:
                rate = slideRate(currentState, edit, chain, nodeId);
                break;
            case EditType::BIND_NS:
                rate = bindNsRate(currentState, edit, chain, nodeId);
                break;
            case EditType::BIND_S:
                rate = bindSRate(currentState, edit, chain, nodeId);
                break;
            case EditType::UNBIND_NS:
                rate = unbindNsRate(currentState, edit, chain, nodeId);
                break;
            case EditType::UNBIND_S:
                rate = unbindSRate(currentState, edit, chain, nodeId);
                break;
            case EditType::SWITCH_SIDE:
                rate = switchSideRate(currentState, edit, chain, nodeId);
                break;
            default:
                break;
        }

        return rate;
    }

    virtual int getEffectRadius() const {
        return transition_dependancy_radius_;
    }

    virtual int getNumStates() const {
        return (2*num_sides_ + 1);
    }
};
