#pragma once
#include <vector>
#include <memory>
#include <unordered_map>

#include "State.h"
#include "Edit.h"

class SubChain;

class Protein {

protected:
    int transition_dependancy_radius_;
    int width_;

    virtual double slideRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const = 0;
    
    virtual double bindNsRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const = 0;

    virtual double bindSRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const = 0;

    virtual double unbindNsRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const = 0;

    virtual double unbindSRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const = 0;

    virtual double switchSideRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const = 0;

public:
    Protein() = default;

    virtual ~Protein() = default;

    virtual double computeRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
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

    virtual int getWidth() const {
        return width_;
    }
};
