#pragma once
#include <vector>
#include <memory>
#include <unordered_map>

#include "State.h"
#include "Edit.h"

class SubChain;

class Protein {

protected:
    int transition_dependency_radius_;
    int width_;

public:
    Protein(int transition_dependency_radius = 0, int width = 1)
        : transition_dependency_radius_(transition_dependency_radius),
          width_(width)
    {}

    virtual ~Protein() = default;

    virtual double computeRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId
    ) const = 0;

    virtual int getEffectRadius() const {
        return transition_dependency_radius_;
    }

    virtual int getWidth() const {
        return width_;
    }
};
