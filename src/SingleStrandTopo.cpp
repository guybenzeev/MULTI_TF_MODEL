#include "SingleStrandTopo.h"

// StateId encoding used here (works for any num_sides_):
// 0                   => FREE
// 1..num_sides_        => NS bound on side i
// num_sides_+1..2*num_sides_ => S bound on side i
void SingleStrandTopo::buildPossibleEditsDict() {
    possibleEditsByState_.clear();

    // Make sure every stateId exists in the dict (even if empty)
    for (int stateId = 0; stateId < getNumStates(); ++stateId) {
        possibleEditsByState_[stateId] = EditList{};
    }

    // FREE state (0): can bind non-specific on any side
    {
        int freeId = 0;
        for (int side = 1; side <= num_sides_; ++side) {
            Edit e{};
            e.type = EditType::BIND_NS;
            e.strandSide = side;
            possibleEditsByState_[freeId].push_back(e);
        }
    }

    // NS states: can slide, bind specific, unbind NS, switch side
    for (int side = 1; side <= num_sides_; ++side) {
        int nsId = side;

        // Slide left
        {
            Edit e{};
            e.type = EditType::SLIDE_LEFT;
            e.strandSide = side;
            possibleEditsByState_[nsId].push_back(e);
        }

        // Slide right
        {
            Edit e{};
            e.type = EditType::SLIDE_RIGHT;
            e.strandSide = side;
            possibleEditsByState_[nsId].push_back(e);
        }

        // Bind specific
        {
            Edit e{};
            e.type = EditType::BIND_S;
            e.strandSide = side;
            possibleEditsByState_[nsId].push_back(e);
        }

        // Unbind non-specific
        {
            Edit e{};
            e.type = EditType::UNBIND_NS;
            e.strandSide = side;
            possibleEditsByState_[nsId].push_back(e);
        }

        // Switch side (only meaningful if num_sides_ > 1, but safe to include)
        for (int newSide = 1; newSide <= num_sides_; ++newSide) {
            if (newSide == side) continue;
            Edit e{};
            e.type = EditType::SWITCH_SIDE;
            e.strandSide = newSide;
            possibleEditsByState_[nsId].push_back(e);
        }
    }

    // S states: can unbind S (back to NS)
    for (int side = 1; side <= num_sides_; ++side) {
        int sId = num_sides_ + side;

        Edit e{};
        e.type = EditType::UNBIND_S;
        e.strandSide = side;
        possibleEditsByState_[sId].push_back(e);
    }
}

double SingleStrandTopo::slideRate(
    const State& currentState,
    const Edit& edit,
    const std::vector<std::unique_ptr<State>>& chain,
    const int nodeId
) const {

    int n = static_cast<int>(chain.size());

    double rate = 0.0;
    
    if(currentState.ns == 0) {
        return 0.0;
    }

    if (edit.type == EditType::SLIDE_RIGHT) {
       if (nodeId + 1 >= n) return 0.0;

        if(chain[nodeId + 1]->free != 1) { // next site must be free
            return 0.0;
        }
        return sliding_rate_;
    }

    if (edit.type == EditType::SLIDE_LEFT) {
        if (nodeId - 1 < 0) return 0.0;

        if(chain[nodeId - 1]->free != 1) { // previous site must be free
            return 0.0;
        }
        return sliding_rate_;
    }

    return rate;
}

double SingleStrandTopo::bindNsRate(
    const State& currentState,
    const Edit& edit,
    const std::vector<std::unique_ptr<State>>& chain,
    const int nodeId
) const{
    if (currentState.free == 0) {
        return 0.0;
    }
    return bind_ns_rate_;
}

double SingleStrandTopo::bindSRate(
    const State& currentState,
    const Edit& edit,
    const std::vector<std::unique_ptr<State>>& chain,
    const int nodeId
) const {
    if(currentState.ns == 0) {
        return 0.0;
    }
    return bind_s_rate_;
}

double SingleStrandTopo::unbindNsRate(
    const State& currentState,
    const Edit& edit,
    const std::vector<std::unique_ptr<State>>& chain,
    const int nodeId
) const {
    if(currentState.ns == 0) {
        return 0.0;
    }
    return unbind_ns_rate_;
}

double SingleStrandTopo::unbindSRate(
    const State& currentState,
    const Edit& edit,
    const std::vector<std::unique_ptr<State>>& chain,
    const int nodeId
) const {
    if(currentState.s == 0) {
        return 0.0;
    }
    return unbind_s_rate_;
}

double SingleStrandTopo::switchSideRate(
    const State& currentState,
    const Edit& edit,
    const std::vector<std::unique_ptr<State>>& chain,
    const int nodeId
) const {
    return 0.0;
}