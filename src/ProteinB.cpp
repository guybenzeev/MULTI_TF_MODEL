#include "ProteinB.h"
#include "SubChain.h"

double ProteinB::slideRate(
    const State& currentState,
    const Edit& edit,
    const std::vector<std::unique_ptr<SubChain>>& chain,
    int nodeId
) const {
    int n = static_cast<int>(chain.size());

    if (!currentState.isBoundNs()) {
        return 0.0;
    }

    if (edit.type == EditType::SLIDE_RIGHT) {
        if (nodeId + 1 >= n) {
            return 0.0;
        }

        if (!chain[nodeId + 1]->getCurrentState().isFree()) {
            return 0.0;
        }

        return sliding_rate_;
    }

    if (edit.type == EditType::SLIDE_LEFT) {
        if (nodeId - 1 < 0) {
            return 0.0;
        }

        if (!chain[nodeId - 1]->getCurrentState().isFree()) {
            return 0.0;
        }

        return sliding_rate_;
    }

    return 0.0;
}

double ProteinB::computeRate(
    const State& currentState,
    const Edit& edit,
    const std::vector<std::unique_ptr<SubChain>>& chain,
    int nodeId
) const {
    switch (edit.type) {
        case EditType::SLIDE_RIGHT:
        case EditType::SLIDE_LEFT:
            return slideRate(currentState, edit, chain, nodeId);

        case EditType::BIND_NS:
            return bind_ns_rate_;

        case EditType::BIND_S:
            return bind_s_rate_;

        case EditType::UNBIND_NS:
            return unbind_ns_rate_;

        case EditType::UNBIND_S:
            return unbind_s_rate_;

        case EditType::SWITCH_SIDE:
            return switch_side_rate_;

        default:
            return 0.0;
    }
}