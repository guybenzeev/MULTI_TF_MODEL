#include "Protein.h"
#include "SubChain.h"

bool Protein::checkForSpaceToBind(
    const Edit& edit,
    const std::vector<std::unique_ptr<SubChain>>& chain,
    int nodeId
) const {
    for (int i = nodeId; i < nodeId + width_; ++i) {
        if (
            i < 0 ||
            static_cast<std::size_t>(i) >= chain.size() ||
            !chain[static_cast<std::size_t>(i)]->getCurrentState().isFree()
        ) {
            return false;
        }
    }

    return true;
}