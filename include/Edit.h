#pragma once
#include <string>
#include "EditType.h"

struct Edit {
    EditType type;     // what kind of edit (BIND_NS, SLIDE_RIGHT, etc.)
    int strandSide;    // which side of the strand (if applicable)
    int protein;      // which protein is involved (if applicable)

    std::string toString() const {
        return editTypeToString(type) + "(side=" + std::to_string(strandSide) + ", protein=" + std::to_string(protein) + ")";
    }

    std::string typeAsString() const {
        return editTypeToString(type);
    }


private:
    static std::string editTypeToString(EditType t) {
        switch (t) {
            case EditType::SLIDE_RIGHT: return "SLIDE_RIGHT";
            case EditType::SLIDE_LEFT:  return "SLIDE_LEFT";
            case EditType::BIND_NS:     return "BIND_NS";
            case EditType::BIND_S:      return "BIND_S";
            case EditType::UNBIND_NS:   return "UNBIND_NS";
            case EditType::UNBIND_S:    return "UNBIND_S";
            case EditType::SWITCH_SIDE:return "SWITCH_SIDE";
            default:                   return "UNKNOWN";
        }
    }
};
