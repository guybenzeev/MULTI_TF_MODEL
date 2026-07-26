#pragma once

#include <string>

#include "EditType.h"

/**
 * @file Edit.h
 * @brief Declares the description of a single chain transition.
 */

/**
 * @struct Edit
 * @brief Describes one protein transition on the DNA strand.
 *
 * Protein identifiers are one-based; zero means that no protein applies (every edit will have a protein apply, this convention is used only to stay consistent with the @ref State.h struct).
 * Strand-side identifiers follow the same convention.
 */
struct Edit {
    EditType type;  /**< Kind of transition. */
    int strandSide; /**< One-based strand side, or zero when not applicable. */
    int protein;    /**< One-based protein identifier, or zero when not applicable. */

    /**
     * @brief Formats all edit fields for diagnostics.
     * @return A human-readable transition description.
     */
    std::string toString() const {
        return editTypeToString(type) + "(side=" + std::to_string(strandSide) + ", protein=" + std::to_string(protein) + ")";
    }

    /**
     * @brief Returns the symbolic name of the transition type.
     * @return A stable uppercase transition name.
     */
    std::string typeAsString() const {
        return editTypeToString(type);
    }

private:
    /**
     * @brief Converts an edit type to its symbolic name.
     * @param t Transition type to convert.
     * @return The uppercase name of @p t.
     */
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
