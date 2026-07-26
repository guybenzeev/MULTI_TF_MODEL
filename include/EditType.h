#pragma once

/**
 * @file EditType.h
 * @brief Declares the transitions supported by the DNA binding model.
 */

/**
 * @enum EditType
 * @brief Identifies a state transition that a protein can perform.
 */
enum class EditType {
    BIND_NS,     /**< Bind non-specifically to a free site. */
    UNBIND_NS,   /**< Unbind from a non-specific site. */
    SWITCH_SIDE, /**< Move to another side of the strand. */
    BIND_S,      /**< Convert a non-specific binding into a specific binding. */
    UNBIND_S,    /**< Convert a specific binding into a non-specific binding. */
    SLIDE_RIGHT, /**< Move one site toward a larger strand index. */
    SLIDE_LEFT   /**< Move one site toward a smaller strand index. */
};
