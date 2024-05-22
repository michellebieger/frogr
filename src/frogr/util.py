"""Utility functions for frogr"""

import typing as t

# function to compute the element factors as a ratio between selected element vs baseline element
# because ATMO only has it as either metallicity OR ratio between selected element vs H
def compute_element_factor(
    metallicity: float, ratios: t.List[t.Tuple[str, float]] = None, ratio_element="O"
) -> t.List[t.Tuple[str, float]]:
    from .io import baseline_element_w_factors

    baseline_factors = baseline_element_w_factors()

    ratio_value = baseline_factors[ratio_element]

    current_ratios = {
        k: v / ratio_value for k, v in baseline_element_w_factors().items()
    }

    new_o = ratio_value * metallicity

    ratios = ratios or []

    for elem, val in ratios:
        current_ratios[elem] = val

    new_ratios = {k: v * new_o for k, v in current_ratios.items()}

    new_ratios["H"] = baseline_factors["H"]
    new_ratios["He"] = baseline_factors["He"]

    new_ratios = {k: v / baseline_factors[k] for k, v in new_ratios.items()}

    return [(k, v) for k, v in new_ratios.items()]
