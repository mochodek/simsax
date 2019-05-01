from math import log

def create_score_functions(sax):
    def match_func(a, b):
        if a == b == "-":
            return -1
        if a == b:
            return 1
        cutpoints = sax.cutpoints[sax.nlevels]
        min_penalty_score = 1 - max(cutpoints[len(cutpoints) // 2], cutpoints[(len(cutpoints) // 2) + 1]) / 2
        symbol_distance = min(max(-1, 1 - sax.symbol_distance(a, b)), min_penalty_score)
        return symbol_distance

    # affine logarithmic gap function
    def gap_func(x, y):  # x is gap position in seq, y is gap length
        if y == 0:  # No gap
            return 0
        elif y == 1:  # Gap open penalty
            return -2
        return - (2 + y / 4.0 + log(y) / 2.0)

    return match_func, gap_func