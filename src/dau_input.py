import numpy as np
import math


def m_create_bp(n, bp_matrix, bp_bias):
    term_bias = []
    if bp_bias != 0:
        term_bias = [{"coefficient": bp_bias, "polynomials": []}]
    terms_single = [{"coefficient": bp_matrix[i, i], "polynomials": [i]}
                    for i in range(n) if bp_matrix[i, i] != 0.0]
    terms_bonded = [{"coefficient": bp_matrix[i, j], "polynomials": [i, j]}
                    for i in range(n) for j in range(i+1, n) if bp_matrix[i, j] != 0.0]
    bp = {"terms": term_bias + terms_single + terms_bonded}
    return bp


def encode_dau_acceptable(tot_sc, num_sc, Wt, W, P, A=1, B=1, C=1):
    '''
    :param tot_sc: total number of possible sidechains
    :param num_sc: list of numbers of possible sidechains of each site
    :param Wt: total weight
    :param W: list of weights of possible sidechains of each site
    :param P: list of prob. of occurance of possible sidechains of each site
    :param A: weight for [one sidechain for one site] penalty terms
    :param B: weight for [weight equaling] penalty terms
    :param C: weight for [probability maximizing] objective terms
    :return: [binary_polynomial_object, penalty_binary_polynomial_object]
    '''

    W_flat = [wi for wc in W for wi in wc]
    P_flat = [pi for pc in P for pi in pc]
    log_P_flat = [math.log(pi, 2) for pi in P_flat]

    penalty_bp_matrix = np.zeros((tot_sc, tot_sc))
    penalty_bp_bias = 0.0
    # Actually bp_matrix can be 1D but we will fix bp_matrix as 2D
    # for now, so that it will be easier to combine penalty_bp_matrix
    # into bp_matrix if we want to!
    bp_matrix = np.zeros((tot_sc, tot_sc))
    bp_bias = 0.0

    # one sidechain for one site -> penalty
    cnt = 0
    for k in num_sc:
        for i in range(cnt, cnt+k):
            penalty_bp_matrix[i, i] -= 1 * A
            for j in range(i+1, cnt+k):
                penalty_bp_matrix[i, j] += 2 * A
        cnt += k
        penalty_bp_bias += 1 * A

    # weight equaling -> penalty
    for i in range(tot_sc):
        penalty_bp_matrix[i, i] += ((W_flat[i] ** 2) - 2 * Wt * W_flat[i]) * B
        for j in range(i+1, tot_sc):
            penalty_bp_matrix[i, j] += (2 * W_flat[i] * W_flat[j]) * B
    penalty_bp_bias += (Wt ** 2) * B

    # probability maximizing -> objective
    for i in range(tot_sc):
        # print(-log_P_flat[i] * C)
        bp_matrix[i, i] += -log_P_flat[i] * C

    binary_polynomial = m_create_bp(tot_sc, bp_matrix, bp_bias)
    penalty_binary_polynomial = m_create_bp(
        tot_sc, penalty_bp_matrix, penalty_bp_bias)

    return (binary_polynomial, penalty_binary_polynomial)


if __name__ == "__main__":

    A = 1  # weight for [one sidechain for one site] penalty terms
    B = 1  # weight for [weight equaling] penalty terms
    C = 1  # weight for [probability maximizing] objective terms

    tot_sc = 4  # total number of sidechains
    num_sc = [2, 2]  # list of numbers of possible sidechains of each site
    Wt = 4  # total weight
    W = [[2, 2], [2, 1]]  # list of weights of possible sidechains of each site
    # list of prob. of occurance of possible sidechains of each site
    P = [[0.2, 0.8], [0.4, 0.6]]

    (bp, penalty_bp) = encode_dau_acceptable(tot_sc, num_sc, Wt, W, P, A, B, C)

    # bp: {'terms': [{'coefficient': 2.321928094887362, 'polynomials': [0]},
    #                   {'coefficient': 0.3219280948873623, 'polynomials': [1]},
    #                   {'coefficient': 1.3219280948873622, 'polynomials': [2]},
    #                   {'coefficient': 0.7369655941662062, 'polynomials': [3]}]}

    # penalty_bp: like bp but for penalty terms
