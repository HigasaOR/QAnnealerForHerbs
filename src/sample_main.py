import argparse
from read_herb_data import read_herb_data
from dau_input import encode_dau_acceptable

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Sample main for testing')
    parser.add_argument('sdf_file', type=str, help="path of the .sdf file")
    parser.add_argument('cIdx_file', type=str, help="path of the .cIdx file")
    parser.add_argument('mw_file', type=str,
                        help="path of the molecular weight file (.txt)")

    args = parser.parse_args()

    data = read_herb_data(args.sdf_file, args.cIdx_file, args.mw_file)

    A = 1
    B = 1
    C = 1

    # randomly choose a weight for now
    Wt = data['possible_molecular_weights'][2]
    # TODO: add for loop to consider all the possible molecular weights

    # because the target weight is the sum of sidechains
    Wt = Wt - data['scaffold_weight']
    # excluding the scaffold weight

    W = []
    P = []

    num_sc = []

    for spos in data['spos_sidechain_info']:
        spos_id = spos['spos_id']
        sidechain_list = spos['sidechain_list']

        tmp_w = []
        tmp_p = []
        for sidechain in sidechain_list:

            prob = sidechain['prob']

            tmp_w.append(sidechain['MW'])
            tmp_p.append(sidechain['prob'])

        W.append(tmp_w)
        P.append(tmp_p)
        num_sc.append(len(tmp_w))

    # W: contains ALL sidechain weights of substituted positions
    # P: contains ALL sidechain probs. of substituted positions


    # run through all possible configurations
    # 'c_' stands for current
    for config in data['configs']:

        num_spos = config[0]  # number of substituted positions
        
        # We are only considering 2 sites here to simplify things
        # TODO: remove this constraint
        if num_spos != 2:
            continue

        sposs = config[1]  # position numbers of substituted positions

        c_W = list(W[i] for i in sposs)  # [[w11, w12, w13], [w21, w22], ...]
        c_P = list(P[i] for i in sposs)
        c_num_sc = list(num_sc[i] for i in sposs)  # [3, 2, ...]

        c_tot_sc = sum(c_num_sc)

        (bp, penalty_bp) = encode_dau_acceptable(
            c_tot_sc, c_num_sc, Wt, c_W, c_P, A, B, C)
        
        print(penalty_bp)

        # do whatever is next

        exit()  # we run for only 1 case
        # TODO: remove this exit()
