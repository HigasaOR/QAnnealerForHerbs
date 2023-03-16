import os
import json
import argparse
from search_scaffold import read_single_scaffold_data, search_id_by_structure
from read_herb_data import load_scaffold_data, read_mw
from dau_input import encode_dau_acceptable


def prep_data(path_info):
    cIdx_db_file = path_info['cIdx_db_file']
    cIdx_folder_path = path_info['cIdx_folder_path']
    sdf_file = path_info['sdf_file']
    mw_file = path_info['mw_file']

    mol, molW = read_single_scaffold_data(sdf_file)

    id = search_id_by_structure(mol, cIdx_db_file)

    cIdx_file_path = os.path.join(cIdx_folder_path, id + '.cIdx')
    data = load_scaffold_data(cIdx_file_path)
    data['scaffold_weight'] = molW
    data['possible_molecular_weights'] = read_mw(mw_file)

    return data


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Pipeline v1')
    parser.add_argument('cIdx_db_file', type=str,
                        help='path of the cIdx database file')
    parser.add_argument('cIdx_folder_path', type=str,
                        help='path of the folder that stores all the cIdx files')
    parser.add_argument('sdf_file', type=str, help='path of the .sdf file')
    parser.add_argument('mw_file', type=str,
                        help='path of the molecular weight file')

    parser.add_argument('config_num', type=int,
                        help='ith config in the cIdx file (start from 0)')
    parser.add_argument(
        'mw_num', type=int, help='ith mw in the molecular weight file (start from 0)')

    parser.add_argument('A', type=float, help='weight A')
    parser.add_argument('B', type=float, help='weight B')
    parser.add_argument('C', type=float, help='weight C')

    args = parser.parse_args()

    ################################################################

    path_info = {
        'cIdx_db_file': args.cIdx_db_file,
        'cIdx_folder_path': args.cIdx_folder_path,
        'sdf_file': args.sdf_file,
        'mw_file': args.mw_file
    }

    config_num = args.config_num
    mw_num = args.mw_num

    A = args.A
    B = args.B
    C = args.C

    ################################################################

    data = prep_data(path_info=path_info)

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

    config = data['configs'][config_num]

    num_spos = config[0]
    sposs = config[1]

    c_W = list(W[i] for i in sposs)  # [[w11, w12, w13], [w21, w22], ...]
    c_P = list(P[i] for i in sposs)
    c_num_sc = list(num_sc[i] for i in sposs)  # [3, 2, ...]
    c_tot_sc = sum(c_num_sc)

    H_weight = 1.007825032
    Wt = data['possible_molecular_weights'][mw_num] - \
        (data['scaffold_weight'] - num_spos * H_weight)

    (bp, penalty_bp) = encode_dau_acceptable(
        c_tot_sc, c_num_sc, Wt, c_W, c_P, A, B, C)

    output_data = {
        'num_sc': c_num_sc,
        'bp': bp,
        'penalty_bp': penalty_bp
    }

    with open("pipeline_v1_output.json", "w") as outfile:
        json.dump(output_data, outfile, indent=2)
