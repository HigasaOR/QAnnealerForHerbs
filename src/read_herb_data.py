import argparse
import json
from rdkit.Chem.Descriptors import MolWt

from search_scaffold import read_single_scaffold_data


def read_mw(data_path):
    ''' return a list of possible molecular weights
    '''
    lines = []
    with open(data_path) as f:
        lines = f.readlines()

    mws = [float(line.strip()) for line in lines]

    return mws


def load_scaffold_data(data_path):
    ''' 
    :param data_path: path of the target .cIdx file
    :return data object containing scaffold information
    '''

    lines = None
    with open(data_path) as f:
        lines = f.readlines()

    data = {}

    lptr = 0

    # read first line
    lptr += 1
    # read second line
    lptr += 1

    # read how many configurations
    cl = lines[lptr]
    data['num_config'] = int(cl.strip())
    lptr += 1

    # read and store the configurations
    configs = []
    for nconf in range(data['num_config']):
        cl = lines[lptr]
        cl_list = cl.split()

        num_spos = int(cl_list[0])
        sposs = [int(cl_list[ns]) for ns in range(1, num_spos+1)]

        configs.append((num_spos, sposs))
        lptr += 1

    data['configs'] = configs

    # read how many substituted positions (spos)
    cl = lines[lptr]
    data['num_spos'] = int(cl.strip())
    lptr += 1

    spos_sidechain_info = []

    for i in range(data['num_spos']):
        sp_sc_info = {}

        cl = lines[lptr]
        lptr += 1
        cl_list = cl.split()
        spos_id, spos_num_sidechain = int(cl_list[0]), int(cl_list[1])

        sp_sc_info['spos_id'] = spos_id

        sc_list = []
        for j in range(spos_num_sidechain):
            sc_info = {}
            cl = lines[lptr]
            lptr += 1
            cl_list = cl.split()
            sc_info['SMILES'] = cl_list[0]
            sc_info['freq'] = float(cl_list[1])
            sc_info['num_valence_electron'] = int(cl_list[2])
            sc_info['prob'] = float(cl_list[3])
            sc_info['MW'] = float(cl_list[4])
            sc_list.append(sc_info)

        sp_sc_info['sidechain_list'] = sc_list
        spos_sidechain_info.append(sp_sc_info)
    
    # sort based on spos id
    spos_sidechain_info.sort(key=lambda x: x['spos_id'])

    data['spos_sidechain_info'] = spos_sidechain_info

    return data


def read_herb_data(sdf_path, cIdx_path, mw_path):
    _, molW = read_single_scaffold_data(sdf_path)
    mw_list = read_mw(mw_path)
    data = load_scaffold_data(cIdx_path)
    data['scaffold_weight'] = molW
    data['possible_molecular_weights'] = mw_list
    return data


def example():
    data = read_herb_data('../data/herbs/1_Cuscuta/scaffold0.sdf',
                          '../data/cIdxs/s0000013029.cIdx', '../data/mw/CA_mw.txt')

    # You can save file for manual check
    with open("s0000013029_data.json", "w") as outfile:
        json.dump(data, outfile, indent=4)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Search scaffold')
    parser.add_argument('sdf_file', type=str, help="path of the .sdf file")
    parser.add_argument('cIdx_file', type=str, help="path of the .cIdx file")
    parser.add_argument('mw_file', type=str,
                        help="path of the molecular weight file (.txt)")

    args = parser.parse_args()

    data = read_herb_data(args.sdf_file, args.cIdx_file, args.mw_file)

    # You can save file for manual check
    with open("s0000000580_data.json", "w") as outfile:
        json.dump(data, outfile, indent=4)
