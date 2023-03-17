import os
import shutil
import argparse
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt


def read_single_scaffold_data(file_path):
    '''
    Read mol from sdf file with rdkit.
    Only one (first) structure is returned.

    :param file_path: path of sdf file
    :return (rdkit Mol object, molecular weight)
    '''

    suppl = Chem.SDMolSupplier(file_path)
    return (suppl[0], ExactMolWt(suppl[0]))


def search_id_by_structure(m, db_file):
    ''' search the id (ex. s0000000001) with corresponding structure
    '''
    print('search start...')
    with open(db_file) as f:
        line = f.readline()
        while line:
            sp = line.split()
            c_sm, c_id = sp[0], sp[1]
            c_mol = Chem.MolFromSmiles(c_sm)
            if c_mol is None:
                c_mol = Chem.MolFromSmiles(c_sm, sanitize=False)
            if c_mol is not None:
                f1 = m.HasSubstructMatch(c_mol)
                f2 = c_mol.HasSubstructMatch(m)
                if f1 and f2:
                    print('cIdx found:', c_id)
                    return c_id
            line = f.readline()

    raise Exception('Search done, data not found.')



def copy_cIdx_by_id(id, cIdx_folder_path):
    cIdx_file_name = id + '.cIdx'
    src_path = os.path.join(cIdx_folder_path, cIdx_file_name)
    dst_path = os.path.normpath('../data/cIdxs')
    shutil.copy2(src_path, dst_path)
    pass


def example():
    mol, molW = read_single_scaffold_data(
        '../data/herbs/1_Cuscuta/scaffold0.sdf')
    # mol: rdkit Mol object
    # molW: molecular weight, not using in this place

    id = search_id_by_structure(
        mol, '../data/db/scaffold_0511.structureIdx.bak')
    print('#################################################################')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Search scaffold')
    parser.add_argument('sdf_file', type=str, help="path of the .sdf file")

    args = parser.parse_args()

    mol, molW = read_single_scaffold_data(args.sdf_file)
    # mol: rdkit Mol object
    # molW: molecular weight, not using in this place

    id = search_id_by_structure(
        mol, '../data/db/scaffold_0511.structureIdx.bak')
    print('#################################################################')
