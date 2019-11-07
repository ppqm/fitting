
import sklearn

import itertools
import functools
import multiprocessing as mp
import os
import subprocess
import time
import copy
import json

import numpy as np
import pandas as pd
from numpy.linalg import norm
from scipy.optimize import minimize

import rmsd
import joblib
import mndo

cachedir = '.pycache'
memory = joblib.Memory(cachedir, verbose=0)

def get_penalty(calc_properties, refs_properties, property_weights, keys=None):

    penalty = 0.0
    n = 0

    return penalty


@memory.cache
def load_data():

    reference = "../dataset-qm9/reference.csv"
    reference = pd.read_csv(reference)

    filenames = reference["name"]
    # energies = reference["binding energy"]

    atoms_list = []
    coord_list = []
    charges = []
    titles = []

    for filename in filenames:

        titles.append(filename)
        charges.append(0)

        filename = "../dataset-qm9/xyz/" + filename + ".xyz"
        atoms, coord = rmsd.get_coordinates_xyz(filename)

        atoms_list.append(atoms)
        coord_list.append(coord)

    offset = 10+100
    to_offset = 110+100

    atoms_list = atoms_list[offset:to_offset]
    coord_list = coord_list[offset:to_offset]
    charges = charges[offset:to_offset]
    titles = titles[offset:to_offset]
    reference = reference[offset:to_offset]

    return atoms_list, coord_list, charges, titles, reference


def minimize_parameters(mols_atoms, mols_coords, reference_properties, start_parameters,
    n_procs=1,
    method="PM3",
    ignore_keys=['DD2','DD3','PO1','PO2','PO3','PO9','HYF','CORE','EISOL','FN1','FN2','FN3','GSCAL','BETAS','ZS']):
    """
    """

    n_mols = len(mols_atoms)

    # Select header
    header = """{:} 1SCF MULLIK PRECISE charge={{:}} iparok=1 jprint=5
nextmol=-1
TITLE {{:}}"""

    header = header.format(method)

    filename = "_tmp_optimizer"
    txt = mndo.get_inputs(mols_atoms, mols_coords, np.zeros(n_mols), range(n_mols), header=header)
    with open(filename, 'w') as f:
        f.write(txt)

    # Select atom parameters to optimize
    atoms = [np.unique(atom) for atom in mols_atoms]
    atoms = list(itertools.chain(*atoms))
    atoms = np.unique(atoms)


    parameters_values = []
    parameters_keys = []
    parameters = {}

    # Select parameters
    for atom in atoms:
        atom_params = start_parameters[atom]

        current = {}

        for key in atom_params:

            if key in ignore_keys: continue

            value = atom_params[key]
            current[key] = value
            parameters_values.append(value)
            parameters_keys.append([atom, key])

        parameters[atom] = current


    # Define penalty func
    def penalty(params, debug=True):

        for param, key in zip(params, parameters_keys):
            parameters[key[0]][key[1]] = param

        mndo.set_params(parameters)

        properties_list = mndo.calculate(filename)
        calc_energies = np.array([properties["energy"] for properties in properties_list])

        diff = reference_properties - calc_energies
        idxs = np.argwhere(np.isnan(diff))
        diff[idxs] = 700.0

        error = np.abs(diff)
        error = error.mean()

        if debug:
            print("penalty: {:10.2f}".format(error))

        return error


    def jacobian(params, dh=0.000001, debug=True):

        # TODO Parallelt

        grad = []

        for i, p in enumerate(params):

            dparams = copy.deepcopy(params)

            dparams[i] += dh
            forward = penalty(dparams, debug=False)

            dparams[i] -= (2.0 * dh)
            backward = penalty(dparams, debug=False)

            de = forward - backward
            grad.append(de/(2.0 * dh))


        grad = np.array(grad)

        if debug:
            nm = np.linalg.norm(grad)
            print("penalty grad: {:10.2f}".format(nm))

        return grad


    start_error = penalty(parameters_values)
    start_error_grad = jacobian(parameters_values)

    quit()

    res = minimize(penalty, parameters_values,
        method="L-BFGS-B",
        jac=jacobian,
        options={"maxiter": 1000, "disp": True})

    parameters_values = res.x
    error = penalty(parameters_values)

    for param, key in zip(parameters_values, parameters_keys):
        parameters[key[0]][key[1]] = param

    end_parameters = parameters

    return end_parameters, error



def learning_curve(
    mols_atoms,
    mols_coords,
    reference_properties,
    start_parameters):

    fold_five = sklearn.model_selection.KFold(n_splits=5, random_state=42, shuffle=True)
    n_items = len(mols_atoms)
    X = list(range(n_items))

    score = []

    for train_idxs, test_idxs in fold_five.split(X):

        train_atoms = [mols_atoms[i] for i in train_idxs]
        train_coords = [mols_coords[i] for i in train_idxs]
        train_properties = reference_properties[train_idxs]

        test_atoms = [mols_atoms[i] for i in test_idxs]
        test_coords = [mols_coords[i] for i in test_idxs]
        test_properties = reference_properties[test_idxs]

        train_parameters, train_error = minimize_parameters(train_atoms, train_coords, train_properties, start_parameters)

        print(train_parameters)


        quit()



    return


def main():


    import argparse
    import sys

    description = """"""

    parser = argparse.ArgumentParser(
        usage='%(prog)s [options]',
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-f', '--format', action='store', help='', metavar='fmt')
    parser.add_argument('-s', '--settings', action='store', help='', metavar='json')
    parser.add_argument('-p', '--parameters', action='store', help='', metavar='json')
    parser.add_argument('-o', '--results_parameters', action='store', help='', metavar='json')
    parser.add_argument('--methods', action='store', help='', metavar='str')

    args = parser.parse_args()

    mols_atoms, mols_coords, mols_charges, titles, reference = load_data()
    ref_energies = reference.iloc[:,1].tolist()
    ref_energies = np.array(ref_energies)

    with open(args.parameters, 'r') as f:
        start_params = f.read()
        start_params = json.loads(start_params)

    # end_params = minimize_parameters(mols_atoms, mols_coords, ref_energies, start_params)
    end_params = learning_curve(mols_atoms, mols_coords, ref_energies, start_params)

    print(end_params)

    quit()

    # TODO select reference

    # TODO prepare input file
    filename = "_tmp_optimizer"
    txt = mndo.get_inputs(atoms_list, coord_list, charges, titles)
    f = open(filename, 'w')
    f.write(txt)
    f.close()

    # TODO prepare parameters
    parameters = np.array([
        -99.,
        -77.,
          2.,
        -32.,
          3.,
    ])
    parameter_keys = [
        ["O", "USS"],
        ["O", "UPP"],
        ["O", "ZP"],
        ["O", "BETAP"],
        ["O", "ALP"],
    ]
    parameter_dict = {}
    parameter_dict["O"] = {}


    # TODO calculate penalty
    # properties_list = mndo.calculate(filename)

    def penalty(params):

        for param, key in zip(params, parameter_keys):
            parameter_dict[key[0]][key[1]] = param

        mndo.set_params(parameter_dict)

        properties_list = mndo.calculate(filename)
        calc_energies = np.array([properties["energy"] for properties in properties_list])

        diff = ref_energies - calc_energies
        idxs = np.argwhere(np.isnan(diff))
        diff[idxs] = 700.0

        error = diff.mean()

        return error


    print(penalty(parameters))


    status = minimize(penalty, parameters,
        method="L-BFGS-B",
        options={"maxiter": 1000, "disp": True})

    print()
    print(status)

    # TODO optimize


    return


if __name__ == "__main__":
    main()
