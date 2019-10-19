
import rmsd
import pandas as pd
import numpy as np
import os
import subprocess

__MNDO__ = "mndo"
__MNDO__ = os.path.expanduser(__MNDO__)

__PARAMETERS_OM2__ = {
    "H": {
        "USS": -11.90627600,
        "ZS": 1.33196700,
        "BETAS": -6.98906400,
        "ALP": 2.54413410,
    },
    "C": {
        "ZP": 1.78753700,
        "BETAS": -18.98504400,
        "BETAP": -7.93412200,
        "ALP": 2.54638000,
    },
    "N": {
        "USS": -71.93212200,
        "UPP": -57.17231900,
        #"ZS": 2.25561400,
        "ZP": 2.25561400,
        #"BETAS": -20.49575800,
        "BETAP": -20.49575800,
    },
    "O": {
        "USS": -99.64430900,
        "UPP": -77.79747200,
        # "ZS": 2.69990500,
        "ZP": 2.69990500,
        # "BETAS": -32.68808200,
        "BETAP": -32.68808200,
        "ALP": 3.16060400,
    },
}

__HEADER__ = """OM2 1SCF MULLIK PRECISE charge={:} iparok=1 jprint=5
nextmol=-1
TITLE {:}"""

__ATOMLINE__ = "{:2s} {:} 0 {:} 0 {:} 0"

__DEBUG__ = "jprint=7"


def get_indexes(lines, pattern):

    idxs = []

    for i, line in enumerate(lines):
        if line.find(pattern) != -1:
            idxs.append(i)

    return idxs


def get_index(lines, pattern):
    for i, line in enumerate(lines):
        if line.find(pattern) != -1:
            return i
    return None


def reverse_enum(L):
    for index in reversed(range(len(L))):
        yield index, L[index]


def get_rev_index(lines, pattern):

    for i, line in reverse_enum(lines):
        if line.find(pattern) != -1:
            return i

    return None


def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def run_mndo_file(filename):

    cmd = __MNDO__
    cmd += " < " + filename

    lines = execute(cmd)

    molecule_lines = []

    for line in lines:

        line = line.strip()
        molecule_lines.append(line)

        if "STATISTICS FOR RUNS WITH MANY MOLECULES" in line:
            raise StopIteration

        if "COMPUTATION TIME" in line:
            yield molecule_lines
            molecule_lines = []


def calculate(filename):

    # stdout, stderr = run_mndo_file(filename)
    # stdout = stdout.decode()
    # stdout = stdout.split("\n")
    #
    # properties = get_properties(stdout)

    return properties


def get_properties(lines):
    """

    get properties of a single calculation

    """

    # TODO UNABLE TO ACHIEVE SCF CONVERGENCE
    # TODO optimization failed

    properties = {}

    # SCF energy
    idx = get_rev_index(lines, "CORE HAMILTONIAN MATR")
    idx -= 9
    line = lines[idx]
    line = line.split()
    value = line[1]
    e_scf = float(value)
    properties["e_scf"] = e_scf

    # Nuclear energy
    idx = get_rev_index(lines, "NUCLEAR ENERGY")
    line = lines[idx]
    line = line.split()
    value = line[2]
    e_nuc = float(value)
    properties["e_nuc"] = e_nuc # ev

    # eisol
    eisol = dict()
    idxs = get_indexes(lines, "EISOL")
    for idx in idxs:
        line = lines[idx]
        line = line.split()
        atom = int(line[0])
        value = line[2]
        eisol[atom] = float(value) # ev


    # Enthalpy of formation
    idx_hof = get_index(lines, "SCF HEAT OF FORMATION")
    line = lines[idx_hof]
    line = line.split("FORMATION")
    line = line[1]
    line = line.split()
    value = line[0]
    value = float(value)
    properties["h"] = value # kcal/mol

    # ionization
    idx = get_rev_index(lines, "IONIZATION ENERGY")
    line = lines[idx]
    value = line.split()[-2]
    e_ion = float(value) # ev
    properties["e_ion"] = e_ion

    # Dipole
    idx = get_rev_index(lines, "PRINCIPAL AXIS")
    line = lines[idx]
    line = line.split()
    value = line[-1]
    value = float(value) # Debye
    properties["mu"] = value

    # # optimized coordinates
    # i = get_rev_index(lines, 'CARTESIAN COORDINATES')
    # idx_atm = 1
    # idx_x = 2
    # idx_y = 3
    # idx_z = 4
    # n_skip = 4
    #
    # if i < idx_hof:
    #     i = get_rev_index(lines, 'X-COORDINATE')
    #     idx_atm = 1
    #     idx_x = 2
    #     idx_y = 4
    #     idx_z = 6
    #     n_skip = 3
    #
    # j = i + n_skip
    # symbols = []
    # coord = []
    #
    # # continue until we hit a blank line
    # while not lines[j].isspace() and lines[j].strip():
    #     l = lines[j].split()
    #     symbols.append(int(l[idx_atm]))
    #     x = l[idx_x]
    #     y = l[idx_y]
    #     z = l[idx_z]
    #     xyz = [x, y, z]
    #     xyz = [float(c) for c in xyz]
    #     coord.append(xyz)
    #     j += 1
    #
    # coord = np.array(coord)
    # properties["coord"] = coord
    # properties["atoms"] = symbols

    # input coords
    idx = get_rev_index(lines, "INPUT GEOMETRY")
    idx += 6
    atoms = []
    coord = []
    j = idx
    idx_atm = 1
    idx_x = 2
    idx_y = 3
    idx_z = 4
    # continue until we hit a blank line
    while not lines[j].isspace() and lines[j].strip():
        l = lines[j].split()
        atoms.append(int(l[idx_atm]))
        x = l[idx_x]
        y = l[idx_y]
        z = l[idx_z]
        xyz = [x, y, z]
        xyz = [float(c) for c in xyz]
        coord.append(xyz)
        j += 1

    # calculate energy
    e_iso = [eisol[a] for a in atoms]
    e_iso = np.sum(e_iso)
    energy = (e_nuc + e_scf - e_iso)

    properties["energy"] = energy

    return properties


def set_params(parameters, cwd=None):
    """
    """

    txt = get_params(parameters)

    if cwd is not None:
        os.chdir(cwd)

    f = open('fort.14', 'w')
    f.write(txt)
    f.close()

    return


def get_params(parameters):
    """
    """

    # TODO andersx has some if-statements in his writer

    txt = ""

    for atomtype in parameters:
        for key in parameters[atomtype]:
            line = "{:8s} {:2s} {:15.11f}\n".format(key, atomtype, parameters[atomtype][key])
            txt += line

    return txt


def get_inputs(atoms_list, coords_list, charges, titles, header=__HEADER__):
    """
    """

    inptxt = ""
    for atoms, coords, charge, title in zip(atoms_list, coords_list, charges, titles):
        txt = get_input(atoms, coords, charge, title, header=header)
        inptxt += txt

    return inptxt


def get_input(atoms, coords, charge, title, header=__HEADER__):
    """
    """

    txt = header.format(charge, title)
    txt += "\n"

    for atom, coord in zip(atoms, coords):
        line = __ATOMLINE__.format(atom, *coord)
        txt += line + "\n"

    txt += "\n"

    return txt


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

    return atoms_list, coord_list, charges, titles


def test_params():

    parameters = {}
    parameters["O"] = {}
    parameters["O"]["USS"] = 666.0
    filename = "_tmp_test_params"

    set_params(parameters)

    atoms = [
    'O',
    'N',
    'C',
    'N',
    'N',
    'H',]

    coords = [
        [ -0.0593325887, 1.2684201211  , 0.0095178503  ],
        [ 1.1946293712 , 1.771776509   , 0.0001229152  ],
        [ 1.9590217387 , 0.7210517427  , -0.0128069641 ],
        [ 1.2270421979 , -0.4479406483 , -0.0121559722 ],
        [ 0.0119302176 , -0.1246338293 , 0.0012973737  ],
        [ 3.0355546734 , 0.7552313348  , -0.0229864829 ],
    ]

    inptxt = get_input(atoms, coords, 0, "testfile")

    f = open(filename, 'w')
    f.write(inptxt)
    f.write(inptxt)
    f.close()

    calculations = run_mndo_file(filename)

    for lines in calculations:
        properties = get_properties(lines)

        print(properties)


    quit()

    idx = get_index(stdout, "USS")
    line = stdout[idx]
    line = line.split()

    value = float(line[-1])

    assert value == 666.0

    return


def main():

    # # Load data
    # atoms_list, coord_list, charges, titles = load_data()
    #
    # # TODO Select data
    #
    # # TODO Set input file
    # txt = get_inputs(atoms_list[10:11], coord_list[10:11], charges[10:11], titles[10:11])
    # f = open("runfile.inp", 'w')
    # f.write(txt)
    # f.close()

    # TODO set params
    set_params(__PARAMETERS_OM2__)

    # TODO Run calculation

    # TODO Parse properties

    return


if __name__ == '__main__':
    test_params()

