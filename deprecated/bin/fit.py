#!/usr/bin/env python2
#
# This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
#
# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# For more information, please refer to <http://unlicense.org>


import multiprocessing as mp
import numpy as np
from scipy.optimize import minimize
from numpy.linalg import norm
import os
from copy import deepcopy
import threading
import subprocess
import csv
import time

from subprocess import Popen, PIPE

EV_TO_HARTREE = 0.036749322
DEBYE_TO_AU = 0.393430307
HARTREE_TO_KCAL_MOL = 627.51

__DIR_PATH__ = os.path.dirname(os.path.realpath(__file__))
__RUN__ = __DIR_PATH__ + "/run_mndo99"

__PWD__ = os.getcwd()
__SCRATCH__ = __PWD__


def parse_reference(filename):

    energy     = dict()
    dipole     = dict()
    ionization = dict()

    with open(filename, 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            if "name" in row:
                continue
            molid = row[0]
            energy[molid] = float(row[1])
            ionization[molid] = float(row[2])
            dipole[molid] = float(row[3])

    return energy, ionization, dipole

__reference_energy__, __reference_ionization__, __reference_dipole__ \
                = parse_reference("/home/andersx/projects/ppqm/reference.csv")


def get_penalty( calc_energy, calc_ionization, calc_dipole):

    epsilon = 0.0001

    penalty  = 0.0
    n = 0

    # From
    W_BINDING    =  1.0 # kcal/mol^-1
    W_IONIZATION = 10.0 # eV^-1 
    W_DIPOLE     = 20.0 # Debye^-1

    for key in calc_energy.keys():

        E_diff = (calc_energy[key] - __reference_energy__[key]) * HARTREE_TO_KCAL_MOL

        I_diff = (calc_ionization[key] -  __reference_ionization__[key]) / EV_TO_HARTREE

        D_diff = (calc_dipole[key] - __reference_dipole__[key]) / DEBYE_TO_AU

        penalty += W_BINDING**2 * E_diff**2 \
                 + W_IONIZATION**2 * I_diff**2 \
                 + W_DIPOLE**2 * D_diff**2

        # print key, calc_energy[key], __reference_energy__[key]
        # print key, calc_ionization[key], __reference_ionization__[key]
        # print key, calc_dipole[key], __reference_dipole__[key]

    return penalty

def parse_master_precise(mndo_output):

    lines = mndo_output.split("\n")

    energy     = dict()
    dipole     = dict()
    ionization = dict()

    mode = "begin"
    e_scf = 0
    e_nuc = 0
    e_iso = 0
    molid = 0


    for i, line in enumerate(lines):

        # if mode == "dipole":
        if "     PRINCIPAL AXIS" in line:

            dipole[molid] = float(line.split()[5]) * DEBYE_TO_AU
            # print "DIPOLE", dipole[molid]
            mode = "begin"


        if mode == "enuc":
            if "NUCLEAR ENERGY" in line:
                e_nuc = float(line.split()[2])
                energy[molid] = (e_nuc + e_scf - e_iso) * EV_TO_HARTREE
                # print "SCF TOTAL ENERGY", energy[molid]

            if "IONIZATION ENERGY" in line:
                ionization[molid] = -1.0 * float(line.split()[2]) * EV_TO_HARTREE
                # print "IONIZATION ENERGY", ionization[molid]

                mode = "begin"

        if mode == "eisol":
            if "TOTAL ENERGY OF THE ATOM (CALC)" in line: 
                tokens = line.split()
                idx = int(tokens[0])
                e = float(tokens[2])

                eisol[idx] = e

            if "  nexmol=-1" in line:
                tokens = lines[i-5].split()
                e_scf = float(tokens[1]) 
                e_iso = np.sum(atoms * eisol)

                tokens = lines[i-1].split()
                molid = tokens[1]
                mode = "enuc"


        if mode == "atoms":

            tokens = line.split()
            if len(tokens) == 5:
                idx = int(tokens[1])
                atoms[idx] += 1.0

            if "****" in line:
                mode = "eisol"
                # print atoms
                eisol = np.zeros((30))

        if mode == "begin":
            if "   NUMBER     NUMBER               (ANGSTROMS)          (ANGSTROMS)          (ANGSTROMS" in line:
                mode = "atoms"
                atoms = np.zeros((30))

    return energy, ionization, dipole



def shell(cmd, shell=False):

    if shell:
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    else:
        cmd = cmd.split()
        p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)

    output, err = p.communicate()

    return output



class Parameters:

    def __init__(self, names, work_dir="."):

        self.names = names
        self.n = len(names)
        self.work_dir = work_dir
        self.reference_energy, self.reference_ionization, self.reference_dipole \
                = parse_reference("/home/andersx/projects/ppqm/reference.csv")

        self.input_files = []


    def setup_master(self):

        for i, filename in enumerate(self.input_files):

            tmp_scr = "scr-"+str(i)
            os.chdir(__SCRATCH__)
            os.mkdir(tmp_scr)
            os.chdir(tmp_scr)

            shell("cp "+__PWD__+"/"+filename+" .")

            os.chdir("..")

        return


    def clean_master(self):

        for i, filename in enumerate(self.input_files):

            tmp_scr = "scr-"+str(i)
            os.chdir(__SCRATCH__)
            shell("rm -r " + tmp_scr)
            os.chdir(__PWD__)

        return


    def run_inputfile(self, i, filename, rmsds, mol_count):

        # print os.getcwd()

        tmp_scr = "scr-"+str(i)
        os.chdir(tmp_scr)

        shell("cp ../fort.14 .")

        cmd = __RUN__ + " " + filename
        print cmd
        out = shell(cmd , shell=True)

        f = open(str(i)+".out", 'w')
        f.write(out)
        f.close()

        calc_energies, calc_ionization, calc_dipole = parse_master_precise(out)

        os.chdir("..")

        print len(calc_energies), cmd

        rmsds[i] = get_penalty(calc_energies, calc_ionization, calc_dipole)
        mol_count[i] = len(calc_energies.keys())


        return



    def run_mndo99_nodisk(self):

        workers = len(self.input_files)

        penalty = mp.Array("d", [0.0 for _ in xrange(workers)])
        mol_count = mp.Array("d", [0.0 for _ in xrange(workers)])

        processes = [mp.Process(target=self.run_inputfile, args=(i, self.input_files[i], penalty, mol_count)) \
                for i in xrange(workers)]

        for p in processes: p.start()
        for p in processes: p.join()

        normalized_penalty = sum(penalty) / sum(mol_count)

        return normalized_penalty



    def write_fort14(self, params):

        output = ""
        for i in range(self.n):

            if "ZP" in self.names[i]:
                zp = self.names[i] + "  " + str(params[i]) + "\n"
                zs = "ZS"+zp[2:]

                output += zs
                output += zp

            elif ("BETAP N" in self.names[i]) or ("BETAP O" in self.names[i]):
                betap = self.names[i] + "  " + str(params[i]) + "\n"
                betas = "BETAS"+betap[5:]

                output += betas
                output += betap

            else:
                output += self.names[i] + "  " + str(params[i]) + "\n"

        f = open( self.work_dir + "/fort.14", "w")
        f.write(output)
        f.close()



    def optimize(self, values):

        self.write_fort14(values)
        # self.run_mndo99()

        penalty = self.run_mndo99_nodisk()

        return penalty

    def jacobian(self, values):

        zenergy = self.optimize(values)
        print "ENERGY: %12.7f" % (zenergy)
        grad = []

        for i, p in enumerate(values):

            dparams = deepcopy(values)

            dh = 0.000001

            dparams[i] += dh
            energy_high = nv.optimize(dparams)

            dparams[i] -= (2.0 * dh)
            energy_low = nv.optimize(dparams)

            de = energy_high - energy_low

            grad.append(de/(2.0 * dh))

            s = nv.names[i]
            # print de
            print "%3i %8s  %15.7f    dE/dP = %22.10f" % \
                    (i+1, s, values[i],  de/dh)

        grad = np.array(grad)

        print "GRADIENT NORM:", norm(grad)

        print
        print " Numpy formatted values at this point:" 
        print "    values = np.array(["
        for v in values:
            print "%20.15f," % v
        print "])"
        print 



        return grad

if __name__ == "__main__":

    import argparse
    import sys

    description = """
Folder dependent. Will look for master files in current directory.
"""

    parser = argparse.ArgumentParser(
                    usage='%(prog)s [options]',
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-q', '--quantum', action='store', help='Which quantum program to parametrize towards', metavar='exe')
    parser.add_argument('-f', '--format', action='store', help='Input file format. Either mop or inp.', metavar='fmt')

    parser.add_argument('-r', '--run', action='store', help='Overwrite the run cmd', metavar='cmd')

    args = parser.parse_args()

    if args.run != "":
        __RUN__ = args.run

    # TODO get format
    input_files = [f for f in os.listdir(__PWD__) if (f.endswith('.inp') and f.startswith('master'))]
    input_files.sort()

    if len(input_files) == 0:
        print "error: no master files"
        print
        print description
        exit()


    # TODO get reference data
    reference_energy, reference_ionization, reference_dipole \
        = parse_reference("/home/andersx/projects/ppqm/reference.csv")


    from mndo import names
    from mndo import values

    nv = Parameters(names)

    nv.input_files = input_files
    nv.setup_master()

    nv.optimize(values)

    # print minimize(nv.optimize, values, jac=nv.jacobian, method="L-BFGS-B",
    #       options={"maxiter": 1000, "disp": True})

    nv.clean_master()
