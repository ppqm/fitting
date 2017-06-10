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

__DIR_PATH__ = os.path.dirname(os.path.realpath(__file__))
__RUN__ = __DIR_PATH__ + "/run_mndo99"

__SCRATCH__ = "scr/"
__PWD__ = os.getcwd()


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

reference_energy, reference_ionization, reference_dipole \
        = parse_reference("/home/andersx/projects/ppqm/reference.csv")

def get_penalty( calc_energies, calc_ionization, calc_dipole):

    epsilon = 0.0001

    rmsd = 0.0
    n = 0

    print calc_energies.keys()

    # for i in range(len(calc)):

    #     if (abs(calc[i]) > epsilon) and \
    #         (abs(self.reference_energy[i]) > epsilon):

    #         rmsd += (calc[i] - self.reference_energy[i])**2
    #         n += 1

    # rmsd /= n
    return np.sqrt(rmsd)

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

        if mode == "dipole":
            if "PRINCIPAL AXIS" in line:

                dipole[molid] = float(line.split()[5]) * DEBYE_TO_AU
                # print "DIPOLE", dipole[molid]
                mode = "begin"


        if mode == "enuc":
            if "NUCLEAR ENERGY" in line:
                e_nuc = float(line.split()[2])
                energy[molid] = (e_nuc + e_scf - e_iso) * EV_TO_HARTREE
                # print "SCF TOTAL ENERGY", energy[molid]

            if "IONIZATION ENERGY" in line:
                ionization[molid] = float(line.split()[2]) * EV_TO_HARTREE
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

    print err

    return output




def run_inputfile(i, filename, rmsds):

    print filename

    tmp_scr = "scr-"+str(i)

    os.chdir(__SCRATCH__)
    os.mkdir(tmp_scr)
    os.chdir(tmp_scr)

    shell("cp "+__PWD__+"/"+filename+" .")

    cmd = __RUN__ + " " + filename

    out = shell(cmd , shell=True)

    f = open(str(i)+".out", 'w')
    f.write(out)
    f.close()

    calc_energies, calc_ionization, calc_dipole = parse_master_precise(out)

    os.chdir("..")
    shell("rm -r "+tmp_scr)

    print len(calc_energies), cmd
    return

    rmsds[i] = get_penalty(calc_energies, calc_ionization, calc_dipole)

    return



def run_mndo99_nodisk():

    workers = 4


    input_files = ["master1.inp",
                   "master2.inp",
                   "master3.inp",
                   "master4.inp"]

    penalty = mp.Array("d", [0.0 for _ in xrange(workers)])

    processes = [mp.Process(target=run_inputfile, args=(i, input_files[i], penalty)) \
            for i in xrange(workers)]

    for p in processes: p.start()
    for p in processes: p.join()

    print penalty
    return penalty



class Parameters:

    def __init__(self, names, work_dir="."):

        self.names = names
        self.n = len(names)
        self.work_dir = work_dir
        self.reference_energy, self.reference_ionization, self.reference_dipole \
                = parse_reference("/home/andersx/projects/ppqm/reference.csv")

        self.output1 = []
        self.output2 = []
        self.output3 = []
        self.output4 = []

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




    def get_penalty(self, calc_energies, calc_ionization, calc_dipole):

        epsilon = 0.0001

        rmsd = 0.0
        n = 0

        print calc_energies.keys()

        # for i in range(len(calc)):

        #     if (abs(calc[i]) > epsilon) and \
        #         (abs(self.reference_energy[i]) > epsilon):

        #         rmsd += (calc[i] - self.reference_energy[i])**2
        #         n += 1

        # rmsd /= n
        return np.sqrt(rmsd)


    def optimize(self, values):

        self.write_fort14(values)
        # self.run_mndo99()

        penalty = run_mndo99_nodisk()

        # calc_energies, calc_ionization, calc_dipole = parse_master_precise(mndo_output)

        # penalty = self.get_penalty(calc_energies, calc_ionization, calc_dipole)

        # print "ENERGY: %12.7f" % (penalty)
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

    from mndo import names
    # from mndo import values_optimized as values
    from mndo import values

    nv = Parameters(names)

    nv.optimize(values)

    exit()

    minimize(nv.optimize, values, jac=nv.jacobian, method="L-BFGS-B", 
            options={"maxiter": 1000, "disp": True})
