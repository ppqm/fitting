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


import numpy as np
from scipy.optimize import minimize
from numpy.linalg import norm
import os
from copy import deepcopy
import threading
import subprocess

from matplotlib import pyplot
import time
import seaborn as sns
import pandas as pd

from subprocess import Popen, PIPE


def shell(cmd, shell=False):

    if shell:
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    else:
        cmd = cmd.split()
        p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)

    output, err = p.communicate()
    return output


output1 = []
output2 = []
output3 = []
output4 = []

def task1():
    global output1
    cmd = ["./run_mndo99", "master1.inp"]
    output1 = subprocess.check_output(cmd)

def task2():
    global output2
    cmd = ["./run_mndo99", "master2.inp"]
    output2 = subprocess.check_output(cmd)

def task3():
    global output3
    cmd = ["./run_mndo99", "master3.inp"]
    output3 = subprocess.check_output(cmd)

def task4():
    global output4
    cmd = ["./run_mndo99", "master4.inp"]
    output4 = subprocess.check_output(cmd)


def run_mndo99_nodisk():

    t1 = threading.Thread(target=task1)
    t2 = threading.Thread(target=task2)
    t3 = threading.Thread(target=task3)
    t4 = threading.Thread(target=task4)

    t1.start()
    t2.start()
    t3.start()
    t4.start()

    t1.join()
    t2.join()
    t3.join()
    t4.join()

    out = output1 + output2 + output3 + output4

    return out

def parse_reference(filename):

    energy = np.zeros((7174))

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    for line in lines:
        tokens = line.split()
        if len(tokens) == 2:
            idx = int(tokens[0])
            energy[idx] = float(tokens[1])

    return energy


class Parameters:

    def __init__(self, names, work_dir="."):

        self.names = names
        self.n = len(names)
        self.work_dir = work_dir
        self.reference_energy = parse_reference("dsgdb7ae2.xyz")

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

    def run_mndo99(self):
        os.system("./compute.sh")

    def get_energies(self, input_file):

        f = open(input_file, "r") 
        lines = f.readlines()
        f.close()

        for line in lines:
            if "SCF BINDING ENERGY" in line:
                energy = float(line[25:42])
                return energy

        return 0.0 # Some dummy number

    def parse_mndo(self):

        energy = np.zeros((7174))
        input_files = [f for f in os.listdir(self.work_dir) if f.endswith('.log')]

        for input_file in input_files:
            e = self.get_energies(self.work_dir + "/" + input_file)
            idx = int(input_file[:4])
            energy[idx] = e

        return energy * 23.0609


    def parse_master(self, filename):

        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    
        energy = np.zeros((7174))

        indexes = []
        energies = []


        for line in lines:
            if "TITLE" in line:
                tokens = line.split()
                indexes.append(int(tokens[1]))
            elif "SCF BINDING ENERGY" in line:
                tokens = line.split()
                energies.append(float(tokens[3]))

        for i, idx in enumerate(indexes):
            energy[idx] = energies[i]

    def parse_master_precise(self, mndo_output):

        lines = mndo_output.split("\n")

        energy = np.zeros((7174))

        mode = "begin"
        e_scf = 0
        e_nuc = 0
        e_iso = 0
        molid = 0

        for i, line in enumerate(lines):

            if mode == "enuc":
                if "NUCLEAR ENERGY" in line:
                    e_nuc = float(line.split()[2])

                    energy[molid] = e_nuc + e_scf - e_iso

                    # print "SCF TOTAL ENERGY", e_nuc + e_scf - e_iso
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
                    molid = int(tokens[1]) 

                    mode = "enuc"


            if mode == "atoms":

                tokens = line.split()
                if len(tokens) == 5:
                    idx = int(tokens[1])
                    atoms[idx] += 1.0

                if "****" in line:
                    mode = "eisol"
                    # print atoms
                    eisol = np.zeros((20))


            if mode == "begin":
                if "   NUMBER     NUMBER               (ANGSTROMS)          (ANGSTROMS)          (ANGSTROMS" in line:
                    mode = "atoms"
                    atoms = np.zeros((20))

        return energy * 23.0609


        

    def get_penalty(self, calc):

        epsilon = 0.0001

        rmsd = 0.0
        n = 0
        for i in range(len(calc)):

            if (abs(calc[i]) > epsilon) and \
                (abs(self.reference_energy[i]) > epsilon):

                rmsd += (calc[i] - self.reference_energy[i])**2
                n += 1

        rmsd /= n
        return np.sqrt(rmsd)


    def optimize(self, values):

        self.write_fort14(values)
        # self.run_mndo99()
        mndo_output = run_mndo99_nodisk()

        calc_energies = self.parse_master_precise(mndo_output)

        penalty = self.get_penalty(calc_energies)

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

    # minimize(nv.optimize, values, method="Powell", 
    # minimize(nv.optimize, values, jac=nv.jacobian, method="Newton-CG", 
    minimize(nv.optimize, values, jac=nv.jacobian, method="L-BFGS-B", 
            options={"maxiter": 1000, "disp": True})


    # ydata = pd.DataFrame(dict({ "Calculated" : er, 
    #                             "Predicted" : e}))



    # rmsd = get_rmsd(er, e)
    # print "RMSD = %6.2f kcal/mol" % rmsd
    # sns.set(style="whitegrid")
    # ax = sns.lmplot(x="Calculated", y="Predicted", data=ydata)
    # ax.set(xlim=[-2500, -500], ylim=[-2500,-500])
    # ax.set(ylabel='PBE0/def2-TZVP HoF [kcal/mol]', xlabel='DFTB3 + ML-correction HoF [kcal/mol]')
    # pyplot.savefig("correlation.png")
