#!/usr/bin/env python2

from subprocess import Popen, PIPE
import multiprocessing as mp
import numpy as np
import os
import re

__DIR_PATH__ = os.path.dirname(os.path.realpath(__file__))
__RUN__ = __DIR_PATH__ + "/run_mopac"

__PWD__ = os.getcwd()
__SCRATCH__ = __PWD__


def shell(cmd, shell=False):

    if shell:
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    else:
        cmd = cmd.split()
        p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)

    output, err = p.communicate()

    return output


ATOM_LIST = [ x.strip() for x in ['h ','he', \
      'li','be','b ','c ','n ','o ','f ','ne', \
      'na','mg','al','si','p ','s ','cl','ar', \
      'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu', \
      'zn','ga','ge','as','se','br','kr', \
      'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag', \
      'cd','in','sn','sb','te','i ','xe', \
      'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy', \
      'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt', \
      'au','hg','tl','pb','bi','po','at','rn', \
      'fr','ra','ac','th','pa','u ','np','pu'] ]


class Mopac:

    def __init__(self, input_files):

        self.workers = len(input_files)
        self.input_files = input_files
        self.eisol = {}

        self.atomic_structure = {}
        self.set_atomic_structure()

        return


    def run_filenames():

        penalty = mp.Array("d", [0.0 for _ in xrange(workers)])

        processes = [mp.Process(target=self.run_filename, \
                        args=(i, self.input_files[i], penalty)) \
                for i in xrange(self.workers)]

        for p in processes: p.start()
        for p in processes: p.join()

        quit()

        print penalty
        return penalty


    def set_atomic_structure(self):
        """ read mop files structure """

        self.atomic_structure = {}

        for filename in self.input_files:

            name = ".".join(filename.split(".")[:-1])

            self.atomic_structure[name] = []
            f = open(filename, 'r')

            for line in f:

                if "TITLE" in line:

                    atoms = []

                    line = f.next()
                    line = line.split()

                    while len(line) > 0:
                        atom = line[0]
                        atom = atom.lower()
                        atoms.append(atom)
                        line = f.next()
                        line = line.split()

                    self.atomic_structure[name].append(atoms)

            f.close()

        return


    def get_eisol(self):

        self.run_filename("dsgdb9nsd_131971.mop")
        out = "dsgdb9nsd_131971.out"
        lines = shell('grep "EISOL" '+out, shell=True)
        lines = lines.split("\n")

        for line in lines:

            line = line.split()
            if len(line) == 0: continue

            atom = int(line[0])
            atom -= 1
            atom = ATOM_LIST[atom]
            energy = float(line[2]) # in ev

            self.eisol[atom] = energy

        return


    def get_properties(self, i, filename, data):

        name = ".".join(filename.split(".")[:-1])
        self.run_filename(filename)

        print self.parse_output(name)

        return


    def run_filename(self, filename):
        shell(__RUN__ + " " + filename)
        return


    def parse_output(self, name):

        # TODO
        # find total energy
        # binding = total_energy - eisol

        energies = shell('grep "TOTAL ENERGY" ' + name + '.out ', shell=True)
        energies = re.findall(r'-*\d+\.\d+', energies)
        energies = [float(energy) for energy in energies]

        for i, energy in enumerate(energies):
            for atom in self.atomic_structure[name][i]:
                energy -= self.eisol[atom]
            energies[i] = energy


        ionizations = shell('grep "IONIZATION POTENTIAL" ' + name + '.out ', shell=True)
        ionizations = re.findall(r'-*\d+\.\d+', ionizations)
        ionizations = [float(ion) for ion in ionizations]

        

        return energies, ionization, dipole



if __name__ == "__main__":


    input_files = [f for f in os.listdir(__PWD__) if f.endswith('.mop')]
    input_files.remove("dsgdb9nsd_131971.mop")
    input_files.sort()

    if len(input_files) > 4:
        print "you crazy. no more than 4!"
        exit()


    mopac = Mopac(input_files)

    print mopac.get_eisol()
    print mopac.get_properties(0, input_files[0], [])



