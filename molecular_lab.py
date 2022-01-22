# Almog Meir, Molecular Lab Simulator, 2022
# This program simulating basic molecular biology lab operations, using main class for
# molecules, methods to manipulate these molecules and functions to perform on test tubes - molecule arrays
# in addition, the program includes virtual test runs of 3SAT algorithm using the lab and enzymatic
# implementation of finite automata

import random

# get_subsequence_index receives as input a subsequence & sequence, and returns the first occurrence index
# of the subsequence in the sequence, -1 if the subsequence doesn't exist on it


def get_subsequence_index(subsequence, sequence):
    for x in range(len(sequence) - len(subsequence) + 1):
        if sequence[x:x + len(subsequence)] == subsequence:
            return x
    return -1

# length_sort sorts a given array of molecules according to their length.
# the array is returned as new array (not in-place sort)


def length_sort(test_tube):
    sorted_mols = sorted(test_tube, key=lambda x: x.leading_length)
    return sorted_mols

# extract receives as input a molecules array and a subsequence, and returns an array including
# only molecules that include the subsequence


def extract(test_tube, subsequence):
    extracted_tube = []
    for molecule in test_tube:
        if get_subsequence_index(subsequence, molecule.leading_strand) != -1:
            extracted_tube.append(molecule)
    return extracted_tube

# amplify receives molecules array, 2 primers for leading and lagging strands, and number n of runs to perform
# the function amplify the amount of molecules that match these primers in a high probability,
# the probability is to simulate more life-realistic scenario of replication


def amplify(test_tube, leading_primer, lagging_primer, n):
    for molecule in test_tube:
        for i in range(n):
            if random.uniform(0, 1) > 0.15:
                molecule.replicate(leading_primer, lagging_primer)

# concatenate receives 2 molecules as input, and returns a new, concatenated molecule
# if the molecules can attach to eah other by a sticky end of either one of them


def concatenate(molecule1, molecule2):
    if molecule1.amount < 1 or molecule2.amount < 1:        # if any of the molecules provided doesn't exists
        return None                                         # any more in the tube, return None
    if random.uniform(0, 1) < 0.03:                         # 3% that the operation will randomly fail,
        return None                                         # regardless of match
    sticky_length1 = molecule1.get_sticky_end_length()
    sticky_length2 = molecule2.get_sticky_end_length()
    if sticky_length1 == 0 and sticky_length2 == 0:         # if both molecules don't have sticky ends,
        return None                                         # they can't concatenate, return None

    # case 1: molecule 1 has lagging sticky end
    if sticky_length1 > 0 and \
            molecule1.lagging_strand[molecule1.leading_length:] == molecule2.lagging_strand[:sticky_length1]:
        concat_mol = Molecule(molecule1.leading_length + molecule2.leading_length)
        concat_mol.set_strand(molecule1.leading_strand + molecule2.leading_strand)
        concat_mol.set_lagging_strand(molecule1.lagging_strand)
        molecule1.amount -= 1
        molecule2.amount -= 1
        return concat_mol

    # case 2: molecule 1 has leading sticky end
    if sticky_length1 < 0 and \
            molecule1.leading_strand[molecule1.lagging_length:] == molecule2.leading_strand[:abs(sticky_length1)]:
        concat_mol = Molecule(molecule1.leading_length)
        concat_mol.set_strand(molecule1.leading_strand)
        concat_mol.set_lagging_strand(molecule1.lagging_strand + molecule2.lagging_strand)
        molecule1.amount -= 1
        molecule2.amount -= 1
        return concat_mol

    # case 3: molecule 2 has lagging sticky end
    if sticky_length2 > 0 and \
            molecule2.lagging_strand[molecule2.leading_length:] == molecule1.lagging_strand[:sticky_length2]:
        concat_mol = Molecule(molecule2.leading_length + molecule1.leading_length)
        concat_mol.set_strand(molecule2.leading_strand + molecule1.leading_strand)
        concat_mol.set_lagging_strand(molecule2.lagging_strand)
        molecule1.amount -= 1
        molecule2.amount -= 1
        return concat_mol

    # case 4: molecule 2 has leading sticky end
    if sticky_length2 < 0 and \
            molecule2.leading_strand[molecule2.lagging_length:] == molecule1.leading_strand[:abs(sticky_length2)]:
        concat_mol = Molecule(molecule2.leading_length)
        concat_mol.set_strand(molecule2.leading_strand)
        concat_mol.set_lagging_strand(molecule2.lagging_strand + molecule1.lagging_strand)
        molecule1.amount -= 1
        molecule2.amount -= 1
        return concat_mol
    return None

# logical_formula structure is a list of lists, each sublist length is n_sat, with integer as variable number,
# negative for negation (as if 1 is v1, -1 is v1' / not v1)

# molecule class, included attributes are length, offset (for restriction enzymes), amount and strand.
# the lagging strand is calculated automatically given the leading strand according to base pairing rules


class Molecule:
    # initialization method, only mandatory field is length, with optional change for offset and amount
    # the molecule object is initialized with empty strands
    def __init__(self, length, offset=0, amount=1):
        self.leading_length = length
        self.lagging_length = length
        self.offset = offset
        self.amount = amount
        self.leading_strand = [None] * length
        self.lagging_strand = [None] * length

    # first option is to randomly generate the strand of the molecule, this is according to it's length
    def random_strand(self):
        bases = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        for base in range(self.leading_length):
            self.leading_strand[base] = random.choice(list(bases.values()))
            self.lagging_strand[base] = list(bases.keys())[list(bases.values()).index(self.leading_strand[base])]

    # second option is to use given strand to set the leading strand of the molecule
    # the method assumes the given array includes only valid bases (A, T, C, G)
    # the method updates the molecule length according to the given strand and completes the lagging strand

    def set_strand(self, strand):
        bases = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        self.leading_length = len(strand)
        self.lagging_length = len(strand)
        self.leading_strand = strand  # set leading strand with the provided strand
        self.lagging_strand = [None] * self.leading_length
        for base in range(self.leading_length):  # update lagging strand
            self.lagging_strand[base] = list(bases.keys())[list(bases.values()).index(self.leading_strand[base])]

    # set_offset allow to change the offset of given molecule
    def set_offset(self, offset):
        self.offset = offset

    # cleave takes a given molecule and enzyme (as molecule object as well) and returns a new molecule
    # if the enzyme can cleave the given molecule, None otherwise
    # the cleave operation consider enzyme's offset to allow sticky end to the new molecule
    def cleave(self, enzyme):
        if self.amount < 1 or enzyme.amount < 1:
            return None
        i = get_subsequence_index(enzyme.leading_strand, self.leading_strand)
        if i == -1:
            return None
        cleaved_mol = Molecule(i + enzyme.leading_length)
        cleaved_mol.set_strand(self.leading_strand[:i + enzyme.leading_length])
        cleaved_mol.lagging_strand = self.lagging_strand[:i + enzyme.leading_length + enzyme.offset]
        cleaved_mol.lagging_length = len(cleaved_mol.lagging_strand)
        # "sticky end" if enzyme has offset
        self.amount -= 1    # Remove one molecule of the type cleaved
        return cleaved_mol

    # replicate checks if the given primers match the head & tail of the given molecule, if so it will
    # double the amount of the provided molecule if there are enough primers,
    # or add the maximum possible given the amount of primers
    def replicate(self, leading_primer, lagging_primer):
        if self.leading_strand[
           self.leading_length - leading_primer.leading_length:] == leading_primer.leading_strand and \
                self.lagging_strand[:lagging_primer.leading_length] == lagging_primer.lagging_strand:
            primers_min_amount = min(leading_primer.amount, lagging_primer.amount)
            self.amount = min(self.amount*2, primers_min_amount + self.amount)
            leading_primer.amount -= min(self.amount / 2, primers_min_amount)
            lagging_primer.amount -= min(self.amount / 2, primers_min_amount)

    # sequence print the given molecule as 2 stranded base pairings
    def sequence(self):
        print(*self.leading_strand)
        print(*self.lagging_strand, "\n")

    # get_sticky_end_length return the difference between leading and lagging strand, if
    # leading is smaller than lagging strand, positive, 0 if equal and negative if
    # lagging is smaller than leading
    def get_sticky_end_length(self):
        return self.lagging_length - self.leading_length

    # set_lagging_strand receives a sequence to set as the lagging strand and set it to current molecule

    def set_lagging_strand(self, strand):
        self.lagging_strand = strand
        self.lagging_length = len(self.lagging_strand)


def main():
    x1 = Molecule(10)
    x2 = Molecule(10)
    e1 = Molecule(3)
    e2 = Molecule(3)
    x4 = Molecule(7)
    strand1 = ['G', 'A', 'A', 'A', 'C', 'A', 'C', 'G', 'T', 'G']
    strand2 = ['T', 'T', 'C', 'A', 'G', 'A', 'C', 'A', 'A', 'C']
    strand3 = ['T', 'C', 'A']
    strand4 = ['A', 'C', 'G']
    strand5 = ['A', 'A', 'C', 'A', 'C', 'G', 'A', 'A']

    x1.set_strand(strand1)
    x2.set_strand(strand2)
    e1.set_strand(strand3)
    e2.set_strand(strand4)
    e1.set_offset(2)
    e2.set_offset(2)

    x1.sequence()
    x2.sequence()
    e1.sequence()
    e2.sequence()

    x2 = x2.cleave(e1)
    x1 = x1.cleave(e2)
    x2.sequence()
    x1.sequence()

    print("x4: ")
    x4.set_strand(strand5)
    x4.sequence()

    x3 = concatenate(x2, x1)
    if x3 is not None:
        print("x3: ")
        x3.sequence()
        print(x3.leading_length, x3.lagging_length, '\n')

    x5 = concatenate(x4, x3)
    if x5 is not None:
        print("x5: ")
        x5.sequence()
        print(x5.leading_length, x5.lagging_length)


if __name__ == '__main__':
    main()
