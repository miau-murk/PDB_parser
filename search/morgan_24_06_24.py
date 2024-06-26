from rdkit import Chem
import numpy as np 
import os

def dict_copy(dict):
    d = {}
    for i in dict.keys():
        d[i] = dict[i].copy()
    return d

def periodic_table(elem, flag="r"):
    PT = {"H": [1.00794, 1], "He": [4.0026, 4], "C": [12.0107, 12], "N": [14.0067, 14], "O": [15.9994, 16],
    "F": [18.9984, 19], "Si": [28.0855, 28], "P": [30.938, 31], "S": [32.066, 32], "Cl": [35.4527, 35.5], "Se": [78.96, 79],
    "Br": [79.904, 80], "I": [126.9045, 127]}
    try:
        if flag == "r":
            return PT[elem][1]
        else:
            return PT[elem][0]
    except:
        return 0

def hybridization(m):
    hyb = {}
    for x in m.GetAtoms():
        hyb[x.GetIdx()+1] = str(x.GetHybridization())
    return hyb

def ZE(m):
    stereo = []
    for b in m.GetBonds():
        stereo.append([b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1,
            str(b.GetBondType()), str(b.GetStereo())])
    return stereo

def stereochemistry(m):
    h = Chem.FindMolChiralCenters(m,force=True,includeUnassigned=True,useLegacyImplementation=True)
    g = {}
    for i in h:
        g[i[0]+1] = i[1]
    return g

def GetRingSystems(mol, includeSpiro=False):
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon and (includeSpiro or nInCommon>1):
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    sys = []
    for sets in systems:
        h = []
        for i in sets:
            h.append(i+1)
        sys.append(h)
    return sys

# function finds the shortest distance between two vertices of a graph
# not in use

# function takes a set of vectors in three-dimensional space as input, if they
# can uniquely form a new basis, it outputs the sign of the determinant of the transition matrix to this new basis
# needs to be improved

def find_coord_orientation(a:list):
    
    n_array = np.array(a)
    det = np.linalg.det(n_array)
    if float(det) < 0:
        return 50
    else:
        return 100

def sorting(a):
	N = len(a)
	for i in range(0, N-1):
		for j in range(0, N-i-1):
			if a[j] == a[j+1]:
				return []
			elif a[j] < a[j+1]:
				a[j], a[j+1] = a[j+1], a[j]
	return a


def contains_hydrogen(mol):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:  # Atomic number 1 corresponds to hydrogen
            return True
    return False



def read_blok(blok, mol):

    ''' Эта функция генерирует два словаря. 1ый содержит координаты, элемент и стереохимию каждого атома
        2ой содержит информацию о связях; пример:
        {1: [[1.324, 4.646, 8.5645], "N", 3], ...}    {1: [[2, 1], [5, 2]], ...} '''

    lines = blok.split("\n")
    name = lines[0]
    
    #inf = lines[3].split()
    #num_atoms, num_bonds = int(inf[0]), int(inf[1])

    num_atoms, num_bonds = mol.GetNumAtoms(), mol.GetNumBonds()

    xyz_cont = lines[4:num_atoms+4]
    bonds = lines[(num_atoms+4):(num_atoms+num_bonds + 4)]

    xyz = {}
    graph = {}
    q = 1
    #graph[1] = {}
    for i in xyz_cont:
        h = i.split()
        xyz[q] = [[float(h[0]), float(h[1]), float(h[2])], h[3], 0]
        q += 1

    for j in bonds:
        k = j.split()
        if int(k[0]) not in graph.keys():
            graph[int(k[0])] = {}
            graph[int(k[0])][int(k[1])] = int(k[2])
        else:
            graph[int(k[0])][int(k[1])] = int(k[2])

    # Adding neighbors for each atom
    N = len(xyz)
    for i in range(1, N + 1):
        if i not in list(graph.keys()):
            cont1 = []
            graph[i] = {}
        else:
            cont1 = list(graph[i].keys())
        for j in graph.keys():
            if i != j:
                cont2 = list(graph[j].keys())
                if  i in cont2 and j not in cont1:
                    graph[i][j] = graph[j][i]
    sorted_graph = dict(sorted(graph.items()))
    for i in sorted_graph.keys():
        s_pp = dict(sorted(sorted_graph[i].items()))
        sorted_graph[i] = s_pp

    return name, xyz, sorted_graph

def remove_hydrogens_without_removeHs(molecule):
    # Creates a new editable molecule
    emol = Chem.EditableMol(Chem.Mol())

    # Creates a mapping of atom indices from the old molecule to the new molecule
    index_map = {}
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() == 1:  # Skip hydrogens, which have atomic number 1
            continue
        new_index = emol.AddAtom(atom)
        index_map[atom.GetIdx()] = new_index

    # Adds bonds between hydrogen atoms in the new molecule
    for bond in molecule.GetBonds():
        begin_atom_idx = bond.GetBeginAtomIdx()
        end_atom_idx = bond.GetEndAtomIdx()

        # Check if either of the bonded atoms is a hydrogen
        if begin_atom_idx in index_map and end_atom_idx in index_map:
            begin_atom_new_idx = index_map[begin_atom_idx]
            end_atom_new_idx = index_map[end_atom_idx]
            emol.AddBond(begin_atom_new_idx, end_atom_new_idx, bond.GetBondType())

    # Gets the molecule from the editable molecule
    new_molecule = emol.GetMol()

   # Updates the atom properties and coordinates if they are available
    new_molecule.UpdatePropertyCache()

    return new_molecule


class MorganTypes:
    sum_of_types = 0
    not_heavy_atoms = ["H", "F", "Cl", "Br", "I", "D"]

    def __init__(self, m):

        if contains_hydrogen(m) == True:    # then you are searching for topology identical atoms in molecule, you shoud # this str
            m = remove_hydrogens_without_removeHs(m) #

        block = Chem.MolToMolBlock(m) 
        

        self._name, self.xyz, self.graph = read_blok(block, m) # dict, dict

        self.hybrid = hybridization(m) # dict
        self.stereo = stereochemistry(m) # [()]
        self.ZE = ZE(m) # list
        self.ring_system = GetRingSystems(m) # [[1, 2, "AROMATIC", "STEREONONE"], ...]
        self.block = block
        self.mgc = {}        

    def find_heavy_atoms(self):
        heavy_cont = []
        for i in self.graph.keys():
            if self.xyz[i][1] != "H":
                heavy_cont.append(i)
        return heavy_cont

    def atomic_configuration(self, atom, m):

        coordsys_cont = []
        neighb = m 
        if self.xyz[neighb][1] == "C":
            for neighb_C in self.graph[neighb].keys(): # # find all neighbors of "C"
                if neighb_C != atom:
                    coordsys_cont.append(neighb_C)
        if len(coordsys_cont) != 3:
            return 0
        else:
            w1 = self.mgc[coordsys_cont[0]]
            w2 = self.mgc[coordsys_cont[1]]
            w3 = self.mgc[coordsys_cont[2]]

            W = [w1, w2, w3]
            WW = sorting(W)
            if WW != []:
                n1 = W.index(WW[0]) # index of the heaviest group
                n2 = W.index(WW[1]) # 
                n3 = W.index(WW[2]) # index of the lightest group

                nn1 = coordsys_cont[n1]
                nn2 = coordsys_cont[n2]
                nn3 = coordsys_cont[n3]

                coord1 = self.xyz[nn1][0]
                coord2 = self.xyz[nn2][0]
                coord3 = self.xyz[nn3][0]

                coord = self.xyz[atom][0]

                comp1 = [coord1[0] - coord[0], coord1[1] - coord[1], coord1[2] - coord[2]]
                comp2 = [coord2[0] - coord[0], coord2[1] - coord[1], coord2[2] - coord[2]]
                comp3 = [coord3[0] - coord[0], coord3[1] - coord[1], coord3[2] - coord[2]]

                return find_coord_orientation([comp1, comp2, comp3])
            
            else:
                return 0

    def calc_morgan_weight(self, k=10, str=0, sts=0, flag_ze=True, flag_indster=False):
        for i in self.xyz.keys():

            if i in list(self.hybrid.keys()):
                if self.hybrid[i] == "SP": 
                    self.mgc[i] = 1
            else: self.mgc[i] = 0

            if i in list(self.hybrid.keys()):
                if self.hybrid[i] == "SP2":
                    self.mgc[i] = 2
            else: self.mgc[i] = 0

            if i in list(self.hybrid.keys()):
                if self.hybrid[i] == "SP3":
                    self.mgc[i] = 3
            else: self.mgc[i] = 0

            if i in list(self.hybrid.keys()):
                if self.hybrid[i] == "SP3D":
                    self.mgc[i] = 4
            else: self.mgc[i] = 0

            if i in list(self.hybrid.keys()):
                if self.hybrid[i] == "SP3D2":
                    self.mgc[i] = 5
            else: self.mgc[i] = 0

            if i in list(self.hybrid.keys()):
                if self.hybrid[i] == "S":
                    self.mgc[i] = 0
            else: self.mgc[i] = 0

            if i in list(self.hybrid.keys()):
                if self.hybrid[i] == "UNSPECIFIED":
                    self.mgc = "ERROR"

                    return self.mgc
                
        ZE_dict = {}
        for cont in self.ZE:
            tt1, tt2 = cont[0], cont[1]
            if tt1 not in ZE_dict:
                ZE_dict[tt1] = []
                ZE_dict[tt1].append(cont[3])
            else:
                ZE_dict[tt1].append(cont[3])
            
            if tt2 not in ZE_dict:
                ZE_dict[tt2] = []
                ZE_dict[tt2].append(cont[3])
            else:
                ZE_dict[tt2].append(cont[3])


        for i in range(k):

            help_dict = dict(self.mgc)
            
            for atom in self.mgc.keys():

                s = help_dict[atom]

                for neigh_atom in self.graph[atom].keys():
                    if self.xyz[neigh_atom][1] != "H":
                        s = s + help_dict[neigh_atom]
                if atom in self.stereo.keys():
                    if self.stereo[atom] == "R":
                        s += str
                    elif self.stereo[atom] == "S":
                        s += sts

                if flag_ze == True:

                    if atom in ZE_dict.keys():
                        if "STEREOZ" in ZE_dict[atom]:
                            s += 150
                        elif "STEREOE" in ZE_dict[atom]:
                            s += 300


                self.mgc[atom] = s
            

            if flag_indster == True: # check this function again; it doesn't work if you trying to check similarity between two different mol files

                help_dict = dict(self.mgc)
                for i in self.mgc.keys():
                    for m in list(self.graph[i].keys()):
                        self.mgc[i] += self.atomic_configuration(i, m)

            for i in self.mgc.keys():
                if self.xyz[i][1] != "H":
                    self.mgc[i] += periodic_table(self.xyz[i][1], flag="r")

        return self.mgc
    
    def find_topology_indentical_atoms(self):
        if self.mgc == {}:
            self.calc_morgan_weight()
        ident_cont = []
        atom_num_cont = list(self.mgc.keys())
        for atom in atom_num_cont:
            atom_weight = self.mgc[atom]
            ids = atom_num_cont.index(atom) + 1
            #for at in atom_num_cont[atom:]:
            for at in atom_num_cont[ids:]:
                at_weight = self.mgc[at]
                if atom_weight == at_weight and self.xyz[atom][1] == self.xyz[at][1]:
                    ident_cont.append((atom, at, self.xyz[atom][1]))
        return ident_cont


    def calc_moran_type(self):
        q = 0
        for i in self.calc_morgan_weight(k=20).keys():
            q += self.calc_morgan_weight(k=20)[i]
        
        return q

