import numpy as np
from Bio import PDB
from Bio.PDB.PDBIO import PDBIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.DSSP import DSSP
from Bio import SeqIO
from Bio.Blast.Applications import NcbipsiblastCommandline
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
import copy
import gudhi
import propka
import time
import os
import sys
import pandas as pd








def atmtyp_to_ele( st ):
    # C,N,O,S,P,H
    st = st.strip()
    if len(st) == 1:
        return st
    elif st[0] == 'H':
        return 'H'
    elif st[0] == 'N':
        return 'N'
    elif st[0] == 'O':
        return 'O'
    elif st[0] == 'S':
        return 'S'
    elif st[0] == 'P':
        return 'P'
    elif st[0] == 'C':
        return 'C'
    else:
        print(st, 'Not in dictionary')
        return

residue_to_one_letter = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "DA ": "A", "DT ": "T", "DC ": "C", "DG ": "G", 
    "A  ": "A", "U  ": "U", "C  ": "C", "G  ": "G",   
    "RA ": "A", "RU ": "U", "RC ": "C", "RG ": "G",
    }

class Atom:
    def __init__(self,atype,resname,chain,resid,coord,charge,radius):
        self.AType = atype
        self.Coord = coord
        self.Charge = charge
        self.Radius = radius
        self.ResName = resname
        self.ResId = resid
        self.Chain = chain
        
ele2index = {'C':0, 'N':1, 'O':2, 'S':3, 'H':4, 'P':3}
ele2index_all = {'C':0, 'N':1, 'O':2, 'S':3, 'P':4, 'H':5}


class ProteinDNAComplex:
    def __init__(self,index,filename,resid,fasta_resid,pdbid,chain,wildresname,mutresname,mutcutnear,bindcutnear,feature_type):
        self.index = index
        self.MutCutNear = mutcutnear
        self.BindCutNear = bindcutnear
        self.PdbId = pdbid
        self.FastaResId = fasta_resid
        self.filename = filename
        self.ResId = resid
        self.ChainProtein = chain 
        self.WildResName = wildresname
        self.MutResName = mutresname
        self.Atoms = []
        self.AtomCoord = []
        self.AtomCharge = []
        self.AtomArea = []
        self.AtomSolvEng = []
        self.Coulomb = []
        self.VDW = []
        self.AtomFeature = []
        self.IndexList = []
        self.Pka = []
        self.PSSM = []
        self.SSFeature = []
        self.NearResCompoFeature = []
        self.TopoFeature = []
        
        self.set_atom_coord_and_charge()
        self.set_group_index_list()
        
        if feature_type=='aux':
            self.set_aux_feature()
        elif feature_type=='topo':
            self.set_topo_feature()
        
        
    def set_atom_coord_and_charge(self):
        f = open(self.filename+'.pqr') 
        contents = f.readlines()
        f.close()
        for line in contents:
            if line[0:4]=='ATOM':
                atom = Atom(atype=atmtyp_to_ele(line[12:16]), resname=residue_to_one_letter[line[17:20]], chain=line[21], resid=int(line[22:26]),
                            coord=[float(line[30:38]),float(line[38:46]),float(line[46:54])],
                            charge=float(line[54:62]),radius=float(line[62:68]) )
                self.Atoms.append(atom)
                self.AtomCoord.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
                self.AtomCharge.append(float(line[54:62]))
                
        self.AtomCoord = np.array(self.AtomCoord)
        self.AtomCharge = np.array(self.AtomCharge)
    
    def set_group_index_list(self):
        # mutate, mutate neighbor, binding protein, binding DNA, all
        # C, N, O, S/P, H, heavy, all
        self.IndexList = [ [ [] for _ in range(7) ] for _ in range(5) ]
        self.IndexList[4].append([])
        
        # mutation site
        IndexMutSite = []
        CoordMutSite = []
        for idx,atom in enumerate(self.Atoms):
            if atom.Chain==self.ChainProtein and atom.ResId==self.ResId:
                IndexMutSite.append(idx)
                CoordMutSite.append(atom.Coord)
        CoordMutSite = np.array(CoordMutSite)
        #print('mutate atom number',len(CoordMutSite))
        
        # mutation neighbor(near residue)
        NearRes = []
        for idx, atom in enumerate(self.Atoms):
            ResID = atom.ResId
            if ResID not in NearRes and atom.Chain==self.ChainProtein and ResID!=self.ResId:
                if np.min(np.linalg.norm(atom.Coord-CoordMutSite, axis=1)) < self.MutCutNear:
                    NearRes.append(ResID)
        NearAtom = []
        for idx, atom in enumerate(self.Atoms):
            if atom.Chain==self.ChainProtein and atom.ResId in NearRes:
                NearAtom.append(idx)
        #print('near atom number',len(NearAtom))
        
        heavy = ['C', 'N', 'O', 'S', 'P']
        # mutation index
        for idx in IndexMutSite:
            self.IndexList[0][ele2index[self.Atoms[idx].AType]].append(idx)
            if self.Atoms[idx].AType in heavy:
                self.IndexList[0][5].append(idx)
            self.IndexList[0][6].append(idx)
        
        # mutation Neighbor index
        for idx in NearAtom:
            self.IndexList[1][ele2index[self.Atoms[idx].AType]].append(idx)
            if self.Atoms[idx].AType in heavy:
                self.IndexList[1][5].append(idx)
            self.IndexList[1][6].append(idx)
        
        # binding site
        ProteinAtom = {}
        DNAAtom = {}
        t = cKDTree(self.AtomCoord)
        binding_Pair = cKDTree.query_pairs(t, self.BindCutNear)
        for (i, j) in binding_Pair:
            if self.Atoms[i].Chain==self.ChainProtein and self.Atoms[j].Chain!=self.ChainProtein:
                if i not in ProteinAtom:
                    ProteinAtom[i] = 1
                if j not in DNAAtom:
                    DNAAtom[j] = 1
            elif self.Atoms[i].Chain!=self.ChainProtein and self.Atoms[j].Chain==self.ChainProtein:
                if j not in ProteinAtom:
                    ProteinAtom[j] = 1
                if i not in DNAAtom:
                    DNAAtom[i] = 1
        
        # binding protein index
        for idx in ProteinAtom.keys():
            self.IndexList[2][ele2index[self.Atoms[idx].AType]].append(idx)
            if self.Atoms[idx].AType in heavy:
                self.IndexList[2][5].append(idx)
            self.IndexList[2][6].append(idx)
        
        # binding dna index
        for idx in DNAAtom.keys():
            self.IndexList[3][ele2index[self.Atoms[idx].AType]].append(idx)
            if self.Atoms[idx].AType in heavy:
                self.IndexList[3][5].append(idx)
            self.IndexList[3][6].append(idx)
        

        # All index
        for idx, atom in enumerate(self.Atoms):
            self.IndexList[4][ele2index_all[self.Atoms[idx].AType]].append(idx)
            if atom.AType in heavy:
                self.IndexList[4][6].append(idx)
            self.IndexList[4][7].append(idx)
        
        for i in range(5):
            for j in range(7):
                self.IndexList[i][j]=np.array(self.IndexList[i][j],int)
                #print(len(self.IndexList[i][j]))
        self.IndexList[4][7] = np.array(self.IndexList[4][7],int)
    
    
    def set_atom_area_and_solveng(self,h):
        os.system('cp -p ../../software/mibpb5/mibpb5 ./')
        os.system('./mibpb5 ' + self.PdbId + ' h=' + str(h))
        
        f = open(self.PdbId + '.englist')
        contents = f.readlines()
        f.close()
        self.AtomSolvEng = [ float(eng) for eng in contents ]
        self.AtomSolvEng = np.array(self.AtomSolvEng)
        
        f = open('partition_area.txt')
        contents = f.readlines()
        f.close()
        self.AtomArea = [ float(item.split()[1]) for item in contents ]
        self.AtomArea = np.array(self.AtomArea)
        os.system('rm ' + self.PdbId + '.englist')
        os.system('rm ' + self.PdbId + '.eng')
        os.system('rm ' + self.PdbId + '.dx')
        os.system('mv partition_area.txt ' + self.PdbId + '.area')
        os.system('rm ' + self.PdbId + '.area')
        os.system('rm intersection_info.txt')
        os.system('rm grid_info.txt')
        os.system('rm bounding_box.txt')
        os.system('rm area_volume.dat')
        os.system('rm mibpb5')
        
    
    def set_atom_coulomb_and_vdw(self,VDW_Cut=10,CLM_Cut=40):
        self.VDW = [ 0 for _ in range(len(self.Atoms)) ]
        self.Coulomb = [ 0 for _ in range(len(self.Atoms)) ]
        t = cKDTree(self.AtomCoord)
        VDW_Pair = cKDTree.query_pairs(t, VDW_Cut)
        CLM_Pair  = cKDTree.query_pairs(t, CLM_Cut)
        
        # VDW
        for (i, j) in VDW_Pair:
            ei    = self.Atoms[i].AType
            ej    = self.Atoms[j].AType
            dis   = np.linalg.norm(self.AtomCoord[i]-self.AtomCoord[j])
            ratio = (self.Atoms[i].Radius+self.Atoms[j].Radius)/dis
            vdw   = np.power(ratio, 12) - 2.*np.power(ratio, 6)
            self.VDW[i] = self.VDW[i] + vdw
            self.VDW[j] = self.VDW[j] + vdw
        self.VDW = np.array(self.VDW)
        
        # Coulomb
        for (i, j) in CLM_Pair:
            ei  = self.Atoms[i].AType
            ej  = self.Atoms[j].AType
            dis = np.linalg.norm(self.AtomCoord[i]-self.AtomCoord[j])
            clb = self.Atoms[i].Charge*self.Atoms[j].Charge/dis
            self.Coulomb[i] = self.Coulomb[i] + clb
            self.Coulomb[j] = self.Coulomb[j] + clb
        self.Coulomb = np.array(self.Coulomb)
    
    
    def construct_atom_level_feature(self):
        
        # surface area
        for i in range(5):
            for j in range(7):
                self.AtomFeature.append(np.sum(self.AtomArea[self.IndexList[i][j]]))
        self.AtomFeature.append(np.sum(self.AtomArea[self.IndexList[4][7]]))
        
        
        # partial charge
        for i in range(5):
            for j in range(7):
                self.AtomFeature.append(np.sum(self.AtomCharge[self.IndexList[i][j]]))
                self.AtomFeature.append(np.sum(np.abs(self.AtomCharge[self.IndexList[i][j]])))
        self.AtomFeature.append(np.sum(self.AtomCharge[self.IndexList[4][7]]))
        self.AtomFeature.append(np.sum(np.abs(self.AtomCharge[self.IndexList[4][7]])))
        
        # coulomb interaction
        for i in range(4):
            for j in [0,1,2,3,5]:
                self.AtomFeature.append(np.sum(self.Coulomb[self.IndexList[i][j]]))
                self.AtomFeature.append(np.sum(np.abs(self.Coulomb[self.IndexList[i][j]])))
        for j in [0,1,2,3,4,6]:
            self.AtomFeature.append(np.sum(self.Coulomb[self.IndexList[4][j]]))
            self.AtomFeature.append(np.sum(np.abs(self.Coulomb[self.IndexList[4][j]])))
        
        # van der waals interaction
        for i in range(4):
            for j in [0,1,2,3,5]:
                self.AtomFeature.append(np.sum(self.VDW[self.IndexList[i][j]]))
        for j in [0,1,2,3,4,6]:
            self.AtomFeature.append(np.sum(self.VDW[self.IndexList[4][j]]))
        
        
        # electrostatic solvation free energy
        for i in range(5):
            for j in range(7):
                self.AtomFeature.append(np.sum(self.AtomSolvEng[self.IndexList[i][j]]))
        self.AtomFeature.append(np.sum(self.AtomSolvEng[self.IndexList[4][j]]))
        
    
    
    def construct_mut_neighbor_composition(self):
        def AAcharge(AA):
            if AA in ['D','E']:
                return -1.
            elif AA in ['R','H','K']:
                return 1.
            else:
                return 0.
        Hydro = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
        PolarAll = ['S','T','N','Q','R','H','K','D','E']
        PolarUncharged = ['S','T','N','Q']
        PolarPosCharged = ['R','H','K']
        PolarNegCharged = ['D','E']
        SpecialCase = ['C','U','G','P']
        AAvolume = {'A': 88.6, 'R':173.4, 'D':111.1, 'N':114.1, 'C':108.5, \
                    'E':138.4, 'Q':143.8, 'G': 60.1, 'H':153.2, 'I':166.7, \
                    'L':166.7, 'K':168.6, 'M':162.9, 'F':189.9, 'P':112.7, \
                    'S': 89.0, 'T':116.1, 'W':227.8, 'Y':193.6, 'V':140.0}
        AAhydropathy = {'A': 1.8, 'R':-4.5, 'N':-3.5, 'D':-3.5, 'C': 2.5, \
                        'E':-3.5, 'Q':-3.5, 'G':-0.4, 'H':-3.2, 'I': 4.5, \
                        'L': 3.8, 'K':-3.9, 'M': 1.9, 'F': 2.8, 'P':-1.6, \
                        'S':-0.8, 'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V': 4.2}
        AAarea = {'A':115., 'R':225., 'D':150., 'N':160., 'C':135., \
                  'E':190., 'Q':180., 'G': 75., 'H':195., 'I':175., \
                  'L':170., 'K':200., 'M':185., 'F':210., 'P':145., \
                  'S':115., 'T':140., 'W':255., 'Y':230., 'V':155.}
        AAweight = {'A': 89.094, 'R':174.203, 'N':132.119, 'D':133.104, 'C':121.154, \
                    'E':147.131, 'Q':146.146, 'G': 75.067, 'H':155.156, 'I':131.175, \
                    'L':131.175, 'K':146.189, 'M':149.208, 'F':165.192, 'P':115.132, \
                    'S':105.093, 'T':119.12 , 'W':204.228, 'Y':181.191, 'V':117.148}
        Groups = [Hydro, PolarAll, PolarUncharged, PolarPosCharged, PolarNegCharged, SpecialCase]
        
        # neighbor residue
        NearRes,NearId = [],[]
        for i in self.IndexList[1][6]:
            ResID = self.Atoms[i].Chain+str(self.Atoms[i].ResId)
            if ResID not in NearId:
                NearId.append(ResID)
                NearRes.append(self.Atoms[i].ResName)
        
        for group in Groups:
            cnt = 0.
            for AA in NearRes:
                if AA in group:
                    cnt = cnt + 1.
            self.NearResCompoFeature.append(cnt)
            self.NearResCompoFeature.append(cnt/max(1., float(len(NearRes))))
        
        Volume = [] 
        Area = []
        Weight = []
        Hydro = []
        Charge = []
        for AA in NearRes:
            Volume.append(AAvolume[AA])
            Hydro.append(AAhydropathy[AA])
            Area.append(AAarea[AA])
            Weight.append(AAweight[AA])
            Charge.append(AAcharge(AA))
        
        if len(NearRes) == 0:
            self.NearResCompoFeature.extend([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.])
        else:
            self.NearResCompoFeature.extend([np.sum(Volume), np.mean(Volume), np.var(Volume)])
            self.NearResCompoFeature.extend([np.sum(Hydro), np.mean(Hydro), np.var(Hydro)])
            self.NearResCompoFeature.extend([np.sum(Area), np.mean(Area), np.var(Area)])
            self.NearResCompoFeature.extend([np.sum(Weight), np.mean(Weight), np.var(Weight)])
        self.NearResCompoFeature.append(np.sum(Charge))
        
    
    def set_residue_pka(self):
        os.system('propka3 ' + self.PdbId + '_protein.pdb')
        f = open(self.PdbId + '_protein.pka')
        contents = f.readlines()
        f.close()
        os.system('rm ' + self.PdbId + '_protein.pka')
        
        line_start = None
        line_end = None
        for i, line in enumerate(contents):
            if 'SUMMARY OF THIS PREDICTION' in line:
                line_start = i+2
            elif line_start is not None and line[0:10] == '----------':
                line_end = i
                break 
        mut_pka = 0
        for i in range(line_start,line_end):
            temp = contents[i].split()
            if len(temp)==5:
                resname,resid,chain,pka,_ = temp
            elif len(temp)==4:
                resname,resid,chain,pka = temp[0][0:3],temp[0][3::],temp[1],temp[2]
                
            if int(resid)==int(self.ResId):
                mut_pka = float(pka)
            else:
                self.Pka.append([resname,float(pka)])
        self.Pka.append(mut_pka)
        
    
    @classmethod
    def get_pka_shift(cls,wildpka,mutpka):
        def extract_N_C(ls):
            res = []
            N,C = 0,0
            for item in ls:
                if item[0].strip()=='N+':
                    N = item[1]
                elif item[0].strip()=='C-':
                    C = item[1]
                else:
                    res.append([ item[0],item[1] ])
            return N,C,res
        wild_res_pka = wildpka[-1]
        mut_res_pka = mutpka[-1]
        wildpka = wildpka[0:-1]
        mutpka = mutpka[0:-1]
        wild_N_pka, wild_C_pka, wildpka = extract_N_C(wildpka)
        mut_N_pka , mut_C_pka, mutpka = extract_N_C(mutpka)
        assert len(wildpka)==len(mutpka)
        pKaGroup = ['ASP','GLU','ARG','LYS','HIS','CYS','TYR']
        resname = [ item[0] for item in wildpka ]
        
        wpka = [ item[1] for item in wildpka ]
        mpka = [ item[1] for item in mutpka ]
        mutpka   = np.array(mpka)
        wildpka  = np.array(wpka)
        defer = mutpka-wildpka
        feature = [ wild_res_pka, mut_res_pka, wild_N_pka, mut_N_pka,wild_C_pka,mut_C_pka,
                    np.min(np.abs(defer)), np.sum(np.abs(defer)), np.max(defer), np.min(defer), np.sum(defer) ]
        
        sum = [0]*len(pKaGroup)
        abs_sum = [0]*len(pKaGroup)
        for j in range(len(wildpka)):
            if resname[j] in pKaGroup:
                index = pKaGroup.index(resname[j])
                abs_sum[index] = abs_sum[index] + np.abs(mutpka[j]-wildpka[j])
                sum[index] = sum[index] + mutpka[j]-wildpka[j]
        feature = feature + sum + abs_sum
        return feature
    
    
    def run_blast(self,e):
        cline = NcbipsiblastCommandline(query=self.PdbId+'_protein.fasta',
                                        db='/mnt/ufs18/rs-046/guowei-search.6/XiangLiu/protein-solubility/swissprot/swissprot',
                                        #db='/mnt/ufs18/rs-046/guowei-search.6/XiangLiu/protein-solubility/refseq_protein/refseq_protein',
                                        #db='/mnt/ufs18/rs-046/guowei-search.6/XiangLiu/protein-solubility/uniref50/uniref50',
                                        num_iterations=3,
                                        evalue=e, # 1e-3
                                        max_target_seqs=500,
                                        out=self.PdbId+'.out',
                                        out_ascii_pssm=self.PdbId+'_protein.pssm',
                                        num_threads=1)
        stdout, stderr = cline()
        
        AAind = {'A':1, 'R':2, 'N':3, 'D':4, 'C':5, 'Q':6, 'E':7, 'G':8, 'H':9, 'I':10, \
                 'L':11,'K':12,'M':13,'F':14,'P':15,'S':16,'T':17,'W':18,'Y':19,'V':20}
        pssmfile = open(self.PdbId + '_protein.pssm')
        lines = pssmfile.read().splitlines()
        number = sum(len(record.seq) for record in SeqIO.parse(self.PdbId+'_protein.fasta', "fasta"))
        pssm_score1 = np.zeros([number, 20])
        pssm_score2 = np.zeros([number, 20])
        pssm_score3 = np.zeros([number, 2])
        
        for idx, line in enumerate(lines[3:3+number]):
            tmp_vec = line.split()
            if int(tmp_vec[0])==self.FastaResId+1:
                assert self.WildResName==tmp_vec[1]
            pssm_score1[idx, :] = list(map(float,tmp_vec[2:22]))
            pssm_score2[idx, :] = list(map(float,tmp_vec[22:42]))
            pssm_score3[idx, :] = list(map(float,tmp_vec[42:]))
        pssmfile.close()
        self.PSSM.append(pssm_score1[self.FastaResId-1,AAind[self.WildResName]-1])
        self.PSSM.append(np.sum(pssm_score1[self.FastaResId-1,:]))
        self.PSSM.append(np.max(pssm_score1[self.FastaResId-1,:]))
        self.PSSM.append(pssm_score2[self.FastaResId-1,AAind[self.MutResName]-1])
        self.PSSM.append(np.sum(pssm_score2[self.FastaResId-1,:]))
        self.PSSM.append(np.max(pssm_score2[self.FastaResId-1,:]))
        self.PSSM.extend(pssm_score3[self.FastaResId-1,:])
        os.system('rm ' + self.PdbId + '.out')
        
        #print('balst feature:',len(self.PSSM))
    
    def set_ss_feature(self):
        #os.system('mkdssp -i ' + self.PdbId + '.pdb -o ' + self.PdbId + '.dssp')
        ss2index = {'H':1, 'E':2, 'G':3, 'S':4, 'B':5, 'T':6, 'I':7, '-':0}
        parser = PDB.PDBParser()
        structure = parser.get_structure("protein", self.PdbId + "_protein.pdb")
        model = structure[0]
        dssp = DSSP(model, self.PdbId + '_protein.pdb',dssp='mkdssp')
        a_key = list(dssp.keys())
        value = dssp[(self.ChainProtein,(' ',self.ResId,' '))]
        self.SSFeature.append(ss2index[value[2]])
        self.SSFeature.extend(value[3::])
        
        
        os.system('../../software/spider/SPIDER2_local/misc/pred_pssm.py ' + self.PdbId + '_protein.pssm')
        spdfile = open(self.PdbId+'_protein.spd3')
        lines = spdfile.read().splitlines()
        spdfile.close()
        line = lines[self.FastaResId+1]
        psi = 0.
        phi = 0.
        pc = 0. 
        pe = 0. 
        ph = 0.
        d0, d1, d2, d3, phi, psi, d4, d5, pc, pe, ph = line.split()
        self.SSFeature.extend([ float(phi), float(psi), float(pc), float(pe), float(ph) ])
        os.system('rm ' + self.PdbId + '_protein.pssm')
        os.system('rm ' + self.PdbId + '_protein.spd3')
        #print('SSfeature:',len(self.SSFeature))
        
        
    
    def set_aux_feature(self):
        # try h for MIBPB
        h_now = 0.95 # 0.6
        while h_now<=1.0:
            try:
                self.set_atom_area_and_solveng(h_now)
                h_now = 1.1
            except Exception as e:
                print('error for ',h_now,e)
                h_now = h_now + 0.02
        
        self.set_atom_coulomb_and_vdw()
        self.construct_atom_level_feature()
        self.construct_mut_neighbor_composition()
        self.set_residue_pka()
        # try evalue for blast
        evalue = 1e-3
        pssm = False
        while pssm==False:
            try:
                self.run_blast(evalue)
                f = open(self.PdbId+'_protein.pssm')
                f.close()
                pssm = True
                break
            except Exception as e:
                print('error for ',evalue,e)
                evalue = evalue * 10
        self.set_ss_feature()
    
    
    def set_rips_h0(self,Cut=12):
        def update_fea(m,idx,N):
            if idx>=N:
                m = m + 1
                return m
            else:
                temp = np.zeros((1,N))
                temp[0,:idx] = 1
                m = m + temp
                return m
        step = 0.5
        bin_num = int(8/step)
        # mutation site, C,N,O
        rips_h0_mut = [ np.zeros((1,bin_num)) ] * 9
        for i in range(3):
            index1 = self.IndexList[0][i]
            atom1 = self.AtomCoord[index1]
            for j in range(3):
                index2 = self.IndexList[1][j]
                atom2 = self.AtomCoord[index2]
                m,n = atom1.shape[0],atom2.shape[0]
                block_dis = cdist(atom1, atom2, metric='euclidean')
                dis_m = np.full((m + n, m + n), 99999.0)
                dis_m[:m,m:] = block_dis
                dis_m[m:,:m] = block_dis.T
                for k in range(m+n):
                    dis_m[k,k] = 0
                rips_complex = gudhi.RipsComplex(distance_matrix=dis_m, max_edge_length=Cut)
                PH = rips_complex.create_simplex_tree().persistence(min_persistence=0.1)
                now = i*3+j
                for item in PH:
                    if item[0]==0:
                        dea = item[1][1]
                        if dea==np.inf:
                           dea = 20 
                        idx = int(dea/step)
                        rips_h0_mut[now] = update_fea(rips_h0_mut[now],idx,bin_num)
        
        # binding site, C,N,O
        rips_h0_bind = [ np.zeros((1,bin_num)) ] * 9
        for i in range(3):
            index1 = self.IndexList[2][i]
            atom1 = self.AtomCoord[index1]
            for j in range(3):
                index2 = self.IndexList[3][j]
                atom2 = self.AtomCoord[index2]
                m,n = atom1.shape[0],atom2.shape[0]
                block_dis = cdist(atom1, atom2, metric='euclidean')
                dis_m = np.full((m + n, m + n), 99999.0)
                dis_m[:m,m:] = block_dis
                dis_m[m:,:m] = block_dis.T
                for k in range(m+n):
                    dis_m[k,k] = 0
                rips_complex = gudhi.RipsComplex(distance_matrix=dis_m, max_edge_length=Cut)
                PH = rips_complex.create_simplex_tree().persistence(min_persistence=0.1)
                now = i*3+j
                for item in PH:
                    if item[0]==0:
                        dea = item[1][1]
                        if dea==np.inf:
                           dea = 20 
                        idx = int(dea/step)
                        rips_h0_bind[now] = update_fea(rips_h0_bind[now],idx,bin_num)
        
        
        
        rips_h0_fea = []
        for i in range(9):
            rips_h0_fea = rips_h0_fea + rips_h0_mut[i].tolist()[0]
        for i in range(9):
            rips_h0_fea = rips_h0_fea + rips_h0_bind[i].tolist()[0]
        self.TopoFeature = self.TopoFeature + rips_h0_fea
        #print('rips h0:',len(rips_h0_fea))
    
    
    def set_rips_l0(self):
        def get_statistic(value):
            if len(value)==0:
                return [0]*7
            else:
                return [ np.max(value),np.min(value),np.sum(value),np.mean(value),np.std(value),np.var(value),len(value) ]
        fil = [ 0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0]
        
        # mutation site, C,N,O
        rips_l0_mut = [ [] ]*9
        for i in range(3):
            index1 = self.IndexList[0][i]
            atom1 = self.AtomCoord[index1]
            for j in range(3):
                index2 = self.IndexList[1][j]
                atom2 = self.AtomCoord[index2]
                m,n = atom1.shape[0],atom2.shape[0]
                block_dis = cdist(atom1, atom2, metric='euclidean')
                dis_m = np.full((m + n, m + n), 99999.0)
                dis_m[:m,m:] = block_dis
                dis_m[m:,:m] = block_dis.T
                now = i*3+j
                for now_fil in fil:
                    L = copy.deepcopy(dis_m)
                    L[L<=now_fil] = -1
                    L[L>0] = 0
                    for k in range(m+n):
                        L[k,k] = - np.sum(L[k,:])
                    eig = np.linalg.eigvalsh(L)
                    eig = eig[eig>1e-8]
                    
                    rips_l0_mut[now] = rips_l0_mut[now] + get_statistic(eig)
        
        # binding site, C,N,O
        rips_l0_bind = [ [] ]*9
        for i in range(3):
            index1 = self.IndexList[2][i]
            atom1 = self.AtomCoord[index1]
            for j in range(3):
                index2 = self.IndexList[3][j]
                atom2 = self.AtomCoord[index2]
                m,n = atom1.shape[0],atom2.shape[0]
                block_dis = cdist(atom1, atom2, metric='euclidean')
                dis_m = np.full((m + n, m + n), 99999.0)
                dis_m[:m,m:] = block_dis
                dis_m[m:,:m] = block_dis.T
                now = i*3+j
                for now_fil in fil:
                    L = copy.deepcopy(dis_m)
                    L[L<=now_fil] = -1
                    L[L>0] = 0
                    for k in range(m+n):
                        L[k,k] = - np.sum(L[k,:])
                    eig = np.linalg.eigvalsh(L)
                    eig = eig[eig>1e-8]
                    
                    rips_l0_bind[now] = rips_l0_bind[now] + get_statistic(eig)
        
        rips_l0_fea = []
        for i in range(9):
            rips_l0_fea = rips_l0_fea + rips_l0_mut[i]
        for i in range(9):
            rips_l0_fea = rips_l0_fea + rips_l0_bind[i]
        self.TopoFeature = self.TopoFeature + rips_l0_fea
        #print('rips l0:',len(rips_l0_fea))
    
    
    def set_alpha_h12(self):
        def sum_max_mean(ls):
            if len(ls)==0:
                return [0,0,0]
            else:
                return [ np.sum(ls),np.max(ls),np.mean(ls) ]
        def max_min(ls):
            if len(ls)==0:
                return [0,0]
            else:
                return [ np.max(ls),np.min(ls) ]
        
        alpha_h = [ ]
        # mutation site, C,N,O
        for i in range(3):
            index1 = self.IndexList[0][i]
            atom1 = self.AtomCoord[index1]
            for j in range(3):
                index2 = self.IndexList[1][j]
                atom2 = self.AtomCoord[index2]
                atom = np.vstack([atom1,atom2])
                alpha_complex = gudhi.AlphaComplex(points=atom)
                PH = alpha_complex.create_simplex_tree().persistence(min_persistence=0.1)
                
                
                birth1,death1,length1 = [],[],[]
                birth2,death2,length2 = [],[],[]
                for item in PH:
                        
                    if item[0]==1:
                        birth1.append( item[1][0])
                        dea = item[1][1]
                        if dea==np.inf:
                            dea = 20
                        death1.append(dea)
                        length1.append(dea-item[1][0])
                    elif item[0]==2:
                        birth2.append( item[1][0])
                        dea = item[1][1]
                        if dea==np.inf:
                            dea = 20
                        death2.append(dea)
                        length2.append(dea-item[1][0])
                    
                alpha_h.extend(sum_max_mean(length1))
                alpha_h.extend(max_min(birth1))
                alpha_h.extend(max_min(death1))
                
                alpha_h.extend(sum_max_mean(length2))
                alpha_h.extend(max_min(birth2))
                alpha_h.extend(max_min(death2))
        
        # binding site, C,N,O
        for i in range(3):
            index1 = self.IndexList[2][i]
            atom1 = self.AtomCoord[index1]
            for j in range(3):
                index2 = self.IndexList[3][j]
                atom2 = self.AtomCoord[index2]
                atom = np.vstack([atom1,atom2])
                alpha_complex = gudhi.AlphaComplex(points=atom)
                PH = alpha_complex.create_simplex_tree().persistence(min_persistence=0.1)
                
                
                birth1,death1,length1 = [],[],[]
                birth2,death2,length2 = [],[],[]
                for item in PH:
                        
                    if item[0]==1:
                        birth1.append( item[1][0])
                        dea = item[1][1]
                        if dea==np.inf:
                            dea = 20
                        death1.append(dea)
                        length1.append(dea-item[1][0])
                    elif item[0]==2:
                        birth2.append( item[1][0])
                        dea = item[1][1]
                        if dea==np.inf:
                            dea = 20
                        death2.append(dea)
                        length2.append(dea-item[1][0])
                    
                alpha_h.extend(sum_max_mean(length1))
                alpha_h.extend(max_min(birth1))
                alpha_h.extend(max_min(death1))
                
                alpha_h.extend(sum_max_mean(length2))
                alpha_h.extend(max_min(birth2))
                alpha_h.extend(max_min(death2))
                
                
        # global
        index = self.IndexList[4][6]
        atom = self.AtomCoord[index]
        alpha_complex = gudhi.AlphaComplex(points=atom)
        PH = alpha_complex.create_simplex_tree().persistence(min_persistence=0.1)
        birth1,death1,length1 = [],[],[]
        birth2,death2,length2 = [],[],[]
        for item in PH:
            if item[0]==1:
                birth1.append( item[1][0])
                dea = item[1][1]
                if dea==np.inf:
                    dea = 20
                death1.append(dea)
                length1.append(dea-item[1][0])
            elif item[0]==2:
                birth2.append( item[1][0])
                dea = item[1][1]
                if dea==np.inf:
                    dea = 20
                death2.append(dea)
                length2.append(dea-item[1][0])
        
        alpha_h.extend(sum_max_mean(length1))
        alpha_h.extend(max_min(birth1))
        alpha_h.extend(max_min(death1))
                
        alpha_h.extend(sum_max_mean(length2))
        alpha_h.extend(max_min(birth2))
        alpha_h.extend(max_min(death2))
        
        
        #print('alpha12:',len(alpha_h))
        self.TopoFeature = self.TopoFeature + alpha_h
    
    def set_topo_feature(self):
        self.set_rips_h0()
        self.set_rips_l0()
        self.set_alpha_h12()
        






