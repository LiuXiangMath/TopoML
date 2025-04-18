import numpy as np
from Bio import PDB
from Bio.PDB.PDBIO import PDBIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
from Feature import ProteinDNAComplex
import torch
import gc
        



residue_to_one_letter = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "DA ": "A", "DT ": "T", "DC ": "C", "DG ": "G", 
    "A  ": "A", "U  ": "U", "C  ": "C", "G  ": "G",   
    "RA ": "A", "RU ": "U", "RC ": "C", "RG ": "G",
    }


def prepare_structure(pdb,chain,resid,mut):
    def get_unique_chain_id(prot):
        ls = []
        for model in prot:
            for chain in model:
                if chain.id not in ls:
                    ls.append(chain.id)
        candidate = [ 'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z' ]
        res = []
        for a in candidate:
            if a not in ls:
                res.append(a)
                if len(res)==2:
                    return res
        print('chain error')
        return
    
    def get_combined_file(filename1,filename2,outfile):
        f = open(filename1)
        protein = f.readlines()
        f.close()
        
        f = open(filename2)
        dna = f.readlines()
        f.close()
        
        f = open(outfile,'w')
        for i in range(len(protein)-2):
            f.write(protein[i])
        for i in range(1,len(dna)-1):
            f.write(dna[i])
        f.write('TER\n')
        f.write('END\n')
        f.close()
        
        
    
    def pdb_to_fasta(pdb_file, output_fasta,resid):
        parser = PDB.PDBParser(QUIET=True)  
        structure = parser.get_structure("protein", pdb_file)
        fasta_records = []
        num = 0
        for model in structure:
            for chain in model:
                seq = []
                chain_id = chain.id
                for residue in chain:
                    if residue.id[1] ==resid:
                        break
                    else:
                        num = num + 1 
                    
        for model in structure:
            for chain in model:
                seq = []
                chain_id = chain.id
                for residue in chain:
                    if residue.id[0] == ' ':  
                        try:
                            one_letter_code = residue_to_one_letter[residue.resname]
                            seq.append(one_letter_code)
                            
                        except KeyError:
                            continue  
                
                sequence = "".join(seq)
                record = SeqRecord(Seq(sequence), id=f"Chain_{chain_id}", description=f"Extracted from {pdb_file}")
                fasta_records.append(record)
        SeqIO.write(fasta_records, output_fasta, "fasta")
        return num
    
    # work directory
    folder = './feature/'
    
    # replace noncanonical 
    NC_AminoA = {'LLP':'LYS', 'M3P':'LYS', 'MSE':'MET', 'F2F':'PHE', 'CGU':'GLU','MYL':'LYS', 'TPO':'THR', 'HSE':'HIS'}
    AminoA = ['ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','SEF','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']
    parser = PDB.PDBParser(QUIET=True)
    prot = parser.get_structure('protein', '../data/PDB/' + pdb + '.pdb')
    model_num = len(prot)
    #chain1,chain2 = get_unique_chain_id(prot)
    iresidue_ids_to_remove = []
    for iresidue in prot[0][chain]:
        if iresidue.resname in NC_AminoA:
            iresidue_id = list(iresidue.id)
            iresidue_id[0] = ' '
            iresidue.id = tuple(iresidue_id)
            iresidue.resname = NC_AminoA[iresidue.resname]
        elif iresidue.resname not in NC_AminoA and iresidue.resname not in AminoA:
            iresidue_ids_to_remove.append(iresidue.id)
                        
    for iresidue_id in iresidue_ids_to_remove:
        prot[0][chain].detach_child(iresidue_id)
    
    io = PDBIO()
    io.set_structure(prot)
    io.save(folder + 'prot.pdb')
    
    # extract protein and dna/rna
    tclfile = open(folder + 'vmd.tcl', 'w')
    tclfile.write('mol new {' + folder + 'prot.pdb} type {pdb} first 0 last 0 step 1 waitfor 1\n')
    tclfile.write('set prot [atomselect top "protein and chain ' + chain + ' and (element C or element N or element O or element S)"]\n')
    tclfile.write('$prot writepdb ' + folder + 'raw_protein.pdb\n')
    tclfile.write('set dna [atomselect top "nucleic and (element C or element N or element O or element P)"]\n')
    tclfile.write('$dna writepdb ' + folder + 'raw_dna.pdb\n')
    tclfile.write('$prot delete\n')
    tclfile.write('$dna delete\n')
    tclfile.write('exit')
    tclfile.close()
    
    os.system('vmd -dispdev text -e ' + folder + 'vmd.tcl')
    os.system('rm ' + folder + 'prot.pdb')
    os.system('rm ' + folder + 'vmd.tcl')
    
    
    
    # copy software
    os.system('cp -p ../software/jackal/jackal_64bit/bin/profix ' + folder)
    os.system('cp -p ../software/jackal/jackal_64bit/bin/scap ' + folder)
    os.system('cp -p ../software/jackal/jackal.dir ' + folder)
    os.chdir(folder)
    
    # profix protein
    os.system('./profix -fix 0 raw_protein.pdb' )
    os.system('mv raw_protein_fix.pdb wild_protein.pdb' )
    os.system('rm raw_protein.pdb')
    
    
    # make mutation
    f = open('prot_scap.list','w')
    f.write(chain + ',' + str(resid) + ',' + mut) 
    f.close()
    os.system('./scap -ini 20 -min 4 ./wild_protein.pdb prot_scap.list')
    os.system('mv ./wild_protein_scap.pdb ./mut_protein.pdb')
    
    os.system('rm scap')
    os.system('rm profix')
    os.system('rm jackal.dir')
    os.system('rm prot_scap.list')
    
    get_combined_file('wild_protein.pdb','raw_dna.pdb','wild.pdb')
    get_combined_file('mut_protein.pdb','raw_dna.pdb','mut.pdb')
    os.system('rm raw_dna.pdb')
    
    
    
    # pdb2pqr
    os.system('pdb2pqr --ff=AMBER --with-ph=7.0 --keep-chain wild.pdb wild.pqr')
    #os.system('pdb2pqr --ff=AMBER --with-ph=7.0 --keep-chain wild_protein.pdb wild_protein.pqr')
    os.system('pdb2pqr --ff=AMBER --with-ph=7.0 --keep-chain mut.pdb mut.pqr')
    
    os.system('rm wild.log')
    os.system('rm mut.log')
    os.system('rm wild.pdb')
    os.system('rm mut.pdb')
    
    # fasta
    fasta_resid1 = pdb_to_fasta('wild_protein.pdb', 'wild_protein.fasta',resid)
    fasta_resid2 = pdb_to_fasta('mut_protein.pdb','mut_protein.fasta',resid)  
    assert fasta_resid1==fasta_resid2
    return fasta_resid1


def restrict_seq_len(fasta,pos,maxl=1300):
    if len(fasta)>maxl:
        half_window = maxl // 2
        start = max(0, pos - half_window)
        end = min(len(fasta), pos + half_window)
        if end - start < maxl:
            if start == 0:
                end = min(len(fasta), start + maxl)
            else:
                start = max(0, end - maxl)
        res = fasta[start:end]
        return res
    else:
        return fasta
def get_seq(filename):
    seq = []
    for record in SeqIO.parse(filename, "fasta"):
        seq.append(str(record.seq))
    return seq[0]


def get_protein_nucleic_feature():
    f = open('../data/PRI-710.txt')
    lines = f.readlines()
    f.close()
    
    feature_path = './feature'
    if not os.path.exists(feature_path):
        os.makedirs(feature_path)
    
    for i in range(2,len(lines)):
        pdb,chain,mut,label,typ = lines[i].strip().split(',')
        residue1,num,residue2 = mut[0], int(mut[1:-1]), mut[-1]
        fea = []
        fasta_resid = prepare_structure(pdb,chain,num,residue2)
        
        # wild-aux
        filename = 'wild'
        resid = num
        pdbid = 'wild'
        MutCutNear = 12
        BindCutNear = 12
        wildprot = ProteinDNAComplex(i,filename,resid,fasta_resid,pdbid,chain,residue1,residue2,MutCutNear,BindCutNear,'aux')
        
        # mut-aux
        filename = 'mut'
        resid = num
        pdbid = 'mut'
        MutCutNear = 12
        BindCutNear = 12
        mutprot = ProteinDNAComplex(i,filename,resid,fasta_resid,pdbid,chain,residue2,residue1,MutCutNear,BindCutNear,'aux')
        
        pka_shift_pos = ProteinDNAComplex.get_pka_shift(wildprot.Pka,mutprot.Pka)
        pka_shift_neg = ProteinDNAComplex.get_pka_shift(mutprot.Pka,wildprot.Pka)
        
        # positive-aux
        aux_feature_pos = wildprot.AtomFeature + wildprot.NearResCompoFeature + wildprot.PSSM + wildprot.SSFeature + \
                          mutprot.AtomFeature + (np.array(mutprot.AtomFeature)-np.array(wildprot.AtomFeature)).tolist() + \
                          mutprot.PSSM + mutprot.SSFeature + (np.array(mutprot.SSFeature)-np.array(wildprot.SSFeature)).tolist() + \
                          pka_shift_pos
        
        # negative-aux
        #aux_feature_neg = mutprot.AtomFeature + mutprot.NearResCompoFeature + mutprot.PSSM + mutprot.SSFeature + \
        #                  wildprot.AtomFeature + (np.array(wildprot.AtomFeature)-np.array(mutprot.AtomFeature)).tolist() + \
        #                  wildprot.PSSM + wildprot.SSFeature + (np.array(wildprot.SSFeature)-np.array(mutprot.SSFeature)).tolist() + \
        #                  pka_shift_neg
        fea = fea + aux_feature_pos
        #print('aux:',len(aux_feature_pos))
        
        
        # wild-topo
        filename = 'wild'
        resid = num
        pdbid = 'wild'
        MutCutNear = 12
        BindCutNear = 12
        wildprot = ProteinDNAComplex(i,filename,resid,fasta_resid,pdbid,chain,residue1,residue2,MutCutNear,BindCutNear,'topo')
        
        # mut-topo
        filename = 'mut'
        resid = num
        pdbid = 'mut'
        MutCutNear = 12
        BindCutNear = 12
        mutprot = ProteinDNAComplex(i,filename,resid,fasta_resid,pdbid,chain,residue2,residue1,MutCutNear,BindCutNear,'topo')
        
        topo_feature_pos = wildprot.TopoFeature + mutprot.TopoFeature
        #topo_feature_neg = mutprot.TopoFeature + wildprot.TopoFeature
        #print('topo:',len(topo_feature_pos))
        fea = fea + topo_feature_pos
        os.system('rm wild_protein.pdb')
        os.system('rm mut_protein.pdb')
        os.system('rm wild.pqr')
        os.system('rm mut.pqr')
        
        
        # esm2
        device = 'cuda:2'
        esm_model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")
        esm_model = esm_model.to(device)
        
        batch_converter = alphabet.get_batch_converter()
        esm_model.eval()
        esm_feature = []
        
        # wild-esm 
        fasta = get_seq('wild_protein.fasta')
        fasta = restrict_seq_len(fasta,num)
        data = [ ('wild',fasta)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_tokens = batch_tokens.to(device)
        
        with torch.no_grad():
            results = esm_model(batch_tokens, repr_layers=[33], return_contacts=True)
        temp = results['representations'][33].cpu().numpy()
        temp = temp.reshape(-1,1280)
        esm_wild = np.mean(temp,axis=0).tolist()
        
        
        
        del batch_tokens, results, temp
        gc.collect()
        torch.cuda.empty_cache()
        
        # mut-esm
        fasta = get_seq('mut_protein.fasta')
        fasta = restrict_seq_len(fasta,num)
        data = [ ('wild',fasta)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_tokens = batch_tokens.to(device)
        with torch.no_grad():
            results = esm_model(batch_tokens, repr_layers=[33], return_contacts=True)
        temp = results['representations'][33].cpu().numpy()
        temp = temp.reshape(-1,1280)
        esm_mut = np.mean(temp,axis=0).tolist()
        
        
        del batch_tokens, results, temp
        gc.collect()
        torch.cuda.empty_cache()
        
        trans_pos = esm_wild + esm_mut
        #trans_neg = esm_mut + esm_wild
        fea = fea + trans_pos
        
        print(len(fea))
        np.save(pdb+chain+mut+'-feature.npy',np.array(fea))
        os.chdir('../..')
        


get_protein_nucleic_feature()       
        
        
        