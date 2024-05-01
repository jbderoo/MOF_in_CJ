# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 17:34:56 2022

@author: jderoo
"""

import os
import sys
import numpy as np
import re
from scipy.spatial.distance import cdist
from tempfile import NamedTemporaryFile
from copy import deepcopy



class ezPDB:


    
    def __init__(self, pdb_path, single_residue=False):
        
        if not os.path.exists(pdb_path):
            print(f'could not find pdb file {pdb_path}\nDouble check your spelling and if the file is in the right place!')
            sys.exit()
                
        self.path = pdb_path
        self.full_path = os.path.abspath(self.path)
        
        atom_num  = []
        resi      = []
        resn      = []
        x         = []
        y         = []
        z         = []
        chain     = []
        atom_type = []
        b         = []
        
        
        three2one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'ASN': 'N', 'GLN': 'Q',
                     'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F',
                     'ALA': 'A', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R',
                     'TRP': 'W', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                     'THI': 'C', 'HIP': 'H', 'HID': 'H', 'MSE': 'M', 'ASH': 'D',
                     'GLH': 'E',
                    }
        
        f = open(pdb_path)
        
        scores = []
        score_me = False
        
        for line in f:
            if line[0:4] == 'ATOM' or line[0:4] == 'HETA':
                atom_num.append(int(line[6:11]))
                resn.append(line[17:20].strip())
                resi.append(int(line[22:26]))
                chain.append(line[21])
                x.append(float(line[31:38]))
                y.append(float(line[39:46]))
                z.append(float(line[47:54]))
                atom_type.append(line[12:16].strip())
                b.append(float(line[60:66]))
                
            elif line[0:4] == 'pose' or score_me == True:
                sc = line.split()[-1]
                score_me = True
                
                if '#END_POSE_ENERGIES_TABLE' in line:
                    score_me = False
                    
                else:
                    scores.append(float(sc))
                                  
        
        f.close()  
        
        atom_num = np.array(atom_num)
        resi     = np.array(resi)
        resn     = np.array(resn)
        x        = np.array(x)
        y        = np.array(y)
        z        = np.array(z)
        chain    = np.array(chain) 
        b        = np.array(b)
        
        atom_coords = np.zeros(shape=[len(x), 3])
        atom_coords[:,0] = x 
        atom_coords[:,1] = y 
        atom_coords[:,2] = z 
        
        seq = ['']*len(np.unique(chain))
        
        chain_track = ''
        chain_idx   = -1
        r_track     = -1
        rn_track    = ''
        unk         = ''
        ri_track    = []
        ro_track    = []
        ro_map      = []
        unique_res  = []
        rdwc        = []
        rdwc_inc    = 0
        prot_len    = 0
        inc         = 1
        b_track     = []
        
        for r in range(len(resn)):
               
            
            if r_track != resi[r] or rn_track != resn[r]:
                
                rdwc_inc += 1
                
                if chain_track != chain[r]:
                    rdwc_inc = 1
                    chain_track = chain[r]
                    chain_idx += 1
                
                b_track.append(b[r])
                r_track  = resi[r]
                rn_track = resn[r]
                ri_track.append(r_track)
                
                try:
                    ro_track.append(ro_map[-1])
                except:
                    ro_track.append(inc)
                unique_res.append(inc)
                
                inc += 1
                
                
                
                #print(f'{inc-1}: added {resi[r]}_{resn[r]}')
                
                if resn[r] in three2one.keys():
                    seq[chain_idx] += three2one[resn[r]]
                    prot_len += 1

                    
                else:
                    unk += resn[r] + '_'
            
                    
            ro_map.append(inc-1)
            rdwc.append(rdwc_inc)
         
        resi_with_chain = np.array([str(resi[i]) + chain[i] for i in range(len(chain))])
        rdwc_with_chain = np.array([str(rdwc[i]) + chain[i] for i in range(len(chain))])
        
        #seq = [j for j in seq if j] 
        
        if len(seq) == 1:
            seq = seq[0]
        
        if b[0] != b[1]:
            b_track = 0
        
        if single_residue == False:
            
            self.xyz         = atom_coords
            self.seq         = seq
            self.unknown_seq = unk[:-1]
            self.total_atoms = max(atom_num)
            self.atom_map    = atom_num
            self.atom_type   = np.array(atom_type)
            self.resi_map    = np.array(resi)
            self.pdb_resi    = np.array(ri_track)
            self.ros_resi    = np.array(ro_track)
            self.uniq_resi   = unique_res
            self.ros_map     = np.array(ro_map)
            self.chain_map   = chain
            self.chains      = np.unique(chain)
            self.resn_map    = resn
            self.protein_len = prot_len
            self.rwc         = resi_with_chain
            self.rdwc        = rdwc_with_chain
            self.b_fac_map   = b
            self.conf        = b_track

            if len(scores) > 0:
                self.score       = scores[0]
                self.per_res_sc  = np.array(scores[1:])
            else:
                self.score       = 0
                self.per_res_sc  = [0]*prot_len
            
        
        else:
            self.xyz         = atom_coords
            self.seq         = seq
            self.unknown_seq = unk[:-1]
            self.total_atoms = max(atom_num)
            self.atom_map    = atom_num
            self.atom_type   = atom_type
            self.resi_map    = resi[0]
            self.pdb_resi    = ri_track[0]
            self.ros_resi    = ro_track[0]
            self.ros_map     = ro_map[0]
            self.chain_map   = chain[0]
            self.chains      = np.unique(chain)
            self.resn_map    = resn[0]
            self.protein_len = prot_len
            self.rwc         = resi_with_chain
            self.rdwc        = rdwc_with_chain
            self.b_fac_map   = b
            self.conf        = b_track
            
            if resn[0] in three2one.keys():
                self.resn    = three2one[resn[0]]
            else:
                self.resn    = resn[0]
            self.uniq_resi   = unique_res
            self.full_path   = 'temporary residue stored in memory'
            self.path        = 'temporary residue stored in memory'
        
     
    def __str__(self):
        statement = f"loaded in {self.path}\nAn updated, simplified version of SHARPEN's PDB object."
        return statement
    
    
    def neighbors(self, res_of_int, contact_cutoff_ub=5, contact_cutoff_lb=0):
        '''
        A function that finds the neighbors of a requested residue. If the pdb has multiple chains, 
        please use the convention {RESIDUE_NUM}{CHAIN} for specifity. For example, to find the neighbors
        of residue 55 on chain A, res_of_int should be 55A
        
        INPUTS:
        pdb_file : type = str
            the pdb that you want to investigate
            
        res_of_int : type = str, list, np.str_
            the residue(s) of interest to find the neighbors. Note that the res_of_int is *not* included
            in the output, since it is always near itself.
            
        OPTIONAL:
            
        contact_cutoff_ub : type = int, float
            the upper bound for the cut off distance, i.e. residues that are farther away than this number
            will not be considered a neighbor
            
        contact_cutoff_lb : type = int, float
            the lower bound for the cut off distance, i.e. residues that are closer than this number 
            will not be considered a neighbor
            
        OUTPUT:
        neighbors : type = list, dict
        if len(res_of_int) == 1, neighbors is a list of the neighbors surrounding the residue of interest.
        if len(res_of_int) > 1, neighbors.keys() is every res_of_int asked for. neighbors[r] are the
        neighbors for residue r.
        '''
        
        if self.path == 'temporary residue stored in memory':
            print('cannot find neighbors of only a single residue!')
            return 0
        
        
        if type(res_of_int) == list or isinstance(res_of_int, np.ndarray):
            if res_of_int[0] != str:
                dumby = []
                for roi in res_of_int:
                    dumby.append(str(roi))
                    
                
                res_of_int = dumby
        
        if type(res_of_int) == int:
            res_of_int = [str(res_of_int)]
            
        if type(res_of_int) == str or isinstance(res_of_int,np.str_):
            res_of_int = [res_of_int]
        
        
        
        chain = self.chain_map
        
        num_chains = len(self.chains)
        
        GM_out = {}
        
        for roi in res_of_int:
            
            
            
            roi_f1 = re.split('(\d+)', roi)
            roi_f2 = [r.strip() for r in roi_f1 if r is not None and r.strip() != '']
    
            
            if len(roi_f2) == 1:
                roi_f2.append(str(chain[0]))
                if num_chains > 1:
                    print(f'Warning: detected {num_chains} chains. Consider using the syntax "55A" if you wanted to request residue 55 to specify and chain A. If this function errors, try this fix first!')
                
            r = int(roi_f2[0])
            
            
            
            res_idx = np.where(r == self.resi_map)[0]
            
            if len(res_idx) == 0:
                print(f'could not find residue {roi}. Its not in your pdb or something broke!')
                sys.exit()
                
    
                
                
            inv_chain = roi_f2[1]
            idx       = np.where(inv_chain == chain[res_idx])[0]
            ridx      = res_idx[idx]
            
            coords = self.xyz[ridx, :]
            dists  = cdist(coords, self.xyz)
            tf_ub  = dists <= contact_cutoff_ub
            tf_lb  = dists >= contact_cutoff_lb
            tf     = tf_ub & tf_lb
            table  = np.where(tf == True)
            my_map = np.unique(table[1])
            
            
            potential_resis = self.resi_map[my_map]
            potential_chain = chain[my_map]
            
            
            potential_combo = [str(potential_resis[i]) + potential_chain[i] for i in range(len(potential_chain))]
            
            out = []
            trip_me = True
            
            for combo in potential_combo:
                if trip_me == True:
                    out.append(combo)
                    trip_me = False
                    
                if combo != out[-1]:
                    out.append(combo)
                    
            out1 = np.array(out)
            idx  = np.where(roi == out1)
            out  = list(np.delete(out1, idx))
            GM_out[roi] = out
                
        if len(GM_out.keys()) == 1:
            GM_out = list(GM_out.values())[0]
        return GM_out
        
    
    def contact_map(self, contact_cutoff=5):
        '''
        creates a contact map of a given pdb. Should work for ligands, no promises!!
        
        INPUT:
        pdb_file : type = str
            the pdb you wish to create a contact map for
            
        OPTIONAL:
        contact_cutoff : type = int, float
            the upper bound for the cut off distance, i.e. residues that are farther away than this number
            will not be considered a neighbor
            
        OUTPUT:
        cont_map : type = 2D numpy array
            the contact map, where a 0 means these residues do not make contact, and 1 means they do make it.
        '''
        
        if self.path == 'temporary residue stored in memory':
            print('cannot find neighbors of only a single residue!')
            return 0
        
        all_res_dup = self.rwc
        res_dup = ''
        all_res = []
        for res in all_res_dup:
            if res != res_dup:
                res_dup = res
                all_res.append(res)
         
        cont_map = np.zeros(shape=[len(all_res), len(all_res)])
        
        
        for i, r in enumerate(all_res):
    
            roi  = self.neighbors(r, contact_cutoff)
    
            
            for val in roi:
                idx = np.where(np.array(all_res) == val)
                cont_map[idx, i] += 1
            
            idx = np.where(np.array(all_res) == r)
            cont_map[idx, i] += 1        
            
        return cont_map
    
    def three2one(self, resi=False, verbose=False):
        '''
        if verbose == True, print the whole dictionary for use in something else.
        Otherwise, three2one acts as a dictionary lookup for a residue 3 letter code.
        '''
        
        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'ASN': 'N', 'GLN': 'Q',
             'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F',
             'ALA': 'A', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R',
             'TRP': 'W', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
             'THI': 'C', 'HIP': 'H', 'HID': 'H', 'MSE': 'M', 'ASH': 'D',
             'GLH': 'E',
             }
        
        if verbose == True:
            print("three2one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'ASN': 'N', 'GLN': 'Q',\n"  \
                  "             'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F',\n"  \
                  "             'ALA': 'A', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R',\n"  \
                  "             'TRP': 'W', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',\n" \
                  "             'THI': 'C', 'HIP': 'H', 'HID': 'H', 'MSE': 'M', 'ASH': 'D',\n"  \
                  "             'GLH': 'E'\n" \
                  "            }   ")
        
        if resi != False:
            r_list = []
  
            for r in resi:               
            
                if resi not in d.keys() and resi != True:
                    print(f'{resi} is not 1 of the 20 canonical amino acids -- cannot convert to single letter!')
            
                else:
                    r_list.append(d[resi])
                    
            return np.array(r_list)    
            
        else:
            return d
        
    def _residue_idx_(self, resi):
        if type(resi) == int:
            idx      = np.where(np.array(self.resi_map) == resi)[0]
            

            
            if len(idx) == 0:
                idx      = np.where(np.array(self.ros_map) == resi)[0]
            
        else:
            
            try:    
                idx      = np.where(np.array(self.resi_map) == int(resi))[0]

            
                if len(idx) == 0:
                    idx      = np.where(np.array(self.ros_map) == resi)[0]
                    
            except:
                
                idx = np.where(self.rwc == resi)[0]
                if len(idx) == 0:
                    print(f"couldn't find {resi}... Double check?")
                    sys.exit()
                    
        
        
        
        return idx
               
           
    def Residue(self, resi):
        '''
        

        Parameters
        ----------
        resi : int, str
            residue of interest 

        Returns
        -------
        residue
            another ezPDB containing only information about the requested residue

        '''


        idx = np.array([], dtype=int)
        
        if type(resi) == int or type(resi) == str:
            resi = [resi]

        for residue in resi:
            idx = np.concatenate((idx, self._residue_idx_(residue)), dtype = int)
        
        
        r = deepcopy(self)
        
        r.atom_type   = self.atom_type[idx]
        r.xyz         = self.xyz[idx, :]
        r.total_atoms = r.xyz.shape[0]
        r.atom_map    = self.atom_map[idx]   
        r.atom_type   = self.atom_type[idx]
        r.resi_map    = self.resi_map[idx]    
        r.pdb_resi    = self.pdb_resi[idx]    
        r.ros_resi    = self.ros_resi[idx]    
        r.ros_map     = self.ros_map[idx]     
        r.chain_map   = self.chain_map[idx]   
        r.chains      = np.unique(r.chain_map)      
        r.resn_map    = self.resn_map[idx]    
        r.protein_len = len(np.unique(r.ros_map))  # TODO make all logic match like this!
        r.rwc         = self.rwc[idx]      
        r.rdwc        = self.rdwc[idx]    
        
        if r.resn_map[0] in self.three2one().keys():
            r.resn    = self.three2one(r.resn_map[0])
        else:
            r.resn    = r.resn_map[0]
            
        r.uniq_resi   = np.unique(r.ros_map)
        r.seq         = r.resn
        r.resi        = np.unique(r.resi_map)
        r.full_path   = 'temporary residue stored in memory'
        r.path        = 'temporary residue stored in memory'
        
        
        '''
        
        Old method of doing it: write to disk, read in. Since been optimized.
        Retained in comments just incase usage of NamedTempFile is needed again.
        
        f1 = open(self.full_path)

        temp = NamedTemporaryFile(delete=False, mode='w+t')
        
        for line in f1:
            if line[0:4] in ['ATOM', 'HETA']:
                if int(line[7:12]) in atoms:
             
                    temp.writelines(line)

        temp.seek(0)
        r = ezPDB(temp.name, single_residue=True)
        temp.close()
        '''
        
        
        
        return r
    
    
    def RosettaResidue(self, resi):
        '''
        The first residue is has a residue number 1, and residues are subsequently numbered
        after that incrementing by 1. This is how RosettaCommons does their ordering.

        Parameters
        ----------
        resi : int, str
            residue of interst. Supports Rosetta numbering of Residue+chain. 
            For example residue 55 on chain A = '55A'

        Returns
        -------
        residue ezPDB object
            returns an ezPDB object for the specified residue number.

        '''
        
        if self.path == 'temporary residue stored in memory':
            print('cannot build a single residue from a single residue!')
            return 0
        
        rwc = self.rwc # np.array([str(self.ros_map[i]) + self.chain_map[i] for i in range(len(self.chain_map))])
        #print(rwc)
        if type(resi) == int:
            idx      = np.where(np.array(self.ros_map) == resi)
            if len(np.unique(np.diff(idx))) > 1:
                print(f'WARNING: found 2 possible residues for number {resi}. Consider using Rosetta numbering scheme for specificity')
                print('for example, to select residue 55 on chain A, resi = "55A"')
                
        
        else:
            
            try: 
              
                idx = np.where(self.ros_map == int(resi))[0]

                
                if len(idx) == 0:
                    idx = np.where(int(resi) == np.array(self.ros_map))[0]
                    rwc = np.array(self.ros_map)

                
                
            except:
  
                idx = np.where(rwc == resi)[0]
                
                if len(idx) == 0:
                    
                    rdwc = self.rdwc
                    idx  = np.where(rdwc == resi)[0]
                    rwc = rdwc
                    

                    
            if len(idx) == 0:
                
                #idx = np.where()
                print(f"couldn't find {resi}... Double check?")
                sys.exit()
                    
        

        
        atoms    = self.atom_map[idx]
        resn     = np.unique(self.resn_map[idx])
        resi     = np.unique(self.resi_map[idx])
        
        f1 = open(self.full_path)

        temp = NamedTemporaryFile(delete=False, mode='w+t')
        
        for line in f1:
            if line[0:4] in ['ATOM', 'HETA']:
                if int(line[7:12]) in atoms:
             
                    temp.writelines(line)
                    
            if f'{resn}_{resi}' in line:
                score = float(line.split()[-1])
                
            else:
                score = 0

        temp.seek(0)
        r = ezPDB(temp.name, single_residue=True)
        temp.close()
        r.score = score
        
        
        return r
    
    def WritePDB(self, file_name, append=False):
        '''
        Write the current ezPDB object to a PDB file

        Parameters
        ----------
        file_name : str (file name / file path)
            The name of the file to be written 
        append : True/False, optional
            If append is True, attach the current ezPDB object to the end of 
            a currently existing PDB file. The default is False.

        Returns
        -------
        A pdb file written to disk.

        '''

        atom_num  = self.atom_map
        atom_type = self.atom_type
        chain     = self.chain_map
        resn      = self.resn_map
        resi      = self.resi_map 
        xyz       = self.xyz
        x         = xyz[:, 0]
        y         = xyz[:, 1]
        z         = xyz[:, 2]
        b         = self.b_factor      
  
        if append == True:
            all_lines = open(file_name, 'r').readlines()
            for i in range(len(all_lines)-1, 0, -1):
                if "ATOM" or "HETATM" in all_lines[i]: 

                    final_atom = int(all_lines[i][6:11])
                    break
                
        else:
            all_lines = ''
            final_atom = 0
            
    
    
        f = open(file_name, 'w')
        f.writelines(all_lines)
            
        
        
        
    
        for i in range(len(x)):
    
            if 'C' in atom_type[i]:
                elem = 'C'
            elif 'O' in atom_type[i]:
                elem = 'O'
            elif 'N' in atom_type[i]:
                elem = 'N'
            elif 'Cu' in atom_type[i]:
                elem = 'Cu'
            elif 'Mg' in atom_type[i]:
                elem = 'Mg'
            elif 'Co' in atom_type[i]:
                elem = 'Co'
            elif 'F' in atom_type[i]:
                elem = 'F'
            elif 'H' in atom_type[i]:
                elem = 'H'
            elif 'S' in atom_type[i]:
                elem = 'S'
            elif 'P' in atom_type[i]:
                elem = 'P'
            elif 'Cl' in atom_type[i]:
                elem = 'Cl'
            elif 'Na' in atom_type[i]:
                elem = 'Na'
            elif 'K' in atom_type[i]:
                elem = 'K'
            else:
                elem = 'D'
            
            if resn[i] in self.three2one().keys():
                atom_start = 'ATOM  '
            else:
                atom_start = 'HETATM'
            
            if len(str(b[i])) > 6:
                bfac = float(str(b[i])[:6])
            else:
                bfac = b[i]
            
            my_line  = atom_start + str(final_atom + atom_num[i]).rjust(5) + '  ' + str(atom_type[i]).ljust(4) + \
                resn[i] + ' ' + chain[i] + str(resi[i]).rjust(4) + '    ' + str(x[i]).rjust(8) + \
                str(y[i]).rjust(8) + str(z[i]).rjust(8) + '  1.00'+ str(bfac).rjust(6) +'      X    ' + elem + '\n'
    
            f.write(my_line)
        f.close()
        
        
        
    def RechainChain(self, start_chain, end_chain):
        '''
        

        Parameters
        ----------
        start_chain : str
            The letter/code of the chain you wish to change.
        end_chain : str
            What you wish to change the name to. 

        Returns
        -------
        changes the chain of the pdb object in memory.

        '''
        
        chain_map = self.chain_map
        idx = np.where(chain_map == start_chain)[0]
       
        
        
        for i in idx:
            
            self.chain_map[i] = end_chain
            self.rwc[i]    = self.rwc [i][:-1] + end_chain
            self.rdwc[i]   = self.rdwc[i][:-1] + end_chain
            
        self.chains = np.unique(self.chain_map)        

        
    def RechainResidue(self, residue_list, end_chain):
        
        idx       = np.array([])
        chain_map = self.chain_map
        
        
        for res in residue_list:
            print(idx)
            print(self._residue_idx_(res))
            idx = np.concatenate((idx, self._residue_idx_(res)))
            
        chain_map[idx] = end_chain
        
        
     
    def SelectChain(self, chain):
        
        idx = np.where(self.chain_map == chain)[0]
        
        if len(idx) <= 1:
            print(f'no chain {chain} found in {self.path} !')
            sys.exit()
        
        r = deepcopy(self)
        
        r.atom_type   = self.atom_type[idx]
        r.xyz         = self.xyz[idx, :]
        r.total_atoms = r.xyz.shape[0]
        r.atom_map    = self.atom_map[idx]   
        r.atom_type   = self.atom_type[idx]
        r.resi_map    = self.resi_map[idx]    
        #r.pdb_resi    = self.pdb_resi[idx]    
        #r.ros_resi    = self.ros_resi[idx]    
        r.ros_map     = self.ros_map[idx]     
        r.chain_map   = self.chain_map[idx]   
        r.chains      = np.unique(r.chain_map)      
        r.resn_map    = self.resn_map[idx]    
        r.protein_len = len(np.unique(r.ros_map))  # TODO make all logic match like this!
        r.rwc         = self.rwc[idx]      
        r.rdwc        = self.rdwc[idx]    
    
            
        r.uniq_resi   = np.unique(r.ros_map)
        i = np.where(self.chains == chain)[0][0]
        r.seq         = self.seq[i]
        r.pdb_resi    = np.unique(r.resi_map)
        r.full_path   = self.full_path[:-4] + f'_chain_{chain}.pdb'
        r.path        = self.path[:-4] + f'_chain_{chain}.pdb'
        r.ros_resi    = np.arange(1, len(r.pdb_resi)+1, 1)
        
        
        return r
    

        
        
        
            
            
            
        

