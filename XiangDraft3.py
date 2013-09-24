from protein import *
from data_test import *

office_Mac_address='/Users/xji3/Dropbox/My Files/BRC/Small Project/Clemens PY code'
initpdb_filename=office_Mac_address+'/input/init_pdbids.txt'
xmlfilename=office_Mac_address+'/input/pdb_search.xml'
ds_dir=office_Mac_address+'/XiangTrialOutput'
initpdb_filename=office_Mac_address+'/input/init_pdbids.txt'
initmono_filename=office_Mac_address+'/input/monomers.txt'
local_pisadir=office_Mac_address+'/XiangTrialOutput/pisa'
initmp_filename=office_Mac_address+'/input/init_mpids.txt'

data = Data(initpdb_filename, xmlfilename, ds_dir)
data.queryPDB(xmlfilename)

    
#---------------------------------------------------------------------------
# Apply additional filters (no membrane proteins, only monomers, ...)
#---------------------------------------------------------------------------
data.filterMembraneProteins(initmp_filename)
data.filterNonMonomers(initmono_filename, local_pisadir)

#---------------------------------------------------------------------------
# Get the AA sequences for the PDB files and remove those with non-unique
# chains
#---------------------------------------------------------------------------
data.fasta4pdb()

remove_id = []
#print 'PDB IDs with more than one chain:\n\t',
n_mtoc = 0
for i in data.pdb_seqs:
    print i
    if len(data.pdb_seqs[i]) > 1:
        n_mtoc += 1
        #print i + ' (' + str(len(self.pdb_seqs[i])) + ')',
    data.pdb_seqs[i] = list(set(data.pdb_seqs[i]))
    if len(data.pdb_seqs[i]) > 1:
        remove_id.append(i)
print '\nTotal:', n_mtoc, '(of', str(len(data.pdbIDs)) + ' remaining PDB IDs)'
        
