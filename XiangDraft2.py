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
data.removeNonUniqueChainPDBs()

#---------------------------------------------------------------------------
# Map PDB-IDs to CCDS-IDs (and keep track of all other mappings)
#---------------------------------------------------------------------------
#data.mapPDB2CCDS()


ccds_table = defaultdict(list)
pdb2ccds = defaultdict(list)

tmp_ccds = []
local_file = 'CCDS.current_short.txt'
f = open(local_file, 'rb')
for line in f:
    if line.strip(): # makes sure empty lines are not included
        s = line.strip().split()
##        print s
##        print len(s)
       
        gene_id = s[3]

        chromosome = s[0]
        ccds_id = s[4]
#        a=raw_input()
        print s[3]
        print s[0]
        print s[4]

##        print chromosome
##        print ccds_id

        # Chromosome will be first in the list, followed by the CCDS-IDs
        if not ccds_table[gene_id]:
            ccds_table[gene_id].append(chromosome)

        if not chromosome==ccds_table[gene_id][0]:
            print
            print 'Attention: A gene id shows up in multiple choromosomes!'
            print 'Data.py readCCDStable function'
            print

        ccds_table[gene_id].append(ccds_id)
        tmp_ccds.append(ccds_id)

##        print ccds_table
##        a=raw_input()

print ccds_table.keys()
pdbID_trial=data.pdbIDs

for i in pdbID_trial:
    for gd in ccds_table.keys():
        for ccds in ccds_table[gd]:
            if ccds[0:4] == 'CCDS':
                if not ccds in pdb2ccds[i]:
                    pdb2ccds[i].append(ccds)


for pdb in pdb2ccds:
    # Create a new instance of protein
    proteins[pdb] = Protein(pdb, data.pdb_seqs[pdb][0], data.out_dir + '/pdb/')

    # Assign CCDS
    for ccds in pdb2ccds[pdb]:
        proteins[pdb].ccds_match.append(PDB_CCDS(ccds, ccds_seqs[ccds])) # create a new instance of PDB_CCDS



        
