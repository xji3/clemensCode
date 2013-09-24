from protein import *
from data_test import *

local_file = 'CCDS_protein.current.faa.gz'
f = gzip.open(local_file, 'rb')

ccds_seqs=defaultdict(list)
readme=False
idx=0
for line in f:
    if line.strip(): # makes sure empty lines are not included
##        print line
       line = line.rstrip('\n')
##        print line
##        a=raw_input()
       if line[0] == '>':
           line = line.split('|')
           key = line[0].replace('>', '')
           readme=True
##           print line
##           print key
       elif key in ccds_seqs and len(ccds_seqs[key]) == (idx+1):
           ccds_seqs[key][idx] = ccds_seqs[key][idx] + line
       else:
           ccds_seqs[key].append(line)
    print key

    print ccds_seqs[key]
    a=raw_input()
        
       

