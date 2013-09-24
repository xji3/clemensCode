from protein import *
from data_test import *

url = 'http://www.uniprot.org/uniprot/?query=organism:9606+AND+database:pdb&format=tab&compress=yes&columns=id,database(PDB)'
local_file = 'uniprot-organism%3A9606+AND+database%3Apdb.tab.gz'
if not os.path.isfile(local_file):
    self.downloadFile(url, local_file)
f = gzip.open(local_file, 'rb')
for line in f:
    if line.strip(): # makes sure empty lines are not included
        tmp = line.strip().split()
        print tmp
        print tmp[1]
        up_id = tmp[0]
        plist = tmp[1].split(';')
        print plist
        plist.pop(-1)
        print plist
        a=raw_input()

f.close()
