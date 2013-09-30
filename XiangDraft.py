from protein import *
from data_test import *


##ccds_seqs=defaultdict(list)
##readme=False
##idx=0
##for line in f:
##    if line.strip(): # makes sure empty lines are not included
####        print line
##       line = line.rstrip('\n')
####        print line
####        a=raw_input()
##       if line[0] == '>':
##           line = line.split('|')
##           key = line[0].replace('>', '')
##           readme=True
####           print line
####           print key
##       elif key in ccds_seqs and len(ccds_seqs[key]) == (idx+1):
##           ccds_seqs[key][idx] = ccds_seqs[key][idx] + line
##       else:
##           ccds_seqs[key].append(line)
##    print key
##
##    print ccds_seqs[key]
##    a=raw_input()
       
       
PDB_id='4A14'
#problem pdbIDs
PDB_id='2V4U' #this protein has no match with the structure fragment from Bio.PDB

pdb_ID=PDB_id

pdb_AA='GAMGLPGAEEAPVRVALRVRPLLPKELLHGHQSCLQVEPGLGRVTLGRDRHFGFHVVLAEDAGQEAVYQACVQPLLEAFFEGFNATVFAYGQTGSGKTYTMGEASVASLLEDEQGIVPRAMAEAFKLIDENDLLDCLVHVSYLEVYKEEFRDLLEVGTASRDIQLREDERGNVVLCGVKEVDVEGLDEVLSLLEMGNAARHTGATHLNHLSSRSHTVFTVTLEQRGRAPSRLPRPAPGQLLVSKFHFVDLAGSERVLKTGSTGERLKESIQINSSLLALGNVISALGDPQRRGSHIPYRDSKITRILKDSLGGNAKTVMIACVSPSSSDFDETLNTLNYASRAQ'

ccds_id='CCDS32325.2'

ccds_AA='MGLEAQRLPGAEEAPVRVALRVRPLLPKELLHGHQSCLQVEPGLGRVTLGRDRHFGFHVVLAEDAGQEAVYQACVQPLLEAFFEGFNATVFAYGQTGSGKTYTMGEASVASLLEDEQGIVPRAMAEAFKLIDENDLLDCLVHVSYLEVYKEEFRDLLEVGTASRDIQLREDERGNVVLCGVKEVDVEGLDEVLSLLEMGNAARHTGATHLNHLSSRSHTVFTVTLEQRGRAPSRLPRPAPGQLLVSKFHFVDLAGSERVLKTGSTGERLKESIQINSSLLALGNVISALGDPQRRGSHIPYRDSKITRILKDSLGGNAKTVMIACVSPSSSDFDETLNTLNYASRAQNIRNRATVNWRPEAERPPEETASGARGPPRHRSETRIIHRGRRAPGPATASAAAAMRLGAECARYRACTDAAYSLLRELQAEPGLPGAAARKVRDWLCAVEGERSALSSASGPDSGIESASVEDQAAQGAGGRKEDEGAQQLLTLQNQVARLEEENRDFLAALEDAMEQYKLQSDRLREQQEEMVELRLRLELVRPGWGGPRLLNGLPPGSFVPRPHTAPLGGAHAHVLGMVPPACLPGDEVGSEQRGEQVTNGREAGAELLTEVNRLGSGSSAASEEEEEEEEPPRRTLHLRRNRISNCSQRAGARPGSLPERKGPELCLEELDAAIPGSRAVGGSKARVQARQVPPATASEWRLAQAQQKIRELAINIRMKEELIGELVRTGKAAQALNRQHSQRIRELEQEAEQVRAELSEGQRQLRELEGKELQDAGERSRLQEFRRRVAAAQSQVQVLKEKKQATERLVSLSAQSEKRLQELERNVQLMRQQQGQLQRRLREETEQKRRLEAEMSKRQHRVKELELKHEQQQKILKIKTEEIAAFQRKRRSGSNGSVVSLEQQQKIEEQKKWLDQEMEKVLQQRRALEELGEELHKREAILAKKEALMQEKTGLESKRLRSSQALNEDIVRVSSRLEHLEKELSEKSGQLRQGSAQSQQQIRGEIDSLRQEKDSLLKQRLEIDGKLRQGSLLSPEEERTLFQLDEAIEALDAAIEYKNEAITCRQRVLRASASLLSQCEMNLMAKLSYLSSSETRALLCKYFDKVVTLREEQHQQQIAFSELEMQLEEQQRLVYWLEVALERQRLEMDRQLTLQQKEHEQNMQLLLQQSRDHLGEGLADSRRQYEARIQALEKELGRYMWINQELKQKLGGVNAVGHSRGGEKRSLCSEGRQAPGNEDELHLAPELLWLSPLTEGAPRTREETRDLVHAPLPLTWKRSSLCGEEQGSPEELRQREAAEPLVGRVLPVGEAGLPWNFGPLSKPRRELRRASPGMIDVRKNPL'

water_localfile_dir='/Users/xji3/Downloads/EMBOSS-6.6.0/emboss'
if os.path.isfile(water_localfile_dir+'/water'):
    water = water_localfile_dir+'/water'
else:
    water = 'water'

out = subprocess.check_output([
                        water,
                        '-stdout',
                        '-auto',
                        '-awidth3=100000',
                        '-asequence=asis:' + pdb_AA,
                        '-bsequence=asis:' + ccds_AA
                        ]).split('\n')
#print out

office_Mac_address='/Users/xji3/clemensCode'
address=office_Mac_address
out_dir=address+'/newDataOutput/'
pdb_file = out_dir +'pdb/'+ pdb_ID + '/' + pdb_ID + '.pdb'
handle = gzip.open(pdb_file + '.gz', "r")
structure = Bio.PDB.PDBParser().get_structure(pdb_ID, handle)

print structure
model = structure[0]

chain_list = [ a.get_id() for a in model ]
print chain_list
chain=model[chain_list[0]]

ppb=Bio.PDB.PPBuilder()
print ppb
bb=ppb.build_peptides(chain)
print bb
tmp = [ str(pp.get_sequence()) for pp in ppb.build_peptides(chain, aa_only=False) ]
print tmp

pdb_structure_AA = ''
print not pdb_structure_AA
assert(not pdb_structure_AA)
pdb_structure_AA = '-' * len(pdb_AA)

edges = [None] * len(tmp)
print edges
multiple=[]

a=tmp[0]
i=0
n = pdb_AA.count(a)
if(n == 1):
    idx = string.find(pdb_AA, a) # return the first element in a's notion in pdb_AA, could be 0 @Xiang
    pdb_structure_AA = pdb_structure_AA[:idx] + a + pdb_structure_AA[idx+len(a):]
    edges[i] = [idx, idx+len(a)]
elif n > 1:
    multiple.append(i)
else:
    print
i=i+1

for i in multiple:
    b = e = -1
    if i == 0 and edges[i+1]:
        b = 0
        e = edges[i+1][0]
    elif i == (len(tmp)-1) and edges[i-1]:
        b = edges[i-1][1]
        e = len(tmp)-1
    elif edges[i-1] and edges[i+1]:
        b = edges[i-1][1]
        e = edges[i+1][0]
    else:
        print 'I don\'t know where to place this fragment. It appears >= 2x in the structure.'
        print 'It will be ignored in the distance matrix:', tmp[i]
        continue
    if self.pdb_AA[b:e].count(tmp[i]) == 1:
        idx = string.find(self.pdb_AA[b:e], tmp[i])
        self.pdb_structure_AA = self.pdb_structure_AA[:b+idx] + tmp[i] + self.pdb_structure_AA[b+idx+len(tmp[i]):]
    else:
        print 'I don\'t know where to place this fragment. It appears >= 2x in the structure.',
        print 'It will be ignored in the distance matrix:', tmp[i]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

iter = out.__iter__()
for line in iter:
    if line.strip() and line.split() and line.split()[0] == 'asis':
        break
aligned_subseq1 = line.split()
pdb_local_alignment_start = int(aligned_subseq1[1]) - 1
s1 = pdb_alignedAA = aligned_subseq1[2] # PDB sequence
iter.next() # skip the vertical bars
aligned_subseq2 = iter.next().split()
assert(len(aligned_subseq1[2]) == len(aligned_subseq2[2]))
ccds_local_alignment_start = int(aligned_subseq2[1]) - 1
s2 = ccds_alignedAA = aligned_subseq2[2] # CCDS sequence
algnmt_length = local_alignment_length = len(aligned_subseq2[2])

s3 = pdb_structure_AA[pdb_local_alignment_start:]
gaps = [ i.start() for i in re.finditer('-', s1) ]
print gaps
string.replace(s3, '-', '')


        
