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

pdb_ID='4A14'
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

a=tmp[0]
i=0
n = pdb_AA.count(a)
if(n == 1):
    idx = string.find(pdb_AA, a) # return the first element in a's notion in pdb_AA, could be 0 @Xiang
    pdb_structure_AA = pdb_structure_AA[:idx] + a + pdb_structure_AA[idx+len(a):]
    edges[i] = [idx, idx+len(a)]
