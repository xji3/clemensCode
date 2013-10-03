from protein import *
from data_test import *
import tempfile
import os


##water_localfile_dir='/Users/xji3/Downloads/EMBOSS-6.6.0/emboss'
##        if os.path.isfile(water_localfile_dir+'/water'):
##            water = water_localfile_dir+'/water'
##        else:
##            water = 'water'
##        return subprocess.check_output([
##                        water,
##                        '-stdout',
##                        '-auto',
##                        '-awidth3=100000',
##                        '-asequence=asis:' + self.pdb_AA,
##                        '-bsequence=asis:' + ccds_AA
##                    ]).split('\n')

msms_localfile_dir='/Users/xji3/clemensCode/aux'
if os.path.isfile(msms_localfile_dir+'/msms'):
    MSMS = msms_localfile_dir+'/msms'
else:
    MSMS = 'msms'


PDB_TO_XYZR="./aux/pdb_to_xyzr"
#pdb_file='/Users/xji3/clemensCode/newDataOutput/pdb/1CS8/1CS8.pdb.gz'
pdb_file='/Users/xji3/clemensCode/tmp/4A14.pdb'

#xyz_tmp=tempfile.mktemp()
if not '/Users/xji3/clemensCode/tmp/':
    tmp_dir=os.mkdir('/Users/xji3/clemensCode/tmp/')
else:
    tmp_dir='/Users/xji3/clemensCode/tmp/'
xyz_tmp=tempfile.mktemp(dir=tmp_dir)+'.xyzr'
PDB_TO_XYZR=PDB_TO_XYZR+" %s > %s"
make_xyz=PDB_TO_XYZR % (pdb_file, xyz_tmp)
if os.system(make_xyz):
    terminate=True
else:
    os.system(make_xyz)

##subprocess.check_output([
##    PDB_TO_XYZR,
##    pdb_file,
##    'a.xyzr'
##    ])

surface_tmp=tempfile.mktemp(dir=tmp_dir)
MSMS=MSMS+" -probe_radius 1.5 -if %s -of %s > "+tempfile.mktemp(dir=tmp_dir)
make_surface=MSMS % (xyz_tmp, surface_tmp)
print make_surface

if os.system(make_surface):
    terminate=True
else:
    os.system(make_surface)

surface_file=surface_tmp+'.vert'

fp=open(surface_file, "rb")
vertex_list=[]
for l in fp.readlines():
    sl=l.split()
#cl: length is 10
#        if not len(sl)==9:
    if not len(sl)==10:
        # skip header
        continue
    vl=map(float, sl[0:3])
    vertex_list.append(vl)
fp.close()

