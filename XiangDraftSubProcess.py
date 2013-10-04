#import itertools, numpy, string, subprocess, re
#from __future__ import division
import subprocess
from collections import defaultdict
import os.path, gzip
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

pdb_file='/Users/xji3/Documents/Naccess/4A14.pdb.gz'

naccess_localfile_dir='/Users/xji3/Documents/Naccess'
if os.path.isfile(naccess_localfile_dir+'/naccess'):
    naccess=naccess_localfile_dir+'/naccess'
else:
    naccess= 'naccess'

print subprocess.check_output([
    naccess,
    pdb_file
    ],cwd=naccess_localfile_dir+'/')

#subprocess.Popen(naccess,cwd=r'./pdb')
##with open('./pdb','w') as outfile:
##    subprocess.check_output([
##        naccess,
##        pdb_file
##        ])
