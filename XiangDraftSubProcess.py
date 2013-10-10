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

naccess_localfile_dir='/Users/xji3/Documents/Naccess'

pdb_file='/Users/xji3/Documents/Naccess/4A14.pdb'
pdb_id='4A14'

if os.path.isfile(naccess_localfile_dir+'/naccess'):
    naccess=naccess_localfile_dir+'/naccess'
else:
    naccess= 'naccess'

print subprocess.check_output([
    'ls',
    '-l'
    ],cwd=naccess_localfile_dir+'/')

if os.path.isfile(naccess_localfile_dir+'/'+pdb_id+'.pdb.gz'):
    print subprocess.check_output([
        'gzip',
        '-d',
        '-f',
        pdb_id+'.pdb.gz',
        ],cwd=naccess_localfile_dir+'/')
    
    print subprocess.check_output([
    naccess,
    pdb_file
    ],cwd=naccess_localfile_dir+'/')
    
elif os.path.isfile(naccess_localfile_dir+'/'+pdb_id+'.pdb'):
    print subprocess.check_output([
    naccess,
    pdb_file
    ],cwd=naccess_localfile_dir+'/')
    
print subprocess.check_output([
    'gzip',
    pdb_id+'.pdb',
    ],cwd=naccess_localfile_dir+'/')






#subprocess.Popen(naccess,cwd=r'./pdb')
##with open('./pdb','w') as outfile:
##    subprocess.check_output([
##        naccess,
##        pdb_file
##        ])
