#!/usr/bin/env python3

# Salmon Quantification - Alignment-based Mode
# Author: James A. Perez

import argparse
import os
import subprocess
import gzip
import shutil
from datetime import datetime
import contextlib
import time
import logging

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description='Salmon quantification in alignment-based mode')
parser.add_argument('aligned_bam', type=str, help='Path to aligned .bam file')
parser.add_argument('transcripts', type=str, help='Path to set of transcripts [*.fa] to quantify')
parser.add_argument('annot', type=str, help='Path to gene annotations for transcripts file')
parser.add_argument('-l', '--lib_type', type=str, default='A',help='library type specification for reads')
parser.add_argument('-o', '--output_dir', type=str, default='.',help='Head output directory to processed files')
parser.add_argument('-p', '--software_path', default=None, help='Path to Salmon executable')
parser.add_argument('--threads', default='12', help='Number of threads')
parser.add_argument('-z','--debug', type=bool, default=False, help='Output cmd to file instead of running')

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

if args.software_path:
    cmd = args.software_path+' quant'
else:
    #cmd = '/opt/salmon-latest_linux_x86_64/bin/salmon quant' #MGI
    cmd = '/usr/local/genome/salmon-0.11.3/bin/salmon quant'

cmd+=' --threads '+args.threads\
    +' -t '+args.transcripts \
    +' -a '+args.aligned_bam+' '\
    +' --libType '+args.lib_type \
    +' --seqBias --gcBias'\
    +' --geneMap '+args.annot \
    +' -o '+args.output_dir

if args.debug==True:
    with open(os.path.join(args.output_dir,'Salmon_cmd.txt'), 'w') as f:
        f.write(cmd)
        quit()

subprocess.check_call(cmd, shell=True, executable='/bin/bash')
