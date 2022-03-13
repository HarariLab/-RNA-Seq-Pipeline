#!/usr/bin/env python3

"""
Spyder Editor

This is a temporary script file.
"""

import argparse
import os
import struct
import subprocess
from datetime import datetime
import contextlib

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

parser = argparse.ArgumentParser(description='Convert FASTQ to UBAM using FastqToSam from Picard.')

#parser.add_argument('infile_array',type=str,
#                   help='Path to alignment .bam files to be processed')
parser.add_argument('fastq', nargs='+', help='FASTQ file')
parser.add_argument('sample_id', type=str,default='null')
parser.add_argument('-t','--read_type', type=str.lower, choices=['se','pe'], default='pe', help='paired-end (pe) OR single-end (se)')
parser.add_argument('-o', '--output_dir', default=os.getcwd(), help='Output directory')
parser.add_argument('-m', '--memory', default='8',type=str, help='Memory in GB')
parser.add_argument('--jar', default='/usr/local/genome/picard_latest/picard.jar', help='Path to Picard.jar')
parser.add_argument('--tmp_dir', default='~/tmp', help='Directory for tmp files')

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
if not os.path.exists(args.tmp_dir):
    os.makedirs(args.tmp_dir)

with cd(args.output_dir):

    cmd='java -jar -Xmx' + args.memory +'g ' + args.jar\
        + ' FastqToSam'
    if args.read_type=='se':
        cmd+=' F1='+args.fastq[0]
    else:
        cmd+=' F1='+args.fastq[0]\
           +' F2='+args.fastq[1]
    cmd+=' O='+args.sample_id+'_unmapped.bam'\
       +' SM='+args.sample_id \
       +' TMP_DIR='+args.tmp_dir
    subprocess.check_call(cmd, shell=True)
