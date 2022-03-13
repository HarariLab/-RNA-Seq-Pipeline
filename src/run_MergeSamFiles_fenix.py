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
import shutil

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

parser = argparse.ArgumentParser(description='Convert FASTQ to UBAM using MergeSamFiles from Picard.')

#parser.add_argument('infile_array',type=str,
#                   help='Path to alignment .bam files to be processed')
parser.add_argument('bams', nargs='+',type=str, help='BAM files')
parser.add_argument('prefix',type=str,default='null', help='prefix for output files')
parser.add_argument('-m', '--memory', default='8',type=str, help='Memory in GB')
parser.add_argument('-o','--output_dir', default=os.getcwd(), help='Output directory')
parser.add_argument('--jar', default='/usr/local/genome/picard_latest/picard.jar', help='Path to Picard.jar')
parser.add_argument('--tmp_dir', default='/data/tmp', help='Directory for tmp files')

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
if not os.path.exists(args.tmp_dir):
    os.makedirs(args.tmp_dir)

with cd(args.output_dir):
    subprocess.check_call('java -jar -Xmx' + args.memory +'g ' + args.jar \
        +' MergeSamFiles I='+' I='.join(args.bams) \
        +' O='+args.prefix+'_merged.bam' \
        +' TMP_DIR='+args.tmp_dir,
            shell=True)

if os.path.exists(args.tmp_dir):
    shutil.rmtree(args.tmp_dir)
