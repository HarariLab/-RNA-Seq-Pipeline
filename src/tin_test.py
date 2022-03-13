#!/usr/bin/env python3
# tin_test.py

import argparse
import simplejson as json
import subprocess
import contextlib
import shutil
import os
import sys
import logging
from datetime import datetime
import time
import shlex
from pytz import timezone, utc

def customTime(*args):
    utc_dt = utc.localize(datetime.utcnow())
    ctz = timezone("US/Central")
    ctz_dt = utc_dt.astimezone(ctz)
    return(ctz_dt.timetuple())
def get_size(path, host=None):
    if isinstance(path, list):
        path=' '.join(path)
    if host:
        cnx = subprocess.Popen(['ssh', host, 'du -sh', path], stdout=subprocess.PIPE)
    else:
        cmd='du -sh '+path
        cnx = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
    return(cnx.communicate()[0].decode('utf-8'))
@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

parser = argparse.ArgumentParser(description='test of tin.py')
parser.add_argument('src', type=str, help='Path to src directory with pipeline scripts')
parser.add_argument('seq_file', type=str, help='Path to local bam in host')
parser.add_argument('cohort', type=str, help='Cohort of the sample ID')
parser.add_argument('sample_id', type=str, help='Sample ID for bam file')
parser.add_argument('output_dir', type=str, help='Path to output directory for all processes')
parser.add_argument('storage_aligned', type=str, help='Path to output directory for aligned files')
parser.add_argument('storage_tin', type=str, help='Path to output directory for tin files')



args = parser.parse_args()

print('sample: '+args.sample_id)
print('seq_file: '+args.seq_file)
print('output destinations:\n'+'\n\t'+args.storage_tin)

# Setup Logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

timestamp = datetime.now().strftime('%Y-%m-%d.%H:%M:%S')

if not os.path.exists(args.cohort+"_logs"):
    os.makedirs(args.cohort+"_logs")

log_file = os.path.join(args.cohort+"_logs",args.sample_id+'_log.txt')

fh = logging.FileHandler(log_file, mode='w')
fh.setLevel(logging.INFO)
logging.Formatter.converter = customTime
display = logging.Formatter("%(asctime)s - %(name)s - %(message)s")
fh.setFormatter(display)
logger.addHandler(fh)

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# 6. Run TIN
logger.info('Running TIN')
cmd = os.path.join(args.src,'tin.py -i') \
    +os.path.join(args.storage_aligned,args.sample_id \
        +'.Aligned.sortedByCoord.out.md.bam ') \
    +' -r /40/pipelines/RNAseq/Dockerfiles/unified_rnaseq_v1/hg19_GencodeCompV19.bed'

logger.info(' cmd * '+cmd)
subprocess.check_call(cmd, shell=True)

logger.info('Finished TIN for: '+args.sample_id)
