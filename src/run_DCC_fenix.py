#!/usr/bin/env python3
# run_DCC_fenix.py

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
#from pathlib import Path

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

parser = argparse.ArgumentParser(description='Run pipeline from DCC JSON file')
parser.add_argument('json', type=str, help='Path to the JSON file')
parser.add_argument('study', type=str, help='Studies from which to align and quantify reads from.')
parser.add_argument('tissue',  type=str, help='Tissue to process.')
parser.add_argument('end', type=str.lower, choices=['pe','se'], help='paired end (PE) or single end (SE)')
parser.add_argument('src', type=str, help='Path to src directory with pipeline scripts')
parser.add_argument('storage_aligned', type=str, help='Path to output directory for aligned files')
parser.add_argument('storage_circ', type=str, help='Path to output directory for circRNA alignments')
parser.add_argument('storage_dcc', type=str, help='Path to output directory for DCC results')
parser.add_argument('storage_dcc_in', type=str, help='Path to directory for DCC input files')
parser.add_argument('repeat_mask', type=str, help='Path to repeat masking file for DCC')
parser.add_argument('annot', type=str, help='Path to annotation file used for STAR, Salmon and DCC')
parser.add_argument('ref', type=str, help='Reference genome fasta')
#parser = argparse.ArgumentParser(description='DCC to call circRNAs')
#parser.add_argument('cohort', type=str, help='Cohort of the sample ID')
#parser.add_argument('seq_file', type=str, help='Path to local bam in host')
#parser.add_argument('sample_id', type=str, help='Sample ID for bam file')
#parser.add_argument('output_dir', type=str, help='Path to output directory for all processes')



args = parser.parse_args()


# Setup Logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

timestamp = datetime.now().strftime('%Y-%m-%d.%H:%M:%S')

if not os.path.exists(args.study+"_logs"):
    os.makedirs(args.study+"_logs")

log_file = os.path.join(args.study+"_logs",args.study+'_log.txt')

fh = logging.FileHandler(log_file, mode='w')
fh.setLevel(logging.INFO)
logging.Formatter.converter = customTime
display = logging.Formatter("%(asctime)s - %(name)s - %(message)s")
fh.setFormatter(display)
logger.addHandler(fh)

output_dir = os.getcwd()
storage_dcc = os.path.join(output_dir,args.storage_dcc,args.tissue)
storage_dcc_in = os.path.join(output_dir,args.storage_dcc_in,args.tissue)
#if not os.path.exists(args.output_dir):
#    os.makedirs(args.output_dir)

if not os.path.exists(storage_dcc):
    os.makedirs(storage_dcc)

if not os.path.exists(storage_dcc_in):
    os.makedirs(storage_dcc_in)

print('output destinations: '+'\n\t'+storage_dcc+'\n\t'+storage_dcc_in)

# Get the list of folders = samples
#p = Path('.')
#print(p)
#sample_id = [x for x in p.iterdir() if x.is_dir()]
#print(sample_id)
samples = os.listdir(path=os.path.join(args.storage_circ,args.tissue))
samples.sort()
#print(samples)
#print(type(samples))
#print(samples.sort())
# Make input files for PE and SE:
#logger.info('DCC input files: unified')
for sample_id in samples:
    if sample_id =='DCC':
        continue
    cmd = 'ls -d '+ os.path.join(args.storage_circ,args.tissue)+'/'+sample_id+'/'+sample_id+'_unified.Chimeric.out.junction >> ' + os.path.join(storage_dcc_in)+'/samplesheet'
    logger.info(' cmd * '+cmd)
    print('sample_id: ' +sample_id)
    print(cmd)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    logger.info('DCC input files: bamfiles')
    cmd = 'ls -d '+ os.path.join(args.storage_aligned,sample_id)+'.Aligned.sortedByCoord.out.md_fixed.bam >> ' + os.path.join(storage_dcc_in)+'/bam_files'

    logger.info(' cmd * '+cmd)
    print(cmd)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

# PE only
if (args.end == "pe"):
    for sample_id in samples:
        #print(sample_id)
        if sample_id =='DCC':
            continue
        logger.info('DCC input files: R1')
        cmd = 'ls -d '+ os.path.join(args.storage_circ,args.tissue)+'/'+sample_id+'/'+sample_id+'_R1.Chimeric.out.junction >> ' + os.path.join(storage_dcc_in)+'/mate1'

        logger.info(' cmd * '+cmd)
        print(cmd)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')

        logger.info('DCC input files: R2')
        cmd = 'ls -d '+ os.path.join(args.storage_circ,args.tissue)+'/'+sample_id+'/'+sample_id+'_R2.Chimeric.out.junction >> ' + os.path.join(storage_dcc_in)+'/mate2'

        logger.info(' cmd * '+cmd)
        print(cmd)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    # Run DCC
    logger.info('Starting DCC for: '+ args.study)
    os.chdir(storage_dcc)
    cmd = 'DCC @' +os.path.join(storage_dcc_in)+'/samplesheet -mt1 @' +os.path.join(storage_dcc_in)+'/mate1 -mt2 @'+os.path.join(storage_dcc_in)+'/mate2 -T 20 -D -R ' +args.repeat_mask +' -an '+args.annot +' -Pi -F -M -Nr 1 1 -fg -k -G -A ' +args.ref +' -B @'+os.path.join(storage_dcc_in)+'/bam_files'

    logger.info(' cmd * '+cmd)
    print(cmd)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')



# SE only
else:
    samples = os.listdir(path=os.path.join(args.storage_circ,args.tissue))
    samples.sort()
    print(samples)
    for sample_id in samples:
        print(sample_id)
        if sample_id =='DCC':
            continue
        logger.info('DCC input files: R1')
        cmd = 'ls -d '+ os.path.join(args.storage_circ,args.tissue)+'/'+sample_id+'/'+sample_id+'_unified.Chimeric.out.junction >> ' + os.path.join(storage_dcc_in)+'/read'

        logger.info(' cmd * '+cmd)
        print(cmd)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    # Run DCC
    logger.info('Starting DCC for: '+ args.study)
    os.chdir(storage_dcc)
    cmd = 'DCC @' +os.path.join(storage_dcc_in)+'/samplesheet -mt1 @' +os.path.join(storage_dcc_in)+'/read -T 10 -D -N -R ' +args.repeat_mask +' -an '+args.annot +' -F -M -Nr 1 1 -fg -k -G -A ' +args.ref +' -B @'+os.path.join(storage_dcc_in)+'/bam_files'

    logger.info(' cmd * '+cmd)
    print(cmd)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

# Rename outfiles
out = os.path.join(args.study+'_'+args.tissue)
outa = os.path.join(out+'_CircCoordinates.txt')
outb = os.path.join(out+'_CircRNACount.txt')
outc = os.path.join(out+'_CircSkipJunctions.txt')
outd = os.path.join(out+'_LinearCount.txt')
print(outa)
os.rename('CircCoordinates', outa)
os.rename('CircRNACount', outb)
os.rename('CircSkipJunctions', outc)
os.rename('LinearCount', outd)

# Remove temp folder
outt = os.getcwd()
outf = os.path.join(outt+'/_tmp_DCC')
print(outf)
shutil.rmtree(outf)
