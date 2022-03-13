#!/usr/bin/env python3
# run_DCC.py

import argparse
#import simplejson as json
import json
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

args = parser.parse_args()

# Pull information from json
with open(args.json) as data:
    paths = json.load(data)

genome = paths.get('genome'); #print('genome: ', end=''); print(genome)
ref = genome.get('ref')
annot = genome.get('annot')
repeat_mask = genome.get('repeat_mask')
#print('annot: repeat_mask: ', end=''); print(annot, end="\t"); print(repeat_mask) # NOTE: How to change the newline default for print to '' or tab
storage = paths.get('storage'); #print('storage: ', end=''); print(storage)
bams = storage.get('bams'); print('bams: ', end=''); print(bams)
outdir = paths.get('outdir'); print('outdir: ', end=''); print(outdir)
storage_bams = os.path.join(outdir, args.study, bams, args.tissue); print('storage_bams: ', end=''); print(storage_bams)

circ = storage.get('circ'); print('circ: ', end=''); print(circ)
storage_circ = os.path.join(outdir, args.study, circ, args.tissue); print('storage_circ: ', end=''); print(storage_circ)

cohorts = paths.get('cohorts')

if args.study not in [s for s in cohorts.keys()]:
   raise ValueError(study+' not in json file. Found: '+', '.join([s for s in cohorts.keys()]))

for cohort, entities in cohorts.items():
    if cohort != args.study:
        continue

    tissues = entities.get('tissues')
    #print('tissues: ', end=''); print(tissues)
    if args.tissue not in [t for t in tissues.keys()]:
        raise ValueError(args.tissue+' not in json file. Found: '+', '.join([t for t in tissues.keys()]))

tokens = tissues.get(args.tissue); #print('tokens: ', end=''); print(tokens)
end = tokens.get('read_type'); #print('end: ', end=''); print(end)

#print('server host: '+paths.get('host'))
local_genome_dir = os.path.join(paths.get('resources'),'genome')
print(local_genome_dir)

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

dcc = storage.get('dcc'); print('dcc: ', end=''); print(dcc)
storage_dcc = os.path.join(outdir, args.study, dcc, args.tissue); print('storage_dcc: ', end=''); print(storage_dcc)

dcc_in = storage.get('dcc_in'); print('dcc_in: ', end=''); print(dcc_in)
storage_dcc_in = os.path.join(outdir, args.study, dcc_in, args.tissue); print('storage_dcc_in: ', end=''); print(storage_dcc_in)

if not os.path.exists(storage_dcc): # TURN FOR MGI
    os.makedirs(storage_dcc) # TURN FOR MGI

if not os.path.exists(storage_dcc_in): # TURN FOR MGI
    os.makedirs(storage_dcc_in) # TURN FOR MGI

print('output destinations: '+'\n\t'+storage_dcc+'\n\t'+storage_dcc_in)

# Get the list of folders = samples
samples = os.listdir(path=os.path.join(storage_circ)) # TURN ON FOR MGI
#tmp = "/40/tmp/Miguel/02.-ProcessedData/06.-circRNA/Hg19" # fenix hack
#samples = os.listdir(path=os.path.join(tmp,args.tissue)) # fenix hack
samples.sort()
#print(samples)

# Make input files for PE and SE:
logger.info('DCC input files: unified')
for sample_id in samples:
    if sample_id =='DCC':
        continue
    cmd = 'ls -d '+ os.path.join(storage_circ)+'/'+sample_id+'/'+sample_id+'_unified.Chimeric.out.junction >> ' + os.path.join(storage_dcc_in)+'/samplesheet'
    logger.info(' cmd * '+cmd)
    print('sample_id: ' +sample_id)
    print(cmd)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    logger.info('DCC input files: bamfiles')
    cmd = 'ls -d '+ os.path.join(storage_bams,sample_id)+'.Aligned.sortedByCoord.out.md_fixed.bam >> ' + os.path.join(storage_dcc_in)+'/bam_files'

    logger.info(' cmd * '+cmd)
    print(cmd)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

# PE only
if (end == "PE"):
    for sample_id in samples:
        #print(sample_id)
        if sample_id =='DCC':
            continue
        logger.info('DCC input files: R1')
        cmd = 'ls -d '+ os.path.join(storage_circ)+'/'+sample_id+'/'+sample_id+'_R1.Chimeric.out.junction >> ' + os.path.join(storage_dcc_in)+'/mate1'

        logger.info(' cmd * '+cmd)
        print(cmd)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')

        logger.info('DCC input files: R2')
        cmd = 'ls -d '+ os.path.join(storage_circ)+'/'+sample_id+'/'+sample_id+'_R2.Chimeric.out.junction >> ' + os.path.join(storage_dcc_in)+'/mate2'

        logger.info(' cmd * '+cmd)
        print(cmd)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    # Run DCC
    logger.info('Starting DCC for: '+ args.study)
    os.chdir(storage_dcc) #TURN ON FOR MGI
    cmd = 'DCC @' +os.path.join(storage_dcc_in)+'/samplesheet -mt1 @' +os.path.join(storage_dcc_in)+'/mate1 -mt2 @'+os.path.join(storage_dcc_in)+'/mate2 -T 20 -D -R ' +repeat_mask +' -an '+annot +' -Pi -F -M -Nr 1 1 -fg -k -G -A ' +ref+' -B @'+os.path.join(storage_dcc_in)+'/bam_files'

    logger.info(' cmd * '+cmd)
    print(cmd)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')



# SE only
else:
    samples = os.listdir(path=os.path.join(storage_circ)) #TURN ON FOR MGI
    #tmp = "/40/tmp/Miguel/02.-ProcessedData/06.-circRNA/Hg19" # fenix hack
    #samples = os.listdir(path=os.path.join(tmp,args.tissue)) # fenix hack
    samples.sort()
    print(samples)
    for sample_id in samples:
        print(sample_id)
        if sample_id =='DCC':
            continue
        logger.info('DCC input files: R1')
        cmd = 'ls -d '+ os.path.join(storage_circ)+'/'+sample_id+'/'+sample_id+'_unified.Chimeric.out.junction >> ' + os.path.join(storage_dcc_in)+'/read'

        logger.info(' cmd * '+cmd)
        print(cmd)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    # Run DCC
    logger.info('Starting DCC for: '+ args.study)
    os.chdir(storage_dcc) #TURN ON FOR MGI
    cmd = 'DCC @' +os.path.join(storage_dcc_in)+'/samplesheet -mt1 @' +os.path.join(storage_dcc_in)+'/read -T 10 -D -N -R ' +repeat_mask +' -an '+annot +' -F -M -Nr 1 1 -fg -k -G -A ' +ref+' -B @'+os.path.join(storage_dcc_in)+'/bam_files'

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
outf = os.path.join(storage_dcc+'/_tmp_DCC')
print(outf)
shutil.rmtree(outf)
