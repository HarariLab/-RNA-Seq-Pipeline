#!/usr/bin/env python3
import argparse
#import simplejson as json
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
import glob

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

parser = argparse.ArgumentParser(description='RNA-seq pipeline from unaligned BAM')
parser.add_argument('host', type=str, help='host server with raw data')
parser.add_argument('src', type=str, help='Path to src directory with pipeline scripts')
parser.add_argument('seq_file', type=str, help='Path to local bam in host')
parser.add_argument('file_type', type=str, help='file type [bam | fastq]')
parser.add_argument('cohort', type=str, help='Cohort of the sample ID')
parser.add_argument('tissues',  type=str, help='Tissue from a cohort to process. Format: tissue1 [tissue2], or comma-separated list') #nargs='+',
parser.add_argument('sample_id', type=str, help='Sample ID for bam file')
parser.add_argument('ref', type=str, help='Reference genome fasta')
parser.add_argument('star_index', type=str, help='Path to directory storing star index files')
parser.add_argument('salmon_index', type=str, help='Path to transcripts fasta file')
parser.add_argument('ref_flat', type=str, help='Path to reference annotation in flat format')
parser.add_argument('ribosomal_int', type=str, help='Ribosomal (rRNA) interval list ')
parser.add_argument('annot', type=str, help='Path to annotation file used for STAR, Salmon and DCC')
parser.add_argument('annot_bed', type=str, help='Path to annotation file in BED foramt for TIN')
parser.add_argument('repeat_mask', type=str, help='Path to repeat masking file for DCC')
parser.add_argument('read_len', type=str, help='Read length used for RNA-seq')
parser.add_argument('read_type', type=str, help='Read type (SE or PE)')
parser.add_argument('adapter_seq', type=str, default=None, help='Adapter sequences')
parser.add_argument('tmp_dir', type=str, help='Path to directory for storing tmp files')
parser.add_argument('output_dir', type=str, help='Path to output directory for all processes')
parser.add_argument('storage_aligned', type=str, help='Path to output directory for aligned files')
parser.add_argument('storage_bams', type=str, help='Path to output directory for marked duplicate bams')
parser.add_argument('storage_qc', type=str, help='Path to output directory for qc files')
parser.add_argument('storage_quant', type=str, help='Path to output directory for quant tables')
parser.add_argument('storage_tin', type=str, help='Path to output directory for TIN files')
parser.add_argument('storage_circ', type=str, help='Path to output directory for circRNA alignments')
parser.add_argument('star_storage_box', type=str, help='Path to the box location for STAR crammed files')
parser.add_argument('circ_storage_box', type=str, help='Path to the box location for circRNA crammed files')
#parser.add_argument('json', type=str, help='Path to the JSON file')
parser.add_argument('cram_move_delete', type=str, help='cram/move/delete bam files')
parser.add_argument('run_circ', type=str, help='run circRNA alignment')
parser.add_argument('remote_storage_name', type=str, help='remote storage name')

args = parser.parse_args()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

#################################################
## 9. cram bam files and move the mto the box
#################################################
if (args.cram_move_delete.lower()=="yes"):

     # 9.0 create directories on the box
     cmd='rclone mkdir '+args.remote_storage_name+':'+args.star_storage_box
     logger.info(' cmd * '+cmd)
     subprocess.check_call(cmd,shell=True)

     if (args.run_circ=="run_circ"):
        cmd='rclone mkdir '+args.remote_storage_name+':'+args.circ_storage_box
        logger.info(' cmd * '+cmd)
        subprocess.check_call(cmd,shell=True)

     # 9. Post-processing: Craming and Delete unnecessary bams
     # 9.1 cram aligned bam files
     cmd='samtools view -C '+os.path.join(args.storage_bams,args.sample_id)+'.Aligned.sortedByCoord.out.md.bam -T '+args.ref \
     +' -o '+os.path.join(args.storage_bams,args.sample_id)+'.Aligned.sortedByCoord.out.md.cram'
     logger.info(' cmd * '+cmd)
     subprocess.check_call(cmd,shell=True)

     # 9.2 cram circRNA aligned bam file
     #if (args.run_circ=="run_circ"):
     #   cmd='samtools view -C '+os.path.join(args.storage_circ,args.sample_id,args.sample_id)+'_unified.Aligned.sortedByCoord.out.bam -T '+args.ref \
     #   +' -o '+os.path.join(args.storage_circ,args.sample_id,args.sample_id) +'_unified.Aligned.sortedByCoord.out.cram'
     #   logger.info(' cmd * '+cmd)
     #   subprocess.check_call(cmd,shell=True)

     # 9.3 move crams to the box
     cmd='rclone copy '+os.path.join(args.storage_bams,args.sample_id)+'.Aligned.sortedByCoord.out.md.cram '+args.remote_storage_name+':'+args.star_storage_box
     logger.info(' cmd * '+cmd)
     print(cmd)
     subprocess.check_call(cmd,shell=True)

     # move circRNA cram
     #if (args.run_circ=="run_circ"):
     #   cmd='rclone copy '+os.path.join(args.storage_circ,args.sample_id,args.sample_id) +'_unified.Aligned.sortedByCoord.out.cram '+args.remote_storage_name+':'+args.circ_storage_box
     #   logger.info(' cmd * '+cmd)
     #   print(cmd)
     #   subprocess.check_call(cmd,shell=True)

     # delete aligned data
     #os.remove(os.path.join(args.storage_aligned,args.sample_id,args.sample_id+'.Aligned.sortedByCoord.out.bam'))
     #os.remove(os.path.join(args.storage_aligned,args.sample_id,args.sample_id+'.Aligned.sortedByCoord.out.bam.bai'))
     #os.remove(os.path.join(args.storage_aligned,args.sample_id,args.sample_id+'.Aligned.toTranscriptome.out.bam'))
     #os.remove(os.path.join(args.storage_bams,args.sample_id+'.Aligned.sortedByCoord.out.md.bam'))
     #os.remove(os.path.join(args.storage_bams,args.sample_id+'.Aligned.sortedByCoord.out.md.bam.bai'))

     #os.remove(os.path.join(args.storage_bams,args.sample_id+'.Aligned.sortedByCoord.out.md.cram'))
     #os.remove(os.path.join(args.storage_circ,args.sample_id,args.sample_id +'_unified.Aligned.sortedByCoord.out.cram'))
