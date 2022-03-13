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
parser.add_argument('tmp_dir', type=str, help='Path to directory for storing tmp files')
parser.add_argument('output_dir', type=str, help='Path to output directory for all processes')
parser.add_argument('storage_aligned', type=str, help='Path to output directory for aligned files')
parser.add_argument('storage_bams', type=str, help='Path to output directory for marked duplicate bams')
parser.add_argument('storage_qc', type=str, help='Path to output directory for qc files')
parser.add_argument('storage_quant', type=str, help='Path to output directory for quant tables')
parser.add_argument('storage_tin', type=str, help='Path to output directory for TIN files')
parser.add_argument('storage_circ', type=str, help='Path to output directory for circRNA alignments')
parser.add_argument('json', type=str, help='Path to the JSON file')

args = parser.parse_args()

print('sample: '+args.sample_id)
print('seq_file: '+args.seq_file)

# Setup Logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# 7. Run circRNA
logger.info('Starting circRNA for: '+ args.sample_id)
#print(args.tissues)
if not os.path.exists(os.path.join(args.storage_circ,args.sample_id)): #,args.tissues
    os.makedirs(os.path.join(args.storage_circ,args.sample_id)) #,args.tissues

# extract the header 
if args.file_type=='fastq' or args.cohort=='msbb':
     cmd = 'samtools view -H '+os.path.join(args.storage_bams,args.sample_id)+'.Aligned.sortedByCoord.out.md.bam > '  \
         + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + ".chimeric.bam"
     subprocess.check_call(cmd, shell=True, executable='/bin/bash')
     logger.info(' cmd * '+cmd)
else:
     # extract header without RG (read groups)
     cmd = 'samtools view -H '+os.path.join(args.storage_bams,args.sample_id)+'.Aligned.sortedByCoord.out.md.bam | grep -v \"RG\" > '  \
             + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + ".chimeric.bam"
     subprocess.check_call(cmd, shell=True, executable='/bin/bash')
     logger.info(' cmd * '+cmd)
     #print (cmd)
 
     # extract read groups header from raw bam
     cmd = 'samtools view -H '+args.seq_file+' | grep \"RG\" >> '  \
             + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + ".chimeric.bam"
     subprocess.check_call(cmd, shell=True, executable='/bin/bash')
     logger.info(' cmd * '+cmd)

# extract the header 
#cmd = 'samtools view -H '+os.path.join(args.storage_bams,args.sample_id)+'.Aligned.sortedByCoord.out.md.bam > '  \
#         + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + ".chimeric.bam"
#subprocess.check_call(cmd, shell=True, executable='/bin/bash')
#logger.info(' cmd * '+cmd)

 
# extract data
cmd = 'samtools view '+os.path.join(args.storage_bams,args.sample_id)+'.Aligned.sortedByCoord.out.md.bam >> '  \
         + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + ".chimeric.bam"
subprocess.check_call(cmd, shell=True, executable='/bin/bash')
logger.info(' cmd * '+cmd)

reverted_bam_dir2 = os.path.join(args.storage_circ,args.sample_id,args.sample_id+'.chimeric_reverted.bam')
cmd=os.path.join(args.src,'run_RevertSam.py ') \
    +' --tmp_dir '+args.tmp_dir \
    +' --output_by_readgroup false' \
    +' -o '+reverted_bam_dir2+' ' \
    + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + ".chimeric.bam"

logger.info(' cmd * '+cmd)
print(cmd)
subprocess.check_call(cmd,shell=True)

if (args.read_type == "PE"):
    # Read 1
    cmd = 'samtools view -h -b -f 0x0040 '+os.path.join(args.storage_circ,args.sample_id,args.sample_id)+'.chimeric_reverted.bam > ' \
    + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + "_1.bam"
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    logger.info(' cmd * '+cmd)

    # Read 2
    cmd = 'samtools view -h -b -f 0x0080 '+os.path.join(args.storage_circ,args.sample_id,args.sample_id)+'.chimeric_reverted.bam > ' \
    + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + "_2.bam"
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    logger.info(' cmd * '+cmd)

    # Align Read 1
    cmd=os.path.join(args.src,'run_STAR_circ_unified_R1.py ') \
    +args.star_index +' ' +os.path.join(args.storage_circ,args.sample_id,args.sample_id +'_R1 ') + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + "_1.bam"

    logger.info(' cmd * '+cmd)
    print(cmd)

    subprocess.check_call(cmd,shell=True)

    # Align Read 2
    cmd=os.path.join(args.src,'run_STAR_circ_unified_R2.py ') \
    +args.star_index +' ' +os.path.join(args.storage_circ,args.sample_id,args.sample_id +'_R2 ')  + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + "_2.bam"

    logger.info(' cmd * '+cmd)
    print(cmd)

    subprocess.check_call(cmd,shell=True)

    # remove the genome from memeory (updated by Abdallah)
    cmd='/opt/STAR-2.7.1a/bin/Linux_x86_64/STAR --genomeLoad Remove --genomeDir'+' '+ args.star_index
    logger.info(' cmd * '+cmd)
    print(cmd)
 
    subprocess.check_call(cmd,shell=True)

#print(args.cohort,args.tissues,args.storage_circ,args.sample_id)
#print(type(args.cohort),type(args.tissues),type(args.storage_circ),type(args.sample_id))

# align chimeric reverted 
cmd=os.path.join(args.src,'run_STAR_circ_unified.py ') \
+args.star_index +' ' +os.path.join(args.storage_circ,args.sample_id,args.sample_id +'_unified ') +os.path.join(args.storage_circ,args.sample_id,args.sample_id) +".chimeric_reverted.bam"

logger.info(' cmd * '+cmd)
print(cmd)

subprocess.check_call(cmd,shell=True)

#python3 /gscmnt/gc2645/wgs/km_test/circRNA/src/pipeline_circ.py /gscmnt/gc2645/wgs/km_test/circRNA/src/pathsCirc.json gtex --tissues Amygdala Cortex -o tmp_star_output



# 9. Post-processing: Delete unnecessary bams
#os.remove(os.path.join(args.storage_aligned,args.sample_id,args.sample_id+'.Aligned.sortedByCoord.out.bam'))
#os.remove(os.path.join(args.storage_aligned,args.sample_id,args.sample_id+'.Aligned.sortedByCoord.out.bam.bai'))
#os.remove(os.path.join(args.storage_aligned,args.sample_id,args.sample_id+'.Aligned.toTranscriptome.out.bam'))

logger.info('Sample ('+args.sample_id+') processed successfully.')
print('Sample ('+args.sample_id+') processed successfully.')


# write completed samples into a file
completed_jobs = "%s" % args.output_dir +'/completed_jobs.txt'
with open(completed_jobs, "a") as outFile:
    outFile.write(args.sample_id)
 

