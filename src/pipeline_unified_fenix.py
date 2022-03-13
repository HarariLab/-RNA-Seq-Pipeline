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
print('output destinations:\n'+'\n\t'.join([args.storage_bams, args.storage_aligned,args.storage_qc,args.storage_quant,args.storage_tin, args.storage_circ]))

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

# 1. Log File Size

logger.info(get_size(args.seq_file))

seq_file = args.seq_file

if args.file_type=='fastq':
    local_raw = os.path.split(args.seq_file)[0]
    fastq1 = os.path.join(args.seq_file+'_1.fastq.gz')
    fastq2 = os.path.join(args.seq_file+'_2.fastq.gz')
    if os.path.exists(fastq1) and os.path.exists(fastq2):
        cmd='python3 ' +os.path.join(args.src,'run_FastqToSam_fenix.py ') \
        +fastq1 + ' ' \
        +fastq2 + ' '\
        +args.sample_id \
        +' --tmp_dir '+args.tmp_dir \
        +' --read_type '+args.read_type \
        +' -o '+local_raw

        logger.info(cmd)
        subprocess.check_call(cmd,shell=True)

        ubams = [os.path.join(local_raw,args.sample_id+'_unmapped.bam')]

    else:
        print('missing unmapped for fastq...exiting')
        quit()

else:

    # 2a. Revert BAM to unaligned BAM (uBAM)

    logger.info('Reverting BAM')

    reverted_bam_dir = os.path.join(args.storage_bams, args.sample_id+'_reverted')
    cmd='python3 ' +os.path.join(args.src,'run_RevertSam_fenix.py ') \
        +' --tmp_dir '+args.tmp_dir \
        +' -o '+reverted_bam_dir+' ' \
        + seq_file

    logger.info(' cmd * '+cmd)
    print(cmd)
    subprocess.check_call(cmd,shell=True)


    # 2b. Retrieve uBAMs

    ubams = [os.path.join(reverted_bam_dir,ubam) for ubam \
            in os.listdir(reverted_bam_dir)]

    logger.info('uBAMs generated: '+', '.join(ubams))

# 2c. Combine unaligned and aligned sequences for msbb
if args.cohort=='msbb':
    local_raw = os.path.split(args.seq_file)[0]
    unmapped_fastq = os.path.join(local_raw,args.sample_id+'.unmapped.fastq.gz')
    if os.path.exists(unmapped_fastq):
        cmd='python3 ' +os.path.join(args.src,'run_FastqToSam_fenix.py ') \
        +os.path.join(local_raw,args.sample_id+'.unmapped.fastq.gz ') \
        +args.sample_id \
        +' --tmp_dir '+args.tmp_dir \
        +' --read_type '+args.read_type \
        +' -o '+local_raw

        logger.info(cmd)
        subprocess.check_call(cmd,shell=True)

        cmd='python3 ' +os.path.join(args.src,'run_MergeSamFiles_fenix.py ') \
        +' '.join([os.path.join(local_raw, args.sample_id+'_unmapped.bam')]+ubams)+' ' \
        +args.sample_id \
        +' --tmp_dir '+args.tmp_dir \
        +' -o '+args.storage_bams

        logger.info(cmd)
        subprocess.check_call(cmd,shell=True)

        msbb_bam = os.path.join(args.storage_bams, args.sample_id+'_merged.bam')
    else:
        print('missing unmapped for fastq...exiting')
        quit()

if os.path.exists(args.tmp_dir):
    shutil.rmtree(args.tmp_dir)

# 3. Align

logger.info('Starting STAR')
star_dir = os.path.join(args.storage_aligned,args.sample_id)
if not os.path.exists(star_dir):
    os.makedirs(star_dir)

cmd ='python3 ' +os.path.join(args.src,'run_STAR_271_fenix.py ') \
    +args.star_index+' ' \
    +args.sample_id+ ' '
if args.cohort=='msbb':
    cmd+=msbb_bam
else:
    cmd+=' '.join(ubams)
cmd+=' -o '+star_dir+' ' \
    +' --readFilesType '+args.read_type

logger.info(' cmd * '+cmd)
subprocess.check_call(cmd, shell=True)

# 4. Post-alignment Picard QC

##  RNASeqMetrcs / AlignmentSumamryMetrics

logger.info('Starting post-alignment QC')
logger.info('Collecting RNA and alignment summary metrics')
cmd = 'python3 ' +os.path.join(args.src,'run_SummaryMetrics_fenix.py ') \
    +os.path.join(args.storage_aligned,args.sample_id,args.sample_id \
        +'.Aligned.sortedByCoord.out.bam ') \
    +args.sample_id +' '\
    +args.ref +' '\
    +args.ref_flat+' '\
    +args.ribosomal_int\
    +' --tmp_dir '+args.tmp_dir \
    +' -o '+args.storage_qc

logger.info(' cmd * '+cmd)
subprocess.check_call(cmd, shell=True)

## Mark Duplicates

logger.info('Marking Duplicates')
cmd ='python3 ' +os.path.join(args.src,'run_MarkDuplicates_unified_fenix.py ') \
    +os.path.join(args.storage_aligned,args.sample_id,args.sample_id \
        +'.Aligned.sortedByCoord.out.bam ') \
    +args.sample_id \
    +' --tmp_dir '+args.tmp_dir \
    +' -o '+args.storage_bams

logger.info(' cmd * '+cmd)
subprocess.check_call(cmd, shell=True)

# 5. Quantify

logger.info('Starting Salmon')
cmd ='python3 ' +os.path.join(args.src,'run_Salmon_fenix.py ') \
    +os.path.join(args.storage_aligned,args.sample_id,args.sample_id \
        +'.Aligned.toTranscriptome.out.bam ') \
    +args.salmon_index+' ' \
    +args.annot \
    +' -o '+os.path.join(args.storage_quant,args.sample_id)

logger.info(' cmd * '+cmd)
subprocess.check_call(cmd, shell=True)

logger.info('Finished quantification of: '+args.sample_id)

#6A.  Index bam for TIN
logger.info('Indexing .md.bam')
cmd = 'samtools index '+ os.path.join(args.storage_bams,args.sample_id) +'.Aligned.sortedByCoord.out.md.bam'

logger.info(' cmd * '+cmd)
subprocess.check_call(cmd, shell=True, executable='/bin/bash')
#cmd = 'samtools index '+args.prefix+'.Aligned.sortedByCoord.out.md.bam'
#subprocess.check_call(cmd, shell=True, )

logger.info('Finished Indexing: '+args.sample_id)

# 6B. Run TIN
logger.info('Running TIN')
cmd ='python ' +os.path.join(args.src,'tin.py -i') \
    +os.path.join(args.storage_bams,args.sample_id \
        +'.Aligned.sortedByCoord.out.md.bam ') \
    +' -r /40/pipelines/RNAseq/Dockerfiles/unified_rnaseq_v1/hg19_GencodeCompV19.bed' # MGI /usr/bin/hg19_GencodeCompV19.bed

logger.info(' cmd * '+cmd)
subprocess.check_call(cmd, shell=True)

if not os.path.exists(args.storage_tin):
    os.makedirs(args.storage_tin)

tin_f1=args.sample_id +'.Aligned.sortedByCoord.out.md.summary.txt'; #print('tin_f1: '+tin_f1)
tin_f2=args.sample_id +'.Aligned.sortedByCoord.out.md.tin.txt'; #print('tin_f2: '+tin_f2)
curr_dir=os.getcwd(); #print('curr_dir: '+curr_dir)
tin_f1_in=os.path.join(curr_dir,tin_f1); #print('tin_f1_in: '+tin_f1_in)
tin_f1_out=os.path.join(args.storage_tin, tin_f1); #print('tin_f1_out: '+tin_f1_out)
tin_f2_in=os.path.join(curr_dir, tin_f2); #print('tin_f2_in: '+tin_f2_in)
tin_f2_out=os.path.join(args.storage_tin, tin_f2); #print('tin_f2_out: '+tin_f2_out)

shutil.move(tin_f1_in, tin_f1_out)
shutil.move(tin_f2_in, tin_f2_out)

logger.info('Finished TIN for: '+args.sample_id)

# 7. Run circRNA
logger.info('Starting circRNA for: '+ args.sample_id)
#print(args.tissues)
if not os.path.exists(os.path.join(args.storage_circ,args.sample_id)): #,args.tissues
    os.makedirs(os.path.join(args.storage_circ,args.sample_id)) #,args.tissues

cmd = 'samtools view -H '+os.path.join(args.storage_bams,args.sample_id)+'.Aligned.sortedByCoord.out.md.bam > '  \
        + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + ".chimeric.bam"
subprocess.check_call(cmd, shell=True, executable='/bin/bash')
logger.info(' cmd * '+cmd)

cmd = 'samtools view '+os.path.join(args.storage_bams,args.sample_id)+'.Aligned.sortedByCoord.out.md.bam | grep \"ch:\" >> '  \
        + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + ".chimeric.bam"
subprocess.check_call(cmd, shell=True, executable='/bin/bash')
logger.info(' cmd * '+cmd)

reverted_bam_dir2 = os.path.join(args.storage_circ,args.sample_id,args.sample_id+'.chimeric_reverted.bam')
cmd='python3 ' +os.path.join(args.src,'run_RevertSam_circ_fenix.py ') \
    +' --tmp_dir '+args.tmp_dir \
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
    cmd='python3 ' +os.path.join(args.src,'run_STAR_circ_unified_R1_fenix.py ') \
    +args.star_index +' ' +os.path.join(args.storage_circ,args.sample_id,args.sample_id +'_R1 ') + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + "_1.bam"

    logger.info(' cmd * '+cmd)
    print(cmd)

    subprocess.check_call(cmd,shell=True)

    # Align Read 2
    cmd='python3 ' +os.path.join(args.src,'run_STAR_circ_unified_R2_fenix.py ') \
    +args.star_index +' ' +os.path.join(args.storage_circ,args.sample_id,args.sample_id +'_R2 ')  + os.path.join(args.storage_circ,args.sample_id,args.sample_id) + "_2.bam"

    logger.info(' cmd * '+cmd)
    print(cmd)

    subprocess.check_call(cmd,shell=True)

#print(args.cohort,args.tissues,args.storage_circ,args.sample_id)
#print(type(args.cohort),type(args.tissues),type(args.storage_circ),type(args.sample_id))
cmd='python3 ' +os.path.join(args.src,'run_STAR_circ_unified_fenix.py ') \
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
