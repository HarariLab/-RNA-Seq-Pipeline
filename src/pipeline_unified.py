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

#logger.info(get_size(args.seq_file))

seq_file = args.seq_file
local_raw = os.path.split(args.seq_file)[0]

if args.file_type=='fastq':
    #local_raw = os.path.split(args.seq_file)[0]
    fastq1 = os.path.join(args.seq_file+'_R1.fastq.gz')
    fastq2 = os.path.join(args.seq_file+'_R2.fastq.gz')
    if os.path.exists(fastq1) and os.path.exists(fastq2):
         cmd=os.path.join(args.src,'run_FastqToSam.py ') \
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
    cmd=os.path.join(args.src,'run_RevertSam.py') \
        +' --tmp_dir '+args.tmp_dir \
        +' -o '+reverted_bam_dir+' ' \
        +' -c '+args.cohort+' ' \
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
    #local_raw = os.path.split(args.seq_file)[0]
    #unmapped_fastq = os.path.join(local_raw,args.sample_id+'.unmapped.fastq.gz')
    unmapped_fastq = glob.glob(local_raw+'/'+args.sample_id+'*.unmapped.*.gz')[0]

    if os.path.exists(unmapped_fastq):
        cmd=os.path.join(args.src,'run_FastqToSam.py ') \
        +unmapped_fastq +' ' \
        +args.sample_id \
        +' --tmp_dir '+args.tmp_dir \
        +' --read_type '+args.read_type \
        +' -o '+local_raw

        logger.info(cmd)
        subprocess.check_call(cmd,shell=True)

        cmd=os.path.join(args.src,'run_MergeSamFiles.py ') \
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

#2d. Combine unaligned and aligned sequences for CommonMind
##if args.cohort=='CommonMind':
##    #local_raw = os.path.split(args.seq_file)[0]
##    cmd=os.path.join(args.src,'run_MergeSamFiles.py ') \
##    +' '.join([os.path.join(local_raw, args.sample_id+'.unmapped.bam')]+ubams)+' ' \
##    +args.sample_id \
##    +' --tmp_dir '+args.tmp_dir \
##    +' -o '+args.storage_bams
 
##    logger.info(cmd)
##    subprocess.check_call(cmd,shell=True)
 
##    msbb_bam = os.path.join(args.storage_bams, args.sample_id+'_merged.bam')
 

# 3. Align

logger.info('Starting STAR')
star_dir = os.path.join(args.storage_aligned,args.sample_id)
if not os.path.exists(star_dir):
    os.makedirs(star_dir)

cmd =os.path.join(args.src,'run_STAR_271.py ') \
    +args.star_index+' ' \
    +args.sample_id+ ' ' \
    +args.adapter_seq+ ' '
if args.cohort in [ 'msbb', 'CommonMind' ]:
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
cmd = os.path.join(args.src,'run_SummaryMetrics.py ') \
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
cmd = os.path.join(args.src,'run_MarkDuplicates_unified.py ') \
    +os.path.join(args.storage_aligned,args.sample_id,args.sample_id \
        +'.Aligned.sortedByCoord.out.bam ') \
    +args.sample_id \
    +' --tmp_dir '+args.tmp_dir \
    +' -o '+args.storage_bams

logger.info(' cmd * '+cmd)
subprocess.check_call(cmd, shell=True)

# 5. Quantify

logger.info('Starting Salmon')
cmd = os.path.join(args.src,'run_Salmon.py ') \
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
cmd = 'tin.py -i ' \
    +os.path.join(args.storage_bams,args.sample_id \
        +'.Aligned.sortedByCoord.out.md.bam ') \
    +' -r '+ args.annot_bed 

logger.info(' cmd * '+cmd)
subprocess.check_call(cmd, shell=True)

if not os.path.exists(args.storage_tin):
    os.makedirs(args.storage_tin)

tin_f1=args.sample_id +'.Aligned.sortedByCoord.out.md.summary.txt'; #print('tin_f1: '+tin_f1)
tin_f2=args.sample_id +'.Aligned.sortedByCoord.out.md.tin.xls'; #print('tin_f2: '+tin_f2)
curr_dir=os.getcwd(); #print('curr_dir: '+curr_dir)
tin_f1_in=os.path.join(curr_dir,tin_f1); #print('tin_f1_in: '+tin_f1_in)
tin_f1_out=os.path.join(args.storage_tin, tin_f1); #print('tin_f1_out: '+tin_f1_out)
tin_f2_in=os.path.join(curr_dir, tin_f2); #print('tin_f2_in: '+tin_f2_in)
tin_f2_out=os.path.join(args.storage_tin, tin_f2); #print('tin_f2_out: '+tin_f2_out)

shutil.move(tin_f1_in, tin_f1_out)
shutil.move(tin_f2_in, tin_f2_out)

logger.info('Finished TIN for: '+args.sample_id)

######################################################################################################################
# 7. Run circRNA
######################################################################################################################
if (args.run_circ=="run_circ"):

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
       +args.star_index +' '+os.path.join(args.storage_circ,args.sample_id,args.sample_id +'_R1 ')+args.adapter_seq+' '+os.path.join(args.storage_circ,args.sample_id,args.sample_id) + "_1.bam"
       logger.info(' cmd * '+cmd)
       print(cmd)

       subprocess.check_call(cmd,shell=True)

       # Align Read 2
       cmd=os.path.join(args.src,'run_STAR_circ_unified_R2.py ') \
       +args.star_index +' '+os.path.join(args.storage_circ,args.sample_id,args.sample_id +'_R2 ')+args.adapter_seq+' '+os.path.join(args.storage_circ,args.sample_id,args.sample_id) + "_2.bam"
       logger.info(' cmd * '+cmd)
       print(cmd)

       subprocess.check_call(cmd,shell=True)

       # remove the genome from memeory
       cmd='/opt/STAR-2.7.1a/bin/Linux_x86_64/STAR --genomeLoad Remove --genomeDir'+' '+ args.star_index
       logger.info(' cmd * '+cmd)

       subprocess.check_call(cmd,shell=True)

   #print(args.cohort,args.tissues,args.storage_circ,args.sample_id)
   #print(type(args.cohort),type(args.tissues),type(args.storage_circ),type(args.sample_id))

   # Align chimeric reverted bam
   cmd=os.path.join(args.src,'run_STAR_circ_unified.py ') \
   +args.star_index +' ' +os.path.join(args.storage_circ,args.sample_id,args.sample_id +'_unified ')+args.adapter_seq+' '+os.path.join(args.storage_circ,args.sample_id,args.sample_id) +".chimeric_reverted.bam" \
   +" --readFilesType "+args.read_type

   logger.info(' cmd * '+cmd)
   print(cmd)

   subprocess.check_call(cmd,shell=True)
######################################################################################################################

#python3 /gscmnt/gc2645/wgs/km_test/circRNA/src/pipeline_circ.py /gscmnt/gc2645/wgs/km_test/circRNA/src/pathsCirc.json gtex --tissues Amygdala Cortex -o tmp_star_output

#################################################
## 8.delete intermmediate bam files 
#################################################
#lete raw data 
RawData_path = os.path.split(args.seq_file)[0]
raw_files = glob.glob(os.path.join(RawData_path, args.sample_id + '*'))
for rf in raw_files:
    os.remove(rf)
 
if args.file_type=='bam':
   if os.path.exists(reverted_bam_dir) and os.path.isdir(reverted_bam_dir):
      shutil.rmtree(reverted_bam_dir)
 

if args.cohort in ['msbb', 'CommonMind']:
   os.remove(msbb_bams)

## remove Transcriptome.out.bam file 
os.remove(os.path.join(args.storage_aligned,args.sample_id,args.sample_id+'.Aligned.toTranscriptome.out.bam'))
 
## remove circRNA bam files 
if (args.run_circ=="run_circ"):
   os.remove(os.path.join(args.storage_circ,args.sample_id,args.sample_id) +".chimeric_reverted.bam")
   os.remove(os.path.join(args.storage_circ,args.sample_id,args.sample_id) + ".chimeric.bam")
   os.remove(os.path.join(args.storage_circ,args.sample_id,args.sample_id) + "_1.bam")
   os.remove(os.path.join(args.storage_circ,args.sample_id,args.sample_id) + "_2.bam")
 
# 8.1 move marked duplicates files to 01.-QC
shutil.move(os.path.join(args.storage_bams,args.sample_id+'.marked_dup_metrics.txt'), os.path.join(args.storage_qc,args.sample_id+'.marked_dup_metrics.txt') )

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
     if (args.run_circ=="run_circ"): 
     	cmd='samtools view -C '+os.path.join(args.storage_circ,args.sample_id,args.sample_id)+'_unified.Aligned.sortedByCoord.out.bam -T '+args.ref \
     	+' -o '+os.path.join(args.storage_circ,args.sample_id,args.sample_id) +'_unified.Aligned.sortedByCoord.out.cram'
     	logger.info(' cmd * '+cmd)
     	subprocess.check_call(cmd,shell=True)

     # 9.3 move crams to the box 
     cmd='rclone copy '+os.path.join(args.storage_bams,args.sample_id)+'.Aligned.sortedByCoord.out.md.cram '+args.remote_storage_name+':'+args.star_storage_box
     logger.info(' cmd * '+cmd)
     print(cmd)
     subprocess.check_call(cmd,shell=True)

     # move circRNA cram
     if (args.run_circ=="run_circ"):
     	cmd='rclone copy '+os.path.join(args.storage_circ,args.sample_id,args.sample_id) +'_unified.Aligned.sortedByCoord.out.cram '+args.remote_storage_name+':'+args.circ_storage_box
     	logger.info(' cmd * '+cmd)
     	print(cmd)
     	subprocess.check_call(cmd,shell=True)

     # delete aligned data
     os.remove(os.path.join(args.storage_aligned,args.sample_id,args.sample_id+'.Aligned.sortedByCoord.out.bam'))
     os.remove(os.path.join(args.storage_aligned,args.sample_id,args.sample_id+'.Aligned.sortedByCoord.out.bam.bai'))
     #os.remove(os.path.join(args.storage_aligned,args.sample_id,args.sample_id+'.Aligned.toTranscriptome.out.bam'))
     os.remove(os.path.join(args.storage_bams,args.sample_id+'.Aligned.sortedByCoord.out.md.bam'))
     os.remove(os.path.join(args.storage_bams,args.sample_id+'.Aligned.sortedByCoord.out.md.bam.bai'))
     os.remove(os.path.join(args.storage_bams,args.sample_id+'.Aligned.sortedByCoord.out.md.cram'))

     if (args.run_circ=="run_circ"):
        os.remove(os.path.join(args.storage_circ,args.sample_id,args.sample_id +'_unified.Aligned.sortedByCoord.out.cram'))


# message to show the completion of the pipeline
logger.info('Sample ('+args.sample_id+') processed successfully.')
print('Sample ('+args.sample_id+') processed successfully.')

# write completed samples into a file
completed_jobs = "%s" % args.output_dir +'/completed_jobs.txt'
with open(completed_jobs, "a") as outFile:
    outFile.write(args.sample_id +'\n')


# 8. Transfer data to fenix/dragon

#echo "Transfer folder 01.-QC ..."
#eteleeb@fenix.psych.wucon.wustl.edu:/40/Cruchaga_Data/bulkRNASeq/$study/pipeline-v1-hg38/02.-ProcessedData/01.-QC

#cmd = 'rsync -ravh '+ os.path.join(args.storage_qc, args.sample_id + '* ' + args.host +  


