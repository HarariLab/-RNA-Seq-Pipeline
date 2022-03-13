#!/usr/bin/env python3
# Author: Francois Aguet
"""
edited:
removed 'delete unneeded files' block
added 'remove fastq's for...' block
check if 'Chimeric.out.sam exists before sorting/indexing
KM: Modified to accomodate circularRNA discovery per Umber
--alignSJoverhangMin 8 /    circ needs 15**

--outSJfilterOverhangMin default: 30 12 12 12
--outSJfilterOverhangMin 15 15 15 15 ***

--seedSearchStartLmax default 50
--seedSearchStartLmax 30  *** will be slightly slower but should catch more circs

--outFilterScoreMin default: 0
--outFilterScoreMin 1		*** more stringent than default; nec

--outFilterMatchNmin default: 0
--outFilterMatchNmin 1		*** more stringent than default; nec

--chimScoreMin default: 0
--chimScoreMin 15		*** more stringent than default; nec

--genomeLoad default=NoSharedMemory
--genomeLoad circRNA: LoadAndKeep - have to use NoSharedMemory for 2-pass mode

--chimOutType: default=['WithinBAM', 'SoftClip']
--chimOutType circRNA: default=['Junctions']
"""

import argparse
import os
import subprocess
import gzip
import shutil
from datetime import datetime
import contextlib

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

parser = argparse.ArgumentParser(description='Run STAR')
parser.add_argument('index', help='Path to STAR index')
parser.add_argument('prefix', help='Prefix for output file names')
parser.add_argument('adapter_seq', help='Adapter sequences')
parser.add_argument('bam', nargs='+', help='BAM input. Format: bam_rg1 [bam_rg2], or comma-separated lists for each if multiple read group bams.')
parser.add_argument('-p', '--software_path', default=None, help='Path to STAR executable')
parser.add_argument('-o', '--outputDir', default='./', help='Output directory')
parser.add_argument('--threads', default='12', help='Number of threads')
parser.add_argument('--annotation_gtf', default=None, help='Annotation in GTF format')
parser.add_argument('--readFilesType', type=str.lower, choices=['pe','se'], help='Input type')
parser.add_argument('--outFilterMultimapNmax', default='20')
parser.add_argument('--alignSJoverhangMin', default='15', help='RNAseq = 8; more stringent than RNAseq') #KM
parser.add_argument('--alignSJDBoverhangMin', default='15')
parser.add_argument('--outFilterMismatchNmax', default='2')
parser.add_argument('--outFilterMismatchNoverLmax', default='0.3')
parser.add_argument('--alignIntronMin', default='20')
parser.add_argument('--alignIntronMax', default='0')
parser.add_argument('--alignMatesGapMax', default='0')
parser.add_argument('--outFilterType', default='Normal')
parser.add_argument('--outFilterScoreMinOverLread', default='0.66')
parser.add_argument('--outFilterMatchNminOverLread', default='0.66')
parser.add_argument('--limitSjdbInsertNsj', default='1000000')
parser.add_argument('--outSAMstrandField', default='None')
parser.add_argument('--outFilterIntronMotifs', default='None', help="Use 'RemoveNoncanonical' for Cufflinks compatibility")
parser.add_argument('--alignSoftClipAtReferenceEnds', default='Yes')
parser.add_argument('--quantMode', default='GeneCounts', nargs='+', help='Outputs read counts, and a BAM with reads in transcriptome coordinates')
parser.add_argument('--outSAMtype', default=['BAM', 'Unsorted'], nargs='+')
parser.add_argument('--outSAMunmapped', default='None', help='Keep unmapped reads in output BAM')
#parser.add_argument('--outSAMattributes', default=['NH', 'HI', 'AS', 'nM', 'NM', 'ch'], nargs='+')
parser.add_argument('--outSAMattributes', default='Standard', help='')
parser.add_argument('--chimSegmentMin', default='15', help='Minimum chimeric segment length; switches on detection of chimeric (fusion) alignments')
parser.add_argument('--chimJunctionOverhangMin', default='15', help='Minimum overhang for a chimeric junction')
parser.add_argument('--chimOutType', default=['Junctions'], help='Critical for circRNA pipeline')
parser.add_argument('--chimMainSegmentMultNmax', default='10', help='')
parser.add_argument('--genomeLoad', default='NoSharedMemory')
parser.add_argument('--sjdbFileChrStartEnd', default='-', help='SJ.out.tab file (e.g., from 1st pass). With this option, only one pass will be run')
parser.add_argument('--limitGenomeGenerateRAM', default='31000000000', help='')
parser.add_argument('--genomeChrBinNbits', default='18', help='')
parser.add_argument('--twopassMode', default='None', help='Use 1 or 2-pass mapping') # edited
parser.add_argument('--outSJfilterOverhangMin', default='15 15 15 15', help='RNAseq = 30 12 12 12') #KM
parser.add_argument('--seedSearchStartLmax', default='30', help='RNAseq = 50; slightly slower but should catch more circs') #KM
parser.add_argument('--outFilterScoreMin', default='1', help='RNAseq = 0; more stringent than RNAseq') #KM
parser.add_argument('--outFilterMatchNmin', default='1', help='RNAseq = 0; more stringent than RNAseq') #KM
parser.add_argument('--chimScoreMin', default='15', help='RNAseq = 0; more stringent than RNAseq') #KM
parser.add_argument('-z', '--debug', type=bool, default=False, help='Print command to output instead of running') # edited

args = parser.parse_args()

if not os.path.exists(args.outputDir):
    os.makedirs(args.outputDir)

if args.software_path:
    starcmd = args.software_path
else:
    starcmd = '/opt/STAR-2.7.1a/bin/Linux_x86_64/STAR' # MGI
    #starcmd = '/usr/local/genome/bin/STAR-2.7.1a' # Fenix

if args.adapter_seq !="none":
    a1 = args.adapter_seq.split(',')[0]
    a2 = args.adapter_seq.split(',')[1]

# set up command
cmd = starcmd+' --runMode alignReads --runThreadN '+args.threads+' --genomeDir '+args.index
#cmd = 'STAR --runMode alignReads --runThreadN '+args.threads+' --genomeDir '+args.index
if args.annotation_gtf is not None:  # only needed if genome index was built w/o annotation
    cmd += ' --sjdbGTFfile '+args.annotation_gtf
if args.sjdbFileChrStartEnd is None:
    cmd += ' --twopassMode '+args.twopassMode
cmd +=' --outFilterMultimapNmax '+args.outFilterMultimapNmax\
    +' --alignSJoverhangMin '+args.alignSJoverhangMin+' --alignSJDBoverhangMin '+args.alignSJDBoverhangMin\
    +' --outFilterMismatchNmax '+args.outFilterMismatchNmax+' --outFilterMismatchNoverLmax '+args.outFilterMismatchNoverLmax\
    +' --alignIntronMin '+args.alignIntronMin+' --alignIntronMax '+args.alignIntronMax+' --alignMatesGapMax '+args.alignMatesGapMax\
    +' --outFilterType '+args.outFilterType\
    +' --outFilterScoreMinOverLread '+args.outFilterScoreMinOverLread+' --outFilterMatchNminOverLread '+args.outFilterMatchNminOverLread\
    +' --outSJfilterOverhangMin '+args.outSJfilterOverhangMin\
    +' --seedSearchStartLmax '+args.seedSearchStartLmax\
    +' --outFilterScoreMin '+args.outFilterScoreMin\
    +' --outFilterMatchNmin '+args.outFilterMatchNmin\
    +' --chimScoreMin '+args.chimScoreMin\
    +' --limitSjdbInsertNsj '+args.limitSjdbInsertNsj\
    +' --readFilesIn '+','.join(args.bam)\
    +' --readFilesType '
if args.readFilesType=='pe':
   cmd += 'SAM PE'
else:
    cmd += 'SAM SE'

if args.adapter_seq !="none":
    cmd += ' --clip3pAdapterSeq '+a1+' '+a2

cmd += ' --readFilesCommand samtools view -h' \
    +' --outFileNamePrefix '+os.path.join(args.outputDir, args.prefix)+'.'\
    +' --outSAMstrandField '+args.outSAMstrandField+' --outFilterIntronMotifs '+args.outFilterIntronMotifs\
    +' --alignSoftClipAtReferenceEnds '+args.alignSoftClipAtReferenceEnds+' --quantMode '+' '+args.quantMode\
    +' --outSAMtype '+' '.join(args.outSAMtype)+' --outSAMunmapped '+args.outSAMunmapped+' --genomeLoad '+args.genomeLoad \
    +' --limitGenomeGenerateRAM '+args.limitGenomeGenerateRAM+' --genomeChrBinNbits '+args.genomeChrBinNbits
if int(args.chimSegmentMin)>0:
    cmd += ' --chimSegmentMin '+args.chimSegmentMin+' --chimJunctionOverhangMin '+args.chimJunctionOverhangMin\
        +' --chimOutType '+' '.join(args.chimOutType)+' --chimMainSegmentMultNmax '+args.chimMainSegmentMultNmax
cmd += ' --outSAMattributes '+args.outSAMattributes
if args.sjdbFileChrStartEnd is not None:
    cmd += ' --sjdbFileChrStartEnd '+args.sjdbFileChrStartEnd

if args.debug==True:
    with open(os.path.join(args.output_dir,args.prefix+'_STAR_cmd.txt'), 'w') as f:
        f.write(cmd)
        quit()

subprocess.check_call(cmd, shell=True, executable='/bin/bash')

with cd(args.outputDir):
    # set permissions
    for r,d,f in os.walk(args.prefix+'._STARpass1'):
        os.chmod(r, 0o755)
    if os.path.exists(args.prefix+'._STARgenome'):
        shutil.rmtree(args.prefix+'._STARgenome')
    if os.path.exists(args.prefix+'._STARtmp'):
        shutil.rmtree(args.prefix+'._STARtmp')


    # sort BAM (use samtools to get around the memory gluttony of STAR)

    cmd = 'samtools sort --threads '+args.threads+' -o '+args.prefix+'.Aligned.sortedByCoord.out.bam '+args.prefix+'.Aligned.out.bam'
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    os.remove(args.prefix+'.Aligned.out.bam')

    cmd = 'samtools index '+args.prefix+'.Aligned.sortedByCoord.out.bam '
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    # sort and index chimeric BAM
    if int(args.chimSegmentMin)>0 and os.path.isfile(args.prefix+'.Chimeric.out.sam'):
        cmd = 'samtools sort --threads '+args.threads+' -o '+args.prefix+'.Chimeric.out.sorted.bam '+args.prefix+'.Chimeric.out.sam '
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        cmd = 'samtools index '+args.prefix+'.Chimeric.out.sorted.bam '
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        os.remove(args.prefix+'.Chimeric.out.sam ')
