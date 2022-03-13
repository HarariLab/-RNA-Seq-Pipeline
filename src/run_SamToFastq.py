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
import logger

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

logger = logging.getLogger('run_SamToFastq')

parser = argparse.ArgumentParser(description='Convert BAM to FASTQ using'\
                                 + 'SamToFastq from Picard.')

#parser.add_argument('infile_array',type=str, 
#                   help='Path to alignment .bam files to be processed')
parser.add_argument('bam_file', type=str, help='BAM file')
parser.add_argument('-p','--prefix',type=str,default='Reads',
                    help='Prefix for output files')
parser.add_argument('-o', '--output_dir',default=os.getcwd(),
                    help='Output directory')
parser.add_argument('-m', '--memory', default='4',type=str, 
                    help='Memory in GB')
parser.add_argument('--jar', default='/opt/picard-tools/picard.jar', 
                    help='Path to Picard.jar')

parser.add_argument('--gzip', type=str.lower, default='1',
                    help='gzip compression level for FASTQs; see "man gzip"')


parser.add_argument('--include_non_pf_reads', type=str.lower, 
                    choices=['true','false'],
                    default='true',
                    help='PF:passed filtering')

parser.add_argument('--include_non_primary_alignments',
                    type=str.lower,choices=['true','false'], 
                    default='false',
                    help='Sets non-primary alignments option')

args = parser.parse_args()

logger.info('Starting SamToFastq')

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
    
with cd(args.output_dir):
    fastq1 = args.prefix+'_1.fastq.gz'
    fastq2 = args.prefix+'_2.fastq.gz'
    fastq0 = args.prefix+'_unpaired.fastq.gz'
    
    subprocess.check_call('mkfifo read1_pipe read2_pipe read0_pipe',shell=True)
    
    subprocess.check_call('gzip -'+args.gzip+' -c < read1_pipe > '+fastq1+' &',shell=True)
    subprocess.check_call('gzip -'+args.gzip+' -c < read2_pipe > '+fastq2+' &',shell=True)
    subprocess.check_call('gzip -'+args.gzip+' -c < read0_pipe > '+fastq0+' &',shell=True)
    
    subprocess.check_call('java -jar -Xmx' + args.memory +'g ' + args.jar\
        +' SamToFastq INPUT='+args.bam_file\
        +' INCLUDE_NON_PF_READS='+ args.include_non_pf_reads\
        +' INCLUDE_NON_PRIMARY_ALIGNMENTS='+args.include_non_primary_alignments\
        +' VALIDATION_STRINGENCY=SILENT'\
        +' FASTQ=read1_pipe SECOND_END_FASTQ=read2_pipe UNPAIRED_FASTQ=read0_pipe',
            shell=True)
    
    subprocess.check_call('rm read1_pipe read2_pipe read0_pipe',shell=True)
    
    with open(fastq0, 'rb') as f0:
        f0.seek(-4,2)
        if struct.unpack('<I',f0.read(4))[0]==0: # empty file
            os.remove(fastq0)
            
logger.info('Finished SamToFastq')
