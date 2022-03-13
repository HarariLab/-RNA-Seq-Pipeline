#!/usr/bin/env python3

import argparse
import subprocess
from datetime import datetime
import os, shutil
from collections import defaultdict
import contextlib

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

parser = argparse.ArgumentParser(description='Collect alignment summary metrics using Picard.')
parser.add_argument('input_bam', type=str, help='BAM file')
parser.add_argument('prefix', type=str, help='Prefix for output files; usually <sample_id>')
parser.add_argument('fasta', type=str, help='Reference sequence fasta')
parser.add_argument('ref_flat', type=str, help='Gene annotations in refFlat form')
parser.add_argument('ribosomal_int', type=str, help='Location of rRNA sequences in genome, in interval_list format')
parser.add_argument('--strand', type=str, default='NONE', help='For strand-specific library prep')
parser.add_argument('--adapter', type=str, default='null',help='List of adapter sequences to use when processing the alignment metrics')
parser.add_argument('-o', '--output_dir', default=os.getcwd(), help='Output directory')
parser.add_argument('-m', '--memory', default='8', type=str, help='Memory, in GB')
parser.add_argument('--max_records_in_ram',
                    type=str,
                    default='150000',
                    help='Specifies the number of records stored in RAM before spilling to disk.')
parser.add_argument('--tmp_dir',
                    type=str,
                    default='~/tmp',
                    help='') # None MGI
parser.add_argument('--validation_stringency', type=str.upper,
                    choices=['STRICT','LENIENT','SILENT'],
                    default='SILENT',
                    help='Validation stringency for all SAM files read by program')
parser.add_argument('--jar', default='/opt/picard-tools/picard.jar', help='Path to Picard jar') # MGI
#parser.add_argument('--jar', default='/usr/local/genome/picard_latest/picard.jar', help='Path to Picard jar') # comment for MGI
args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
if not os.path.exists(args.tmp_dir):
    os.makedirs(args.tmp_dir)

with cd(args.output_dir):
    cmd='java -jar -Xmx'+args.memory+'g '+args.jar \
        +' CollectRnaSeqMetrics I='+args.input_bam \
        +' O='+args.prefix+'.RNA_Metrics.txt' \
        +' REF_FLAT='+args.ref_flat \
        +' RIBOSOMAL_INTERVALS='+args.ribosomal_int \
        +' STRAND_SPECIFICITY='+args.strand \
        +' TMP_DIR='+args.tmp_dir \
        +' MAX_RECORDS_IN_RAM='+args.max_records_in_ram \
        +' VALIDATION_STRINGENCY='+args.validation_stringency
    subprocess.check_call(cmd, shell=True)

    if os.path.exists(args.tmp_dir):
        shutil.rmtree(args.tmp_dir)
    if not os.path.exists(args.tmp_dir):
        os.makedirs(args.tmp_dir)

    cmd='java -jar -Xmx'+args.memory+'g '+args.jar \
        +' CollectAlignmentSummaryMetrics I='+args.input_bam \
        +' ADAPTER_SEQUENCE='+args.adapter \
        +' O='+args.prefix+'.Summary_Metrics.txt' \
        +' TMP_DIR='+args.tmp_dir \
        +' R='+args.fasta \
        +' MAX_RECORDS_IN_RAM='+args.max_records_in_ram \
        +' VALIDATION_STRINGENCY='+args.validation_stringency
    print(cmd)
    subprocess.check_call(cmd, shell=True)

    if os.path.exists(args.tmp_dir):
        shutil.rmtree(args.tmp_dir)

print('Finished QC Summary Metrics')
