#!/usr/bin/env python3
import argparse
import os
import struct
import subprocess
import shutil
from datetime import datetime
import contextlib
import logging

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description='Revert Bam to unaligned'\
                                 + 'read groups from Picard.')

parser.add_argument('bam_file', type=str, help='BAM file')
parser.add_argument('-o', '--output_dir', type=str, help='Output bam name')
parser.add_argument('-m', '--memory', default='8',type=str,
                    help='Memory in GB')
parser.add_argument('--jar', default='/opt/picard-tools/picard.jar',
                    help='Path to Picard.jar')
#parser.add_argument('--output_by_readgroup', type=str.lower, #changed from upper
#                    choices=['true','false'],
#                    default='false',
#                    help='Separates unaligned bams by read group')
parser.add_argument('--tmp_dir', type=str.lower,
                    default='/data/tmp', 
                    help='Directory for tmp files during reversion')
parser.add_argument('--sanitize', type=str.upper,
                    choices=['true','false'],
                    default='true',
                    help='Discard reads to produce consistent output BAM')
parser.add_argument('--validation_stringency', type=str.upper,
                    choices=['STRICT','LENIENT','SILENT'],
                    default='LENIENT',
                    help='Validation stringency for all SAM files read by program')
parser.add_argument('--max_records_in_ram',
                    type=str,
                    default='150000',
                    help='Specifies the number of records stored in RAM before spilling to disk.')
parser.add_argument('--max_discard_fraction',
                    type=str,
                    default='0.005',
                    help='Sets max fraction of reads to be discarded by SANITIZE=TRUE')
parser.add_argument('--attribute_to_clear',
                    nargs='+',
                    type=str,
                    default=['NH','HI','nM','NM','ch'],
                    help='Format: NH [HI nM...] or comma-separated list for each BAM attribute')
args = parser.parse_args()

#if not os.path.exists(args.output_dir):
#    os.makedirs(args.output_dir)
if not os.path.exists(args.tmp_dir):
    os.makedirs(args.tmp_dir)


cmd='java -jar -Xmx' + args.memory +'g ' + args.jar \
        +' RevertSam I='+args.bam_file \
        +' OUTPUT_BY_READGROUP=false' \
        +' SANITIZE='+args.sanitize \
        +' MAX_DISCARD_FRACTION='+args.max_discard_fraction \
        +' MAX_RECORDS_IN_RAM='+args.max_records_in_ram \
        +' VALIDATION_STRINGENCY='+args.validation_stringency
if(args.attribute_to_clear):
        cmd+=' ATTRIBUTE_TO_CLEAR=' \
        +' ATTRIBUTE_TO_CLEAR='.join(args.attribute_to_clear) \
        +' TMP_DIR='+args.tmp_dir
cmd+=' O='+args.output_dir

logger.info(cmd)
subprocess.check_call(cmd, shell=True)

if os.path.exists(args.tmp_dir):
    shutil.rmtree(args.tmp_dir)
