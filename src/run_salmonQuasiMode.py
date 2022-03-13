#!/usr/bin/env python3
import argparse
import simplejson as json
import subprocess
import contextlib
import os
import sys
sys.path.insert(0,'/home/perezj/.local/misc')
import logging
from loggerConfig import customTime 

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

def getValue(data,path,sep='/'):
    keys = path.split(sep) 
    if(len(keys)==1): 
        key = keys.pop(0)
        return(data.get(key))
    else:
        key = keys.pop(0)
        data = data.get(key)
        return(getValue(data,sep.join(keys),sep))

parser = argparse.ArgumentParser(description='Run Salmon in Quasi-Alignment Quantification Mode')
parser.add_argument('json', type=str, help='Path the JSON file')
args = parser.parse_args()


paths = eval(open(args.json, 'r').read())
reads = paths.get('reads')

files = [s for s in os.listdir(reads.get('dir')) if os.path.isfile(os.path.join(reads.get('dir'),s))]

sampleSet = [s for s in files if os.path.splitext(s)[1]==reads.get('suffix')]

if('subset' in reads):
    with open(getValue(reads,'subset/file'),'r') as fh:
        string = fh.read().strip()
        sampleSet = [s.strip()+reads.get('suffix') for s in string.split(getValue(reads,'subset/sep'))] 

if not os.path.exists(paths.get('outDir')):
    os.makedirs(paths.get('outDir'))

logging.basicConfig(filename=paths.get('logFile'),
                    format="%(asctime)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S",
                    filemode='w',
                    level=logging.INFO)

logging.Formatter.converter = customTime

logging.info('Started Analysis: '+paths.get('analysisName'))
logging.info('description: '+paths.get('description'))

for idx, readFn in enumerate(sampleSet):

    readsPath = os.path.join(reads.get('dir'),readFn)
    sampleName = os.path.splitext(readFn)[0]
    outDir = os.path.join(paths.get('outDir'),sampleName)

    if not os.path.exists(outDir):
        os.makedirs(outDir)

    with cd(paths.get('outDir')):
        logging.info('Processing '+sampleName)
    
        logging.info('Starting SamToFastq')
    
        cmd=getValue(paths,'pipeline/SamToFastq')\
        +' --memory 5'\
        +' -o '+sampleName\
        +' --prefix '+sampleName+' '\
        +readsPath 
        
        logging.info(' * command: '+cmd)
        subprocess.check_call(cmd,shell=True)

        logging.info('Starting Salmon')
    
        cmd=getValue(paths,'pipeline/Salmon')\
        +' --threads 12'\
        +' --libType A'\
        +' -i '+getValue(paths,'genome/index')\
        +' -1 '+os.path.join(sampleName,sampleName+'_1.fastq.gz')\
        +' -2 '+os.path.join(sampleName,sampleName+'_2.fastq.gz')\
        +' -o '+sampleName

        logging.info(' * command: '+cmd)
        subprocess.check_call(cmd,shell=True)

        if(os.path.exists(os.path.join(sampleName,sampleName+'_1.fastq.gz'))):
            os.remove(os.path.join(sampleName,sampleName+'_1.fastq.gz'))
            os.remove(os.path.join(sampleName,sampleName+'_2.fastq.gz'))

logging.info('Finished Salmon Quantification')
