#!/usr/bin/env python3
import argparse
import simplejson as json
import subprocess
import contextlib
import shutil
import os
import sys
from datetime import datetime
import shlex

def file_exists_remote(host,path):
    """Test if a file exists at path on a host accessible with SSH."""
    status = subprocess.call(
            ['ssh', host, 'test -f {}'.format(shlex.quote(path))])
    if status == 0:
        return True
    if status == 1:
        return False
def get_remote_samples(bam_dir, host=None, sep=None, ext='.bam', idx=None, use_full_path=False):
    if host:
        cmd='ssh '+host+' ls '+bam_dir
        cnx = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
        dirlist = [f.strip() for f in cnx.communicate()[0].decode('utf-8').split('\n') if f]
        bam_list = [os.path.splitext(f.strip())[0] for f in dirlist if os.path.splitext(f)[1]==ext] 
    else:
        bam_list = [os.path.splitext(f.strip())[0] for f in os.listdir(bam_dir) if os.path.splitext(f)[1]==ext]

    if use_full_path:
        return(bam_list)
    if not idx:
        idx=0 
    return([f.split(sep)[idx] for f in bam_list])

def get_quantified_samples(tissue_dir):
    cnx = subprocess.Popen(shlex.split('find '+tissue_dir+' -name "quant.sf"'), stdout=subprocess.PIPE)
    dirlist = [f.strip() for f in cnx.communicate()[0].decode('utf-8').split('\n') if f]
    return([f.split('/')[-2] for f in dirlist]) 
 
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

parser = argparse.ArgumentParser(description='Run pipeline from RNA-Seq JSON file')
parser.add_argument('json', type=str, help='Path to the JSON file')
parser.add_argument('study', type=str, help='Studies from which to align and quantify reads from.')
parser.add_argument('tissues', nargs='+', type=str, help='Tissue from a tissue to process. Format: tissue1 [tissue2], or comma-separated list')
parser.add_argument('-d','--download', action='store_true', default=False,help='Download bams only (i.e., no processing)')
parser.add_argument('-s','--subset',type=str, default=None,help='File with sample subset to process')
parser.add_argument('--check_status', action='store_true',help='Only show the number of samples processed from total')


args = parser.parse_args()

with open(args.json) as data:
    paths = json.load(data)

cohorts = paths.get('cohorts')

if args.study not in [s for s in cohorts.keys()]:
   raise ValueError(study+' not in json file. Found: '+', '.join([s for s in cohorts.keys()]))


print('server host: '+paths.get('host'))

for cohort, entities in cohorts.items():
    if cohort != args.study:
        continue 

    print('study: '+cohort)

    tissues = entities.get('tissues')
    for i in args.tissues: 
        if i not in [t for t in tissues.keys()]:
            raise ValueError(i+' not in json file. Found: '+', '.join([t for t in tissues.keys()]))

    for tissue, data in tissues.items():
        if tissue not in args.tissues:
            continue
        print('tissue: '+tissue)
        missing = [i for i in ['dir','tokens','id_sep','read_type','read_len'] if i not in data.keys()]
        if missing:
            raise ValueError('incomplete data in json file. Missing: '+', '.join(missing)) 
            

        output_dir = os.path.join(paths.get('outdir'),cohort)
        remote_raw_dir = data.get('dir')
        local_raw_dir = os.path.join(output_dir, getValue(paths,'storage/raw'),tissue)
        ext = data.get('tokens')[0].get('ext')
        file_type  = data.get('tokens')[0].get('type') 

        
        # Check for incomplete samples
#        sample_set = get_remote_samples(remote_raw_dir, paths.get('host'), data.get('id_sep'), ext)
#
#        if os.path.exists(os.path.join(output_dir, getValue(paths,'storage/quant'))):
#            quantified_set = get_quantified_samples(os.path.join(output_dir, getValue(paths,'storage/quant')))
#        else:
#            quantified_set = []
#
#        print('Sample size: '+str(len(sample_set)))
#
#        incomplete_set = [i for i in sample_set if i not in quantified_set]
#
#        print(str(len(incomplete_set))+'/'+str(len(sample_set))+' samples not quantified.')
#        if args.check_status:
#            continue
#        
        # 1. Download input files
        if args.download:
            print('Downloading samples...')
    
            if not os.path.exists(local_raw_dir):
                os.makedirs(local_raw_dir)
    
            cmd='rsync -avh "'+paths.get('host')+':'+remote_raw_dir+'/ " '+local_raw_dir
            print(cmd) 
            subprocess.check_call(cmd, shell=True)
            continue 
     
    
        # 2. Download Genome/Annotation Files
#        print('Downloading salmon/star index and annotation files')
#
        local_genome_dir = os.path.join(paths.get('outdir'),'genome')
#
#        if not os.path.exists(local_genome_dir):
#            os.makedirs(local_genome_dir)
#
#        remote_star_index = [t.get('dir') for t in getValue(paths,'genome/star_index') \
#                     if t.get('len')==data.get('read_len')][0]
#        remote_salmon_index = getValue(paths, 'genome/salmon_index')
#        annot = getValue(paths, 'genome/annot')
#        ref = getValue(paths, 'genome/ref')
#        ref_flat = getValue(paths, 'genome/ref_flat')
#        ribosomal_int = getValue(paths, 'genome/ribosomal_int')
#
#        cmd='rsync -avh "'\
#            +paths.get('host')\
#            +':'+' '.join([remote_star_index, remote_salmon_index, annot, ref_flat, ribosomal_int]) \
#            +'" '+local_genome_dir
#        
#        print(cmd) 
#        subprocess.check_call(cmd, shell=True)
#    
#        # 3. Get input files

        raw = os.listdir(local_raw_dir) 
        tokens = data.get('tokens')
        read_groups = [s.strip() for s in raw if  os.path.splitext(s)[1]==ext]
        if file_type == 'fastq':
            read_groups = [r.split(data.get('id_sep'))[0] for r in read_groups ] 
            read_groups = set(read_groups)
            print(read_groups)

        ## Process a select subset of samples ##
        if(args.subset):
            with open(args.subset,'r') as fh:
                string = fh.read().strip()
                subset = [s.strip() for s in string.split('\n')]
                unprocessed_list = [s for s in read_groups if s in subset]
                if not unprocessed_list: # empty
                    print('No samples in subset present for tissue "'+tissue+'"')
        else:
            unprocessed_list = read_groups 

#            unprocessed_list = [i for i in read_groups if i.split(data.get('id_sep'))[0] in incomplete_set] 
#            print(str(len(unprocessed_list))+'/'+str(len(read_groups))+' samples to process.')
# 
        for i, rg in enumerate(unprocessed_list):

            job_name = cohort+"_"+tissue[:2]+"_"+str(i+1)
            file_type = file_type
            host = paths.get('host')
            src = paths.get('src')
            sample_id = rg.split(data.get('id_sep'))[0]
            path_to_local_seq_file= os.path.join(local_raw_dir, rg)
            ref = os.path.join(local_genome_dir, os.path.split(getValue(paths,'genome/ref'))[1])
            star_index = os.path.join(local_genome_dir,data.get('read_len'))
            salmon_index = os.path.join(local_genome_dir, os.path.split(getValue(paths,'genome/salmon_index'))[1])
            ref_flat = os.path.join(local_genome_dir, os.path.split(getValue(paths,'genome/ref_flat'))[1])
            ribosomal_int = os.path.join(local_genome_dir, os.path.split(getValue(paths,'genome/ribosomal_int'))[1])
            annot = os.path.join(local_genome_dir, os.path.split(getValue(paths,'genome/annot'))[1])
            read_len = data.get('read_len')
            read_type = data.get('read_type')
            tmp_dir = os.path.join(paths.get('tmp'),sample_id) 
            storage_aligned = os.path.join(output_dir, getValue(paths,'storage/aligned'),sample_id) 
            storage_bams = os.path.join(output_dir, getValue(paths,'storage/bams')) 
            storage_qc = os.path.join(output_dir,  getValue(paths,'storage/qc'),sample_id) 
            storage_quant = os.path.join(output_dir, getValue(paths,'storage/quant'), sample_id) 

            cmd="bsub -Is -q research-hpc -J "+job_name \
            +" -M 50000000 " \
            +" -R 'select [mem>50000 && gtmp > 15] rusage[mem=50000, gtmp=15]'" \
            +" -a 'docker(apollodorus/brain-eqtl:rnaseq2)'" \
            +" /bin/bash -c '"+os.path.join(src,"pipeline.py")+" " \
            +" ".join([host,src,path_to_local_seq_file, file_type, cohort, sample_id, ref, star_index, salmon_index, \
              ref_flat, ribosomal_int, annot, read_len, read_type, tmp_dir, output_dir, storage_aligned, storage_bams, storage_qc, storage_quant])+"'"

#            print(cmd)
            subprocess.Popen(cmd,shell=True)
