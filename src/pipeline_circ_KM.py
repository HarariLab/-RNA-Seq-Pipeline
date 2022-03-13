#!/usr/bin/env python3
import argparse
import json
import subprocess
import contextlib
import glob
from datetime import datetime
import concurrent.futures
import os
import time

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

def bsub(cmd, mem=50, gtmp=12, docker_image="de1malo1bonum/rnaseq_km:unified.0.1.1", job_name='circRNA', debug=False): #apollodorus/star:2.7.1

    mem = str(mem)
    gtmp=str(gtmp)
    bsub = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false" \
    + " bsub -Is -q research-hpc" \
    + " -J " + job_name \
    + " -M "+mem+"000000" \
    + " -R 'select[mem>"+mem+"000 && gtmp > "+gtmp+"] rusage[mem="+str(mem)+"000, gtmp="+str(gtmp)+"]'" \
    + " -a 'docker("+docker_image+")'" \
    + " /bin/bash -c '"+cmd+"'"

    if debug:
        print('* '+bsub)
        quit()

    return(bsub)

# intended outdir: /gscmnt/gc2645/wgs/RNAseq/circRNA/data
parser = argparse.ArgumentParser(description='STAR alignment for circular RNA-seq.')
parser.add_argument('json', type=str, help='Path to JSON file with input paths.')
parser.add_argument('study', type=str, help='')
parser.add_argument('-m', '--max_workers', type=int, default=500, help='')
parser.add_argument('-s', '--subset', type=str, help='')
parser.add_argument('-t', '--tissues', type=str,  help='') #nargs='+',
parser.add_argument('-o', '--output_dir', type=str, default='.', help='Path to output directory for all processes')
parser.add_argument('--debug', action='store_true', help='')

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

with open(args.json) as fp:
    paths = json.load(fp)

if args.subset:
    with open(args.subset) as fp:
        subset = fp.read().strip().split('\n')

def run_star(arg_obj):

    exc_dir = os.path.dirname(os.path.abspath(__file__))
    r1_pass = ' '.join([os.path.join(exc_dir, 'run_STAR_circ_unified_R1.py'),
        arg_obj.genome,
        arg_obj.sample_id+'_R1',
        arg_obj.r1,
        ' -o ' + os.path.join(arg_obj.outdir)])

    r2_pass = ' '.join([os.path.join(exc_dir, 'run_STAR_circ_unified_R2.py'),
        arg_obj.genome,
        arg_obj.sample_id+'_R2',
        arg_obj.r2,
        ' -o ' + os.path.join(arg_obj.outdir)])

    pair_pass = ' '.join([os.path.join(exc_dir, 'run_STAR_circ_unified.py'),
        arg_obj.genome,
        arg_obj.sample_id+'_unified',
        arg_obj.pair,
        ' -o ' + os.path.join(arg_obj.outdir)])

    if args.debug:

        [print(' * ' + cmd) for cmd in (r1_pass, r2_pass, pair_pass)]
        quit()

    for cmd in (r1_pass, r2_pass, pair_pass):
        subprocess.check_call(bsub(cmd), shell=True)

    return(0)

class arg():
    def __init__(self, sample_id, r1_bam, r2_bam, pair_bam, genome_dir, outdir):

        self.sample_id = sample_id
        self.r1 = r1_bam
        self.r2 = r2_bam
        self.pair = pair_bam
        self.genome = genome_dir
        self.outdir = outdir

for study_name, tissue_dict in paths.get("cohorts").items():
    #print(args.cohort)
    if study_name.lower() != args.study.lower():
        continue
    print('tissue_dict: ', end=''); print(tissue_dict)
    print('tissue_dict tissues: ', end=''); print(tissue_dict.get("tissues"))
    print('tissue_dict items: ', end=''); print(tissue_dict.get("tissues").items())
    for tissue_name, input_dict in tissue_dict.get('tissues'): #.items() # YOU ARE HERE

        if args.tissues and (tissue_name.lower() not in [_.lower() for _ in args.tissues]):
            continue

        r1_bams = glob.glob(os.path.join(input_dict.get('input'), '*_1.bam'))
        r2_bams = glob.glob(os.path.join(input_dict.get('input'), '*_2.bam'))
        pair_bams = glob.glob(os.path.join(input_dict.get('input'), '*reverted.bam'))

        # sort all bam lists inplace by sample id
        def sid(x):
            return(os.path.split(x)[1].split('.chimeric')[0])

        [l.sort(key=lambda x: sid(x)) for l in (r1_bams, r2_bams, pair_bams)]

        num_jobs = 0
        futures = list()

        executor = concurrent.futures.ProcessPoolExecutor(max_workers=args.max_workers)
        print('r1_bams: ', end=''); print(r1_bams)
        print('r2_bams: ', end=''); print(r2_bams)
        print('pair_bams: ', end=''); print(pair_bams)
        for r1, r2, pair in zip(r1_bams, r2_bams, pair_bams):
            if args.subset and (sid(r1) not in subset):
                continue

            assert sid(r1) == sid(r2) == sid(pair)

            outdir = os.path.join(args.output_dir, study_name, tissue_name, sid(r1))
            star_args = arg(sid(r1), r1, r2, pair, input_dict.get('star_genome'), outdir)

            if args.debug:
                run_star(star_args)
                break

            futures.append(executor.submit(run_star, star_args))
            num_jobs = num_jobs + 1

    print('[ {} ] Submitted {} jobs for tissue "{}".'.format(datetime.now().strftime('%b %d %H:%M:%S'), num_jobs, tissue_name))

    concurrent.futures.wait(futures)

    print('[ {} ] Jobs finished.'.format(datetime.now().strftime('%b %d %H:%M:%S')))
