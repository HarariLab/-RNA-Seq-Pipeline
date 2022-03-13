#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser(description='Get Average Read Lengths from STAR *Log.final.out files')
parser.add_argument('headDir',type=str,help='Path to head STAR output directory. Assumes sample names are set as directory names in this tree.')
#parser.add_argument('query',type=str,help='the label of the line for the corresponding metric (e.g., "Number of splices: Total")')
#parser.add_argument('-s','--sep',type=str, default='.',help='Separator between the sample name and suffix for star output (e.g., "_" in "/sample_Log.final.out"')
args = parser.parse_args()

def cd(path):
    saved_path = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(saved_path)

dirSet = [s for s in os.listdir(args.headDir) if s.startswith('H_VY')]

rows=[]
samples=[]
for (i,d) in enumerate(dirSet):

    lines = []
    samples.append(d.split('_star_v4')[0])
    searchFile = os.path.join(args.headDir,d,samples[i]+'_rRNA_interval.RNA_Metrics')

    with open(searchFile,'r') as fh:
       for (j,row) in enumerate(fh):
           lines.append(row.split())
           if(len(row.split()) <= 0): next
           elif(row.split()[1]=='METRICS'):
               s_idx = j+1
    if(i==0):
        header=[col.strip() for col in lines[s_idx][:-3]]
    rows.append(lines[s_idx+1])

#ids = pd.Series([_id for _id in samples])
df = pd.DataFrame(columns=header, data=rows, index=samples)
df.to_csv(sys.stdout, sep='\t', index=True)

