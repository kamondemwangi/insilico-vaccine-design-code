
from __future__ import division
import pandas as pd
import random
import os
import sys
import csv
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

import time
import argparse
from os import path
from datetime import datetime
import itertools

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def get_arguments(tcl_file, bcl_file, r=10):
    """
    fetch input parameters from the argument parser

    parameters:
        - tcl_file (csv file): file having the T-cell epitopes
        - bcl_file (csv file): file having the B-cell epitopes
        - r (int): number of random epitopes
    returns:
        - params (dict): dictionary with random size as key and values as
        input files
    """
    params = dict()
    params[r] = [os.path.abspath(tcl_file), os.path.abspath(bcl_file)]
    return params


def generate_combinations(params, tcell_linker='AAY', bcell_linker='GPGPG'):

    """
    generate a specified number of random combinations of epitopes given
    input files (both T-cell and B-cell) epitopes files

    parameters:
        - params (dict): dictionary with random size as key and values as
        input files
        - tcell_linker (str): string characters for linking the T-cell epitopes
        - bcell_linker (str): string characters for linking the B-cell epitopes

    returns:
        - outdir (dir): directory to write the outputs
    """

    tcl_dict = dict() # dictionary to store random t-cell epitope combinations
    bcl_dict = dict() # dictionary to store random b-cell epitope combinations
    epitope_dict = dict()
    epitopes_list = list()

    # read the input files and extract epitopes as a list
    for randomsize, files in params.items():
        for f in files:
            data = pd.read_csv(f)
            epitopes = data.iloc[:,1].tolist()
            epi_len = int(len(''.join(epitopes))/len(epitopes)) # length of each epitope
            celltype = os.path.splitext(os.path.basename(f))[0]
            epitope_dict[celltype+'_'+str(epi_len)] = epitopes


    # randomly shuffle the epitope sequences using the random size
    headers_list = list()
    for ct, ep in epitope_dict.items():
        for i in range(randomsize):
            non_random_seq = ''.join(ep)
            random.shuffle(ep)
            epitope = ''.join(ep)
            epitopes_list.append(epitope)
            headers_list.append(ct)

    tcell_epitopes = epitopes_list[0:randomsize]
    tcell_headers = headers_list[0:randomsize]
    bcell_epitopes = epitopes_list[randomsize:]
    bcell_headers = headers_list[randomsize:]

    for i, j in enumerate(zip(tcell_headers, tcell_epitopes)):
        seqid = '>'+''+j[0]+'_'+str(i+1)
        epitope = j[1]
        if not epitope in tcl_dict:
            tcl_dict[epitope] = [seqid]
        else:
            tcl_dict[epitope].append(seqid)


    # create output directory
    tcell_bs = j[0].rsplit('_',2)[0]
    outdir = os.path.join(BASE_DIR, 'output')
    if not os.path.exists(outdir):
        os.makedirs(outdir)



    tcell_out = os.path.join(outdir, tcell_bs)+'_raw_shuffled_epitopes.fa'
    tcell_linked = os.path.join(outdir, tcell_bs)+'_linked_shuffled_epitopes.fa'
    tcell_list = list()
    with open(tcell_out, 'w') as f_obj, open(tcell_linked, 'w') as l_obj:
        for seq, id in tcl_dict.items():
            epi_len = int(id[0].rsplit('_',2)[1])
            sub = list(map(''.join, zip(*[iter(seq)]*epi_len)))
            linked_seq = "{}".format(tcell_linker).join(sub)
            tcell_list.append(linked_seq)
            print("writing {} epitopes to {}".format(''.join(id), os.path.basename(tcell_out)))
            f_obj.write(''.join(id))
            f_obj.write('\n')
            if len(seq) > 1:
                f_obj.write('\n'.join(seq[i:i+60] for i in range(0,len(seq), 60)))
                f_obj.write('\n')
            else:
                f_obj.write(seq)
                f_obj.write('\n')

            # write linked epitopes
            print("\nwriting {} epitopes to {}".format(''.join(id), os.path.basename(tcell_linked)))
            l_obj.write(''.join(id))
            l_obj.write('\n')
            if len(linked_seq) > 1:
                l_obj.write('\n'.join(linked_seq[i:i+60] for i in range(0,len(linked_seq), 60)))
                l_obj.write('\n')
            else:
                l_obj.write(linked_seq)
                l_obj.write('\n')



    for i, j in enumerate(zip(bcell_headers, bcell_epitopes)):
        seqid = '>'+''+j[0]+'_'+str(i+1)
        epitope = j[1]
        if not epitope in bcl_dict:
            bcl_dict[epitope] = [seqid]
        else:
            bcl_dict[epitope].append(seqid)

    bcell_bs = j[0].rsplit('_',2)[0]
    bcell_out = os.path.join(outdir, bcell_bs)+'_raw_shuffled_epitopes.fa'
    bcell_linked = os.path.join(outdir, bcell_bs)+'_linked_shuffled_epitopes.fa'
    bcell_list = list()
    with open(bcell_out, 'w') as f_obj, open(bcell_linked, 'w') as l_obj:
        for seq, id in bcl_dict.items():
            epi_len = int(id[0].rsplit('_',2)[1])
            sub = list(map(''.join, zip(*[iter(seq)]*epi_len)))
            linked_seq = "{}".format(bcell_linker).join(sub)
            bcell_list.append(linked_seq)
            print("writing {} epitopes to {}".format(''.join(id), os.path.basename(bcell_out)))
            f_obj.write(''.join(id))
            f_obj.write('\n')
            if len(seq) > 1:
                f_obj.write('\n'.join(seq[i:i+60] for i in range(0,len(seq), 60)))
                f_obj.write('\n')
            else:
                f_obj.write(seq)
                f_obj.write('\n')

            # write linked epitopes
            print("\nwriting {} epitopes to {}".format(''.join(id), os.path.basename(bcell_linked)))
            l_obj.write(''.join(id))
            l_obj.write('\n')
            if len(linked_seq) > 1:
                l_obj.write('\n'.join(linked_seq[i:i+60] for i in range(0,len(linked_seq), 60)))
                l_obj.write('\n')
            else:
                l_obj.write(linked_seq)
                l_obj.write('\n')



    # concatenate the B-cell and T-cell epitopes using linkers

    linker = 'AKFVAAWTLKAAAEAAAK'
    out = os.path.join(outdir, 'btcell_epitopes.fa')
    with open(out, 'w') as f_obj:
        for i, epitope in enumerate(zip(bcell_list, tcell_list)):
            linked_btcell = "{}".format(linker)+epitope[0]+"{}".format(tcell_linker)+(epitope[1])
            seqid = '>'+str(i+1)
            f_obj.write(seqid)
            f_obj.write('\n')
            if len(linked_btcell) > 1:
                f_obj.write('\n'.join(linked_btcell[i:i+60] for i in range(0,len(linked_btcell), 60)))
                f_obj.write('\n')
            else:
                f_obj.write(linked_btcell)
                f_obj.write('\n')




# argument parsers / command line options
parser=argparse.ArgumentParser()
helpstr = """python generate_epitope_combinations.py [options]"""
required_group = parser.add_argument_group('required arguments')
required_group.add_argument('-t', '--tcell', help="T-cell epitope file path")
required_group.add_argument('-b', '--bcell', help="B-cell epitope file path")
parser.add_argument('-r', '--randomsize', default=10, type=int, help="number of randomized epitopes")
parser.add_argument('-tcl', '--tclinker', default='AAY', type=str, help="T-cell epitopes linker")
parser.add_argument('-bcl', '--bclinker', default='GPGPG', type=str, help="B-cell epitopes linker")
args=parser.parse_args()

# open input file:
if args.tcell != None:
    tcell_f = args.tcell
else:
    print("Please specify the  T-cell epitopes file!\n")
    sys.exit(2)

if args.bcell != None:
    bcell_f = args.bcell
else:
    print("Please specify the  B-cell epitopes file!\n")
    sys.exit(2)

if args.randomsize != None and type(args.randomsize) == int:
    size = args.randomsize
else:
    print("specify the random size number\n")

if args.tclinker != None and type(args.tclinker) == str:
    tlinker = args.tclinker.upper()
else:
    print("specify the T-cell epitopes linker\n")

if args.bclinker != None and type(args.bclinker) == str:
    blinker = args.bclinker.upper()
else:
    print("specify the B-cell epitopes linker\n")

# run the function
fs = get_arguments(tcell_f, bcell_f, size)
g = generate_combinations(fs, tlinker, blinker)
