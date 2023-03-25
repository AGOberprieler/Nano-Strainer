#!/usr/bin/env python

# Author: Tankred Ott
# Basic usage in simple mode (input: one FASTA file): 7_vsearch_testclusteringthresholds.py simple /[path/to/inputfile].fasta /[path/to/outputfolder]/ [all CTs to be tested, separated by spaces, e.g.: --ct 0.50 --ct 0.75] -t [number of threads]
# Basic usage in multi mode (input: one directory of FASTA files): 7_vsearch_testclusteringthresholds.py multi "/[path/to/inputfolder]/*.fasta" /[path/to/outputfolder]/ [all CTs to be tested, separated by spaces, e.g.: --ct 0.50 --ct 0.75] -t [number of threads]
# NOTE: For the input path in multi mode, the double quotation marks are required
# NOTE: VSEARCH is called in line 95 ("VSEARCH_CALL"); replace "/path/to/vsearch" with your own path of the VSEARCH application. Parameters for the VSEARCH clustering process can be modified by altering the command

import argparse
import os
import glob

parser = argparse.ArgumentParser('Run "vsearch --cluster-fast" for multiple clustering thresholds (ids) and/or multiple inputs.')
subparsers = parser.add_subparsers(dest='subcommand')
simple_parser = subparsers.add_parser(
    'simple',
    help='Run "vsearch --cluster-fast" for multiple clustering thresholds (ids).'
)
simple_parser.add_argument(
    'in_file',
    type=str,
    help='Path to input fasta file containing sequences to be clustered.'
)
simple_parser.add_argument(
    'out_dir',
    type=str,
    help='''
        Path outputs will be written to, will be generated if it does not exist.
        For each clustering threshold (id), a single directory within
        out_dir will be generated.
    '''
)
simple_parser.add_argument(
    '-c', '--ct',
    action='append',
    type=float,
    help='''
        Clustering thresholds (ids).
    ''',
    default=[]
)
simple_parser.add_argument(
    '-t', '--threads',
    type=int,
    help='''
        Number of threads
    ''',
    default=1
)

multi_parser = subparsers.add_parser(
    'multi',
    help='''
        Run "vsearch --cluster-fast" for multiple inputs and multiple clustering thresholds (ids).

        Inputs are passed as shell-expandable string (e.g. "*fasta").
        Outputs will be written to subdirectories within the output directory, where each subdirectory
        contains the vsearch results for a single input. The output directories are named according to
        the base names of the corresponding inputs sans file extension (e.g. my_file.fasta outputs will
        be written to <out_dir>/my_file/).
    '''
)
multi_parser.add_argument(
    'in_files',
    type=str,
    help='Shell-expandable string, e.g. "*fasta". ATTENTION: the double quotes are required!'
)
multi_parser.add_argument(
    'out_dir',
    type=str,
    help='''
        Path outputs will be written to, will be generated if it does not exist.
        For each input, a single directory within out_dir will be generated.
    '''
)
multi_parser.add_argument(
    '-c', '--ct',
    action='append',
    type=float,
    help='''
        Clustering thresholds (ids).
    ''',
    default=[]
)
multi_parser.add_argument(
    '-t', '--threads',
    type=int,
    help='''
        Number of threads
    ''',
    default=1
)

VSEARCH_CALL = '/path/to/vsearch --cluster_fast {} --centroids {} --clusters {} --msaout {} --uc {} --blast6out {} --id {} --log {} --threads {} --clusterout_id --clusterout_sort --consout {} --qmask none --strand both'

def simple_subcommand(in_file, out_dir, cts, threads):
    if not os.path.isfile(in_file):
        err = 'File "{}" does not exist or is a directory.'.format(in_file)
        raise FileNotFoundError(err)
    
    if len(cts) < 1:
        print('Please provide at least one clustering threshold (-c/--ct)')

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    for ct in cts:
        ct_str = str(round(ct * 100))

        out_dir_cur = os.path.join(out_dir, ct_str)
        if not os.path.exists(out_dir_cur):
            os.mkdir(out_dir_cur)

        vsearch_call = VSEARCH_CALL.format(
            in_file,
            os.path.join(out_dir_cur, 'centroids.fasta'),
            os.path.join(out_dir_cur, 'cluster'),
            os.path.join(out_dir_cur, 'clusteralignments.fasta'),
            os.path.join(out_dir_cur, 'uclust.txt'),
            os.path.join(out_dir_cur, 'blast6.txt'),
            ct,
            os.path.join(out_dir_cur, 'log.txt'),
            threads,
            os.path.join(out_dir_cur, 'clusterconsensus.fasta'),
        )

        os.system(vsearch_call)

def multi_subcommand(in_files, out_dir, cts, threads):
    in_file_paths = glob.glob(in_files)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    for in_file_path in in_file_paths:
        base_name = os.path.basename(in_file_path)
        base_name = os.path.splitext(base_name)[0]
        out_dir_cur = os.path.join(out_dir, base_name)
        simple_subcommand(
            in_file=in_file_path,
            out_dir=out_dir_cur,
            cts=cts,
            threads=threads
        )

if __name__ == "__main__":
    args = parser.parse_args()
    # print(args)

    subcommand = args.subcommand
    out_dir = args.out_dir
    cts = args.ct
    threads = args.threads

    if subcommand == 'simple':
        in_file = args.in_file
        simple_subcommand(
            in_file=in_file,
            out_dir=out_dir,
            cts=cts,
            threads=threads,
        )
    elif subcommand == 'multi':
        in_files = args.in_files
        multi_subcommand(
            in_files=in_files,
            out_dir=out_dir,
            cts=cts,
            threads=threads,
        )

    