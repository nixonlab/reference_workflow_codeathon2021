#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO


def gbk2fna(args):
    _ = SeqIO.write(
        SeqIO.parse(args.infile, 'genbank'),
        args.outfile,
        "fasta"
    )

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Convert Genbank flat file to FASTA file')
    parser.add_argument('infile',
                        nargs='?',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        help="Input Genbank flat file. Default: stdin."
                        )
    parser.add_argument('outfile',
                        nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="Output FASTA file. Default: stdout."
                        )
    gbk2fna(parser.parse_args())
