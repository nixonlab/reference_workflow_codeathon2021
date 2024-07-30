#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import pandas as pd

from utils import GTFCOLS, convert_attstr

def get_metadata(args):
    outpre = args.out_prefix
    if not outpre.endswith('/'):
        outpre += '/' if os.path.isdir(outpre) else '.'

    intfeats = ['gene', 'transcript']
    if args.intfeats:
        intfeats = list(pd.Series(intfeats + args.intfeats)
                        .drop_duplicates()
                        )

    print(f'Feature types to process: {", ".join(intfeats)}.', file=sys.stderr)

    allfeats = pd.read_csv(
        args.infile,
        sep='\t',
        comment='#',
        names=GTFCOLS,
        converters={'start':int, 'end':int, 'attributes': convert_attstr},
    )
    print("Loaded GTF", file=sys.stderr)

    # Subset gene, transcript, and other features
    feattypes = {}
    for ftype in intfeats:
        _feats = allfeats[allfeats['feature']==ftype]
        _feats = pd.concat(
            [
                _feats.drop(columns='attributes'),
                _feats['attributes'].apply(pd.Series)
            ],
            axis=1
        )
        feattypes[ftype] = _feats

    del allfeats

    # make output table for other feature types
    ofeats = [f for f in intfeats if f not in ['gene', 'transcript']]
    for ftype in ofeats:
        _outfn = f'{outpre}metadata.{ftype}_features.txt.gz'
        _feats = feattypes.pop(ftype)
        if ftype == 'exon':
            _feats.set_index('exon_id').fillna('').to_csv(_outfn, sep='\t')
        else:
            _feats.fillna('').to_csv(_outfn, sep='\t', index=False)
        del _feats
        print(f"Wrote metadata for '{ftype}' to {_outfn}", file=sys.stderr)

    # Gene features
    genefeats = feattypes.pop('gene')

    # metadata.gene_features
    _outfn = f'{outpre}metadata.gene_features.txt.gz'
    genefeats.set_index('gene_id').fillna('').to_csv(_outfn, sep='\t')
    print(f"Wrote metadata for 'gene' to {_outfn}", file=sys.stderr)

    # gene ID -> gene name
    _outfn = f'{outpre}metadata.gid_gname.txt.gz'
    gid_gname = genefeats[['gene_id', 'gene_name']].drop_duplicates()
    if gid_gname['gene_id'].duplicated().sum():
        sys.stderr.write('WARNING: gene_id is repeated in gid_gname\n')
    gid_gname.dropna().to_csv(_outfn, sep='\t', header=False, index=False)
    print(f"  Wrote 'gid_gname' to {_outfn}", file=sys.stderr)
    del gid_gname

    # gene ID -> gene type
    _outfn = f'{outpre}metadata.gid_gtype.txt.gz'
    gid_gtype = genefeats[['gene_id', 'gene_type']].drop_duplicates()
    if gid_gtype['gene_id'].duplicated().sum():
        sys.stderr.write('WARNING: gene_id is repeated in gid_gtype\n')
    if gid_gtype['gene_type'].isnull().sum():
        sys.stderr.write('WARNING: gene_type is NULL in gid_gtype\n')
    gid_gtype.dropna().to_csv(_outfn, sep='\t', header=False, index=False)
    print(f"  Wrote 'gid_gtype' to {_outfn}", file=sys.stderr)
    del gid_gtype

    # gene ID -> hgnc
    _outfn = f'{outpre}metadata.gid_hgnc.txt.gz'
    gid_hgnc = genefeats[['gene_id', 'hgnc_id']].drop_duplicates()
    if gid_hgnc['gene_id'].duplicated().sum():
        sys.stderr.write('WARNING: gene_id is repeated in gid_hgnc\n')
    # if gid_hgnc['hgnc_id'].isnull().sum():
    #     sys.stderr.write('WARNING: hgnc_id is NULL in gid_hgnc\n')
    gid_hgnc.dropna().to_csv(_outfn, sep='\t', header=False, index=False)
    print(f"  Wrote 'gid_hgnc' to {_outfn}", file=sys.stderr)
    del gid_hgnc

    # Done with genefeats
    del genefeats

    # Transcript features
    _outfn = f'{outpre}metadata.tx_features.txt.gz'
    txfeats = feattypes.pop('transcript')
    txfeats.set_index('transcript_id').fillna('').to_csv(_outfn, sep='\t')
    print(f"Wrote metadata for 'transcript' to {_outfn}", file=sys.stderr)

    # transcript ID -> transcript name
    _outfn = f'{outpre}metadata.tid_tname.txt.gz'
    tid_tname = txfeats[['transcript_id', 'transcript_name']].drop_duplicates()
    if tid_tname['transcript_id'].duplicated().sum():
        sys.stderr.write('WARNING: transcript_id is repeated in tid_tname\n')
    tid_tname.dropna().to_csv(_outfn, sep='\t', header=False, index=False)
    print(f"  Wrote 'tid_tname' to {_outfn}", file=sys.stderr)
    del tid_tname

    # transcript ID -> transcript type
    _outfn = f'{outpre}metadata.tid_ttype.txt.gz'
    tid_ttype = txfeats[['transcript_id', 'transcript_type']].drop_duplicates()
    if tid_ttype['transcript_id'].duplicated().sum():
        sys.stderr.write('WARNING: transcript_id is repeated in tid_ttype\n')
    if tid_ttype['transcript_type'].isnull().sum():
        sys.stderr.write('WARNING: transcript_type is NULL in tid_ttype\n')
    tid_ttype.dropna().to_csv(_outfn, sep='\t', header=False, index=False)
    print(f"  Wrote 'tid_ttype' to {_outfn}", file=sys.stderr)
    del tid_ttype

    # transcript ID -> hgnc
    _outfn = f'{outpre}metadata.tid_hgnc.txt.gz'
    tid_hgnc = txfeats[['transcript_id', 'hgnc_id']].drop_duplicates()
    if tid_hgnc['transcript_id'].duplicated().sum():
        sys.stderr.write('WARNING: transcript_id is repeated in tid_hgnc\n')
    # if tid_hgnc['hgnc_id'].isnull().sum():
    #     sys.stderr.write('WARNING: hgnc_id is NULL in tid_hgnc\n')
    tid_hgnc.dropna().to_csv(_outfn, sep='\t', header=False, index=False)
    print(f"  Wrote 'tid_hgnc' to {_outfn}", file=sys.stderr)
    del tid_hgnc

    # transcript_id to gene_id
    _outfn = f'{outpre}metadata.tid_gid.txt.gz'
    tid_gid = txfeats[['transcript_id', 'gene_id']].drop_duplicates()
    if tid_gid['transcript_id'].duplicated().sum():
        sys.stderr.write('WARNING: transcript_id is repeated in tid_gid\n')
    if tid_gid['gene_id'].isnull().sum():
        sys.stderr.write('WARNING: gene_id is NULL in tid_gid\n')
    tid_gid.set_index('transcript_id').fillna('')\
        .to_csv(_outfn, sep='\t', header=False)
    print(f"  Wrote 'tid_gid' to {_outfn}", file=sys.stderr)

    # gene_id to transcript_id list
    _outfn = f'{outpre}metadata.gid_tid.txt.gz'
    gid_tid = tid_gid.groupby('gene_id')['transcript_id'].apply(list)
    gid_tid.apply(lambda x:','.join(map(str,x)))\
        .to_csv(_outfn, sep='\t', header=False)
    print(f"  Wrote 'gid_tid' to {_outfn}", file=sys.stderr)

    return



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get metadata from GTF')
    parser.add_argument('infile',
                        help="Input GTF."
                        )
    parser.add_argument('out_prefix',
                        help = "Output file prefix."
                        )
    parser.add_argument('--intfeats', action='append',
                        help = '''Create metadata tables for these feature 
                                  types (can specify multiple times). Tables 
                                  for 'gene' and 'transcript' 
                                  features are always enabled, other feature
                                  types may include: ['exon', 'CDS', 'UTR', 
                                  'start_codon', 'stop_codon', 
                                  'Selenocysteine'].
                                  ''')
    args = parser.parse_args()
    get_metadata(args)

## args = parser.parse_args()