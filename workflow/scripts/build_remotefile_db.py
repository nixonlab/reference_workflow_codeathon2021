#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re
import requests
import json
import yaml

def gencode_url(gencode_release):
    base_url = 'http://ftp.ebi.ac.uk/pub/databases/gencode'
    if re.match('^M\d+', gencode_release):
        return f'{base_url}/Gencode_mouse/release_{gencode_release}'
    elif re.match('\d+', gencode_release):
        return f'{base_url}/Gencode_human/release_{gencode_release}'

def build_remote_db(
        remote_files: dict[str, dict],
        remote_filesets: dict[str, dict],
        gencode_versions: list[str],
):

    db = {}

    """ Individual remote files """
    db = {**db, **remote_files}

    """ Remote file sets """
    for n, d in remote_filesets.items():
        print(f'Adding remote information for {n}', file=sys.stderr)
        _baseurl = d['baseurl']
        r = requests.get(d['checksum_url'])
        _iter = re.findall(r"(?P<_md5>[a-fA-F\d]{32})\s+(?P<_fn>\S+)", r.text)
        for _md5, _fn in _iter:
            _fn0 = re.sub(r"^\./", "", _fn) # removed leading dotslash
            _fnL = re.sub(r"/", ".", _fn0)  # replace slash with dot
            assert _fnL not in db, f'ERROR: filename conflict "{_fnL}"'
            db[_fnL] = {
                'url': f'{_baseurl}/{_fn0}',
                'md5': _md5,
            }

    """ Add GENCODE versions """
    for gr in gencode_versions:
        print(f'Adding remote information for GENCODE v{gr}', file=sys.stderr)
        _baseurl = gencode_url(gr)
        r = requests.get(f'{_baseurl}/MD5SUMS')
        _iter = (row.split() for row in r.text.split('\n') if row)
        for _md5, _fn in _iter:
            db[_fn] = {
                'url': f'{_baseurl}/{_fn}',
                'md5': _md5,
            }

    return db


def main(args):
    data = yaml.safe_load(args.remote_yaml)
    outdb = build_remote_db(
        data['remote_files'],
        data['remote_filesets'],
        data['gencode_versions']
    )
    json.dump(outdb, args.remotedb_json, indent=4)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Build remote file database')
    parser.add_argument('remote_yaml',
                        nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin,
                        help="Input YAML file. Default: stdin.")
    parser.add_argument('remotedb_json',
                        nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="Output JSON file. Default: stdout.")
    main(parser.parse_args())
