#!/usr/bin/env python3
"""Verify exact public inputs; optionally download missing files without replacing any."""
import argparse
import json
import urllib.request
from data_utils import ROOT,digest
p=argparse.ArgumentParser();p.add_argument('--download',action='store_true');args=p.parse_args()
for row in json.loads((ROOT/'data/input_manifest.json').read_text()):
    path=ROOT/row['file']
    if not path.exists():
        if not args.download: raise SystemExit('Missing input: '+row['file']+'; run scripts/verify_inputs.py --download')
        if not row['url'].startswith(('https://ftp.ncbi.nlm.nih.gov/','https://data.broadinstitute.org/')):
            raise SystemExit('Unexpected download host')
        path.parent.mkdir(parents=True,exist_ok=True)
        partial=path.with_suffix(path.suffix+'.partial')
        with urllib.request.urlopen(row['url'],timeout=120) as src,partial.open('wb') as dst:
            while block:=src.read(8*1024*1024):dst.write(block)
        if digest(partial)!=row['sha256']:raise SystemExit('Downloaded checksum mismatch; retained partial file')
        partial.rename(path)
    if path.stat().st_size!=row['bytes'] or digest(path)!=row['sha256']:
        raise SystemExit('Input checksum mismatch; source preserved: '+row['file'])
    print('Verified:',row['file'])
manifest=json.loads((ROOT/'data/reference/provenance.json').read_text())
if digest(ROOT/'data/reference/supplementary_donor_metadata.tsv')!=manifest['derived_tsv_sha256']:
    raise SystemExit('Supplementary donor metadata checksum mismatch')
print('Verified: supplementary donor metadata excerpt')
