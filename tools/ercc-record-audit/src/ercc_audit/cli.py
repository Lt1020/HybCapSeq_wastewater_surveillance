import argparse
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from . import __version__
from .core import InputError, audit_samples, fit_ercc

SAMPLE_FIELDS=['sample_id','biological_id','replicate','condition','sample_role',
    'rna_ng_ul','dna_ng_ul','library_input_ul','recorded_sample_mass_ng',
    'ercc_volume_ul','ercc_dilution_factor','ercc_lot','total_effective_bases']
STANDARD_FIELDS=['lot_id','ercc_id','length_nt','stock_copies_per_ul']
OBS_FIELDS=['sample_id','ercc_id','aligned_bases','covered_bases','quant_include','exclusion_reason']

def read_tsv(path, required):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        fields=reader.fieldnames or []
        if len(fields)!=len(set(fields)):
            raise InputError(f'{path}: duplicate column names')
        missing=set(required)-set(fields)
        if missing:raise InputError(f'{path}: missing columns {sorted(missing)}')
        rows=[]
        for line,row in enumerate(reader,2):
            if None in row or any(value is None for value in row.values()):
                raise InputError(f'{path}:{line}: invalid column count')
            if any(value.strip() for value in row.values()):rows.append(row)
    if not rows:raise InputError(f'{path}: empty table')
    return rows

def write_tsv(path,rows,fields=None):
    fields=fields or list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t')
        writer.writeheader();writer.writerows(rows)

def destination(path):
    out=Path(path)
    if out.exists() and any(out.iterdir()):
        raise InputError(f'{out}: output directory must be new or empty')
    out.mkdir(parents=True,exist_ok=True)
    return out

def main(argv=None):
    parser=argparse.ArgumentParser(description='ERCC-only calibration and record audit; no pathogen quantification')
    parser.add_argument('--version',action='version',version=__version__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('init',help='write empty input templates');p.add_argument('--outdir',required=True)
    p=sub.add_parser('audit-records',help='validate nucleic-acid input records')
    p.add_argument('--samples',required=True);p.add_argument('--outdir',required=True)
    p=sub.add_parser('fit-ercc',help='fit only synthetic ERCC calibration data')
    for field in ('samples','standards','observations','outdir'):p.add_argument('--'+field,required=True)
    args=parser.parse_args(argv)
    try:
        if args.command=='init':
            out=destination(args.outdir)
            for name,fields in [('samples',SAMPLE_FIELDS),('standards',STANDARD_FIELDS),('observations',OBS_FIELDS)]:
                write_tsv(out/(name+'.tsv'),[],fields)
            print('Created three empty TSV templates.');return 0
        samples=read_tsv(args.samples,SAMPLE_FIELDS)
        audited=audit_samples(samples)
        if args.command=='fit-ercc':
            standards=read_tsv(args.standards,STANDARD_FIELDS)
            observations=read_tsv(args.observations,OBS_FIELDS)
            detail,fits=fit_ercc(samples,standards,observations)
        out=destination(args.outdir)
        write_tsv(out/'sample_audit.tsv',audited)
        if args.command=='fit-ercc':
            write_tsv(out/'ercc_detail.tsv',detail)
            write_tsv(out/'ercc_fits.tsv',fits)
        files=[getattr(args,k) for k in ('samples','standards','observations') if hasattr(args,k)]
        manifest=dict(tool='ercc-record-audit',version=__version__,python=sys.version,
            created_utc=datetime.now(timezone.utc).isoformat(),command=vars(args),
            inputs=[dict(path=str(Path(p).resolve()),sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest()) for p in files],
            scope='synthetic ERCC statistics only; no organism absolute quantities',
            assumptions={'RNA_molecular_mass_g_mol_nt':340,'minimum_distinct_levels':4,
                'automatic_outlier_deletion':False,'LoD_or_LoQ_estimation':False})
        (out/'run_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
        print(f'Wrote {out}');return 0
    except (InputError,OSError,UnicodeError) as e:
        print(f'ERROR: {e}',file=sys.stderr);return 2

