#!/usr/bin/env python3
"""PTPRC model-level truncation sensitivity using real IsoQuant tables and PolyASite 2.0.

This is a scoped sensitivity analysis, not a claim that endpoint support proves a
full-length molecule. It compares all positive PTPRC reference/model counts with
nested model-end filters and records missing/low-support states explicitly.
"""
from __future__ import annotations
import csv, gzip, hashlib, json, math, re, time
from collections import defaultdict
from pathlib import Path
import numpy as np

import os
ROOT=Path(os.environ.get('SCTHREAD_PROJECT_ROOT','/gpfs/home/fuzc/project/scTHREAD/NAR_database'))
CLASS=ROOT/'docs/processing_logic_20260909/run_classification.tsv'
POLYA=Path(os.environ.get('SCTHREAD_POLYASITE_ATLAS','/gpfs/home/fuzc/project/Genome_reference/endpoint_atlas/hg38/atlas.clusters.2.0.GRCh38.96.bed.gz'))
OUT=ROOT/'results/truncation_sensitivity_20260913'
OUT.mkdir(parents=True,exist_ok=True)
TARGET_GENE='ENSG00000081237'  # human PTPRC in the registered GRCh38 annotation
TARGET_NAME='PTPRC'


def sha(path:Path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for c in iter(lambda:f.read(1024*1024),b''):h.update(c)
 return h.hexdigest()

def attrs(text:str):
 out={}
 for k,v in re.findall(r'(\w+)\s+"([^"]+)"',text): out[k]=v
 # tags are repeated and are parsed separately below
 out['tags']=re.findall(r'tag\s+"([^"]+)"',text)
 return out

def load_pas(path:Path):
 by={}
 with gzip.open(path,'rt') as f:
  for line in f:
   if not line or line.startswith('#'): continue
   x=line.rstrip('\n').split('\t')
   if len(x)<6: continue
   chrom=x[0] if x[0].startswith('chr') else 'chr'+x[0]
   try: start=int(x[1]); end=int(x[2])
   except ValueError: continue
   by.setdefault((chrom,x[5]),[]).append((start,end,x[3]))
 for key in by: by[key].sort()
 return by

def pas_match(pas, chrom, strand, pos0, win=25):
 candidates=pas.get((chrom,strand),())
 # small linear scan over nearby sorted entries; PTPRC has few clusters
 for start,end,pid in candidates:
  if start > pos0 + win: break
  if end >= pos0 - win and start <= pos0 + win: return True,pid
 return False,''

def parse_reference_ptprc_ids(gtf:Path):
 ids=set(); opener=gzip.open if gtf.suffix=='.gz' else open
 with opener(gtf,'rt') as f:
  for line in f:
   if not line or line.startswith('#'): continue
   x=line.rstrip('\n').split('\t')
   if len(x)<9 or x[2] != 'transcript': continue
   a=attrs(x[8])
   if a.get('gene_id') == TARGET_GENE and a.get('transcript_id'): ids.add(a['transcript_id'].split('.',1)[0])
 return ids

def parse_ptprc_models(gtf:Path,pas):
 models={}; current={}
 opener=gzip.open if gtf.suffix=='.gz' else open
 with opener(gtf,'rt') as f:
  for line in f:
   if not line or line.startswith('#'): continue
   x=line.rstrip('\n').split('\t')
   if len(x)<9 or x[2] not in {'transcript','exon'}: continue
   a=attrs(x[8]); gid=a.get('gene_id','')
   if gid!=TARGET_GENE: continue
   tid=a.get('transcript_id','')
   if not tid: continue
   rec=models.setdefault(tid,{'transcript_id':tid,'gene_id':gid,'chrom':x[0],'strand':x[6],'start':None,'end':None,'source':x[1],'tags':set(),'tsl':''})
   try: st=int(x[3]); en=int(x[4])
   except ValueError: continue
   rec['start']=st if rec['start'] is None else min(rec['start'],st)
   rec['end']=en if rec['end'] is None else max(rec['end'],en)
   rec['tags'].update(a.get('tags',[]))
   if 'transcript_support_level' in a: rec['tsl']=a['transcript_support_level']
 out={}
 for tid,m in models.items():
  # GTF is 1-based inclusive; convert terminal coordinate to 0-based for BED comparison.
  terminal=(m['end']-1) if m['strand']=='+' else (m['start']-1)
  ok,pid=pas_match(pas,m['chrom'],m['strand'],terminal,25)
  tsl=m['tsl'].split()[0] if m['tsl'] else ''
  ref=tid.startswith('ENST')
  hc_ref=ref and (tsl in {'1','2'} or 'basic' in m['tags']) and not ({'mRNA_end_NF','cds_end_NF'} & m['tags'])
  m.update({'terminal_0based':terminal,'pas_supported':ok,'pas_id':pid,'is_reference':ref,'hc_reference_model':hc_ref,'tags':'|'.join(sorted(m['tags']))})
  out[tid]=m
 return out

def parse_counts(path:Path):
 vals={}
 with path.open() as f:
  rd=csv.DictReader(f,delimiter='\t')
  for r in rd:
   tid=r.get('feature_id','').split('.',1)[0]
   try: v=float(r.get('count','0'))
   except (TypeError,ValueError): continue
   if v>0: vals[tid]=vals.get(tid,0.0)+v
 return vals

with CLASS.open(newline='') as f:
 runs=[r for r in csv.DictReader(f,delimiter='\t') if r.get('study') in {'GSE276974','GSE307660'}]
pas=load_pas(POLYA)
reference_gtf_candidates=[
 Path(os.environ.get('SCTHREAD_REFERENCE_ROOT','/gpfs/home/fuzc/project/Genome_reference')+'/10x_ref/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'),
 ROOT/'references/10x_refdata-gex-GRCh38-2024-A.genes.gtf.gz',
]
reference_gtf=next((x for x in reference_gtf_candidates if x.is_file()), None)
if reference_gtf is None:
 raise SystemExit('MISSING_REGISTERED_GRCH38_REFERENCE_GTF')
reference_ptprc_ids=parse_reference_ptprc_ids(reference_gtf)
model_rows=[]; usage_rows=[]; run_rows=[]; missing=[]
for r in runs:
 run=r['run']; outdir=Path(r['isoquant_output_dir']); stem=r['isoquant_output_stem'] or run
 gtf=Path(r['isoquant_cli_source'])
 counts_path=outdir/f'{stem}.transcript_counts.tsv'
 if not gtf.is_file() or not counts_path.is_file():
  missing.append({'run':run,'study':r['study'],'gtf':str(gtf),'counts':str(counts_path),'reason':'missing_model_or_count_file'}); continue
 models=parse_ptprc_models(gtf,pas)
 counts=parse_counts(counts_path)
 # Keep PTPRC transcripts/models in the model GTF. Counts not represented in the model GTF are retained as unmapped audit rows.
 target={tid:v for tid,v in counts.items() if tid in models}
 unmapped={tid:v for tid,v in counts.items() if tid in reference_ptprc_ids and tid not in models}
 if not target:
  missing.append({'run':run,'study':r['study'],'gtf':str(gtf),'counts':str(counts_path),'reason':'no_positive_ptprc_models'}); continue
 for tid,m in models.items():
  if tid not in target: continue
  base={'run':run,'study':r['study'],'transcript_id':tid,'count':target[tid],'is_reference':m['is_reference'],'hc_reference_model':m['hc_reference_model'],'pas_supported':m['pas_supported'],'pas_id':m['pas_id'],'terminal_0based':m['terminal_0based'],'tags':m['tags'],'tsl':m['tsl']}
  model_rows.append(base)
 filters={
  'F0_all_positive_models':lambda m: True,
  'F1_pas_supported_models':lambda m: bool(m['pas_supported']),
  'F2_high_confidence_reference_models':lambda m: bool(m['hc_reference_model'] and m['pas_supported']),
  'F3_pas_supported_novel_models':lambda m: bool((not m['is_reference']) and m['pas_supported'] and target.get(m['transcript_id'],0)>=3),
 }
 for fname,keep in filters.items():
  selected=[(tid,m,target[tid]) for tid,m in models.items() if tid in target and keep(m)]
  total=sum(v for _,_,v in selected)
  for tid,m,v in selected:
   usage_rows.append({'run':run,'study':r['study'],'filter':fname,'transcript_id':tid,'count':v,'usage':(v/total if total else None),'pas_supported':m['pas_supported'],'is_reference':m['is_reference'],'hc_reference_model':m['hc_reference_model']})
  run_rows.append({'run':run,'study':r['study'],'filter':fname,'n_positive_models':len(selected),'total_count':total,'n_unmapped_positive_models':len(unmapped),'unmapped_positive_count':sum(unmapped.values()),'scorable':bool(total>0)})
for row in model_rows:
 row['model_support_class']='reference_hc' if row['hc_reference_model'] else ('reference_other' if row['is_reference'] else 'novel_run_local')
fields=list(run_rows[0]) if run_rows else ['run','study','filter','n_positive_models','total_count','n_unmapped_positive_models','unmapped_positive_count','scorable']
with (OUT/'ptprc_model_usage_by_run.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(run_rows)
with (OUT/'ptprc_model_rows.tsv').open('w',newline='') as f:
 fields=list(model_rows[0]) if model_rows else ['run','study','transcript_id','count','is_reference','hc_reference_model','pas_supported','pas_id','terminal_0based','tags','tsl','model_support_class']
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(model_rows)
with (OUT/'ptprc_model_usage.tsv').open('w',newline='') as f:
 fields=list(usage_rows[0]) if usage_rows else ['run','study','filter','transcript_id','count','usage','pas_supported','is_reference','hc_reference_model']
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(usage_rows)
with (OUT/'ptprc_model_missing.tsv').open('w',newline='') as f:
 fields=list(missing[0]) if missing else ['run','study','gtf','counts','reason']
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(missing)
summary={'generated_at':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'script':str(ROOT/'scripts/run_ptprc_model_truncation_sensitivity_20260913.py'),'target_gene':TARGET_GENE,'target_name':TARGET_NAME,'studies':['GSE276974','GSE307660'],'runs_in_scope':len(runs),'runs_scorable':len({x['run'] for x in run_rows}),'missing_or_unscorable':len(missing),'filters':list(filters),'polyasite_path':str(POLYA),'polyasite_sha256':sha(POLYA),'window_bp_each_side':25,'gtf_coordinate_conversion':'GTF 1-based inclusive terminal coordinate converted to 0-based for BED comparison','model_level_note':'PAS support is endpoint plausibility, not proof of full-length molecule. F2 is a high-confidence reference-model subset.'}
(OUT/'ptprc_model_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(summary,indent=2,ensure_ascii=False))
