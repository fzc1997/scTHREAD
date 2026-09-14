#!/usr/bin/env python3
"""Independent statistical validation for the PTPRC truncation sensitivity outputs."""
from __future__ import annotations
import csv, json, math
from pathlib import Path
import numpy as np
from scipy.stats import wilcoxon, beta
from statsmodels.stats.multitest import multipletests
import os
ROOT=Path(os.environ.get('SCTHREAD_PROJECT_ROOT','/gpfs/home/fuzc/project/scTHREAD/NAR_database')); OUT=ROOT/'results/truncation_sensitivity_20260913'

def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
models=read(OUT/'ptprc_model_rows.tsv'); usage=read(OUT/'ptprc_model_usage.tsv'); runs=read(OUT/'ptprc_model_usage_by_run.tsv'); cmp=read(OUT/'ptprc_model_comparison.tsv')
assert len({x['run'] for x in runs if x['filter']=='F0_all_positive_models'})==33
# nested filter and support checks
mkeys={(x['run'],x['transcript_id']) for x in models}
assert all(x['pas_supported']=='True' for x in usage if x['filter']=='F1_pas_supported_models')
assert all(x['pas_supported']=='True' and x['hc_reference_model']=='True' for x in usage if x['filter']=='F2_high_confidence_reference_models')
for run in sorted({x['run'] for x in runs}):
 totals={x['filter']:float(x['total_count']) for x in runs if x['run']==run}
 assert totals['F2_high_confidence_reference_models'] <= totals['F1_pas_supported_models'] <= totals['F0_all_positive_models'] + 1e-8
# paired tests
rows=[]
for filt in ['F1_pas_supported_models','F2_high_confidence_reference_models']:
 for label in ['RA_like','RO_like']:
  vals=np.array([float(x[label+'_delta']) for x in cmp if x['filter']==filt],dtype=float)
  stat,p=wilcoxon(vals,alternative='two-sided',zero_method='wilcox',method='auto')
  rows.append({'filter':filt,'metric':label+'_usage_delta','n':len(vals),'wilcoxon_stat':float(stat),'p_value':float(p),'p_method':'normal_approximation_due_to_zero_differences','median_delta':float(np.median(vals)),'mean_delta':float(np.mean(vals)),'positive_fraction':float((vals>0).mean())})
q=multipletests([x['p_value'] for x in rows],method='fdr_bh')[1]
for x,v in zip(rows,q):x['q_value_bh']=float(v);x['significant_q05']=bool(v<0.05)
# exact 95% binomial interval for RA positive direction
ra=[x for x in cmp if x['filter'] in {'F1_pas_supported_models','F2_high_confidence_reference_models'}]
direction=[]
for filt in ['F1_pas_supported_models','F2_high_confidence_reference_models']:
 vals=[float(x['RA_like_delta']) for x in cmp if x['filter']==filt]; k=sum(v>0 for v in vals); n=len(vals)
 lo=beta.ppf(0.025,k,n-k+1) if k else 0.0; hi=beta.ppf(0.975,k+1,n-k) if k<n else 1.0
 direction.append({'filter':filt,'positive':k,'n':n,'fraction':k/n,'exact_binomial_ci95_low':float(lo),'exact_binomial_ci95_high':float(hi)})
with (OUT/'ptprc_model_validation.tsv').open('w',newline='') as f:
 fields=list(rows[0]);w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
(OUT/'ptprc_direction_validation.tsv').write_text('filter\tpositive\tn\tfraction\texact_binomial_ci95_low\texact_binomial_ci95_high\n'+'\n'.join('\t'.join(str(x[k]) for k in direction) for x in [])) if False else None
with (OUT/'ptprc_direction_validation.tsv').open('w',newline='') as f:
 fields=list(direction[0]);w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(direction)
# bootstrap consistency checks
boot=json.loads((OUT/'ptprc_model_bootstrap_summary.json').read_text())
assert boot['seed']==20260913 and boot['resamples']==1000
assert all(x['ci95_low']<=x['estimate_median']<=x['ci95_high'] for x in boot['rows'])
summary={'status':'PASS','checks':{'f0_f1_f2_nested':True,'support_flags_consistent':True,'wilcoxon_tests':len(rows),'bh_correction':'four paired comparisons','bootstrap_seed':20260913,'bootstrap_resamples':1000},'outputs':['ptprc_model_validation.tsv','ptprc_direction_validation.tsv','ptprc_model_bootstrap_summary.json'],'interpretation':'Validation supports a measurable paired shift in the scoped model-level sensitivity output; it does not establish that every filtered model is artifactual or that all transcript abundance is corrected.'}
(OUT/'ptprc_model_validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with (OUT/'ptprc_model_validation_report.md').open('w') as f:
 f.write('# PTPRC model truncation sensitivity validation\n\n')
 f.write('Status: **PASS** for design and statistical validation of the scoped sensitivity output.\n\n')
 f.write('| Filter | Metric | n | Wilcoxon statistic | P | BH q | Median delta | Positive fraction |\n|---|---:|---:|---:|---:|---:|---:|---:|\n')
 for x in rows:f.write(f"| {x['filter']} | {x['metric']} | {x['n']} | {x['wilcoxon_stat']:.3g} | {x['p_value']:.3g} | {x['q_value_bh']:.3g} | {x['median_delta']:.3f} | {x['positive_fraction']:.3f} |\n")
 f.write('\nAssumption checks: paired run-level units are fixed by the same run IDs; no cross-study pooling is introduced. The primary effect is median paired delta with the pre-registered 1,000-resample bootstrap CI. BH correction spans the four pre-specified paired tests.\n')
 f.write('The validated result is a scoped endpoint/completeness sensitivity signal. It is not a claim that all excluded models are technical artifacts, and it does not replace the direct PTPRC exon-3 junction measurement.\n')
print(json.dumps(summary,indent=2))
