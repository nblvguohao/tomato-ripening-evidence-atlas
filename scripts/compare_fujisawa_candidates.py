"""Deterministic factual prior-list comparison from compact, traced inputs."""
from pathlib import Path
import argparse
import collections
import csv
import hashlib
import itertools
import json

ROOT=Path(__file__).resolve().parents[1]
def read(p):
    with p.open(newline='') as f: return list(csv.DictReader(f))
def save(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--inputs',type=Path,default=ROOT/'data/processed/fujisawa2013')
    parser.add_argument('--candidates',type=Path,default=ROOT/'results/hr_revision_20260831/rin_candidate_evidence_revised.csv')
    parser.add_argument('--out',type=Path,default=ROOT/'results/hr_prior_comparison_20260831')
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    prior=read(args.inputs/'dataset2_targets.csv'); binding=read(args.inputs/'dataset1_gene_index.csv')
    peaks=read(args.inputs/'candidate_peak_rows.csv');cand=read(args.candidates)
    assert len(prior)==len({r['gene_id'] for r in prior})==241
    assert collections.Counter(r['original_direction'] for r in prior)=={'positive':137,'negative':104}
    assert len(binding)==len({r['gene_id'] for r in binding})==1200
    assert len(cand)==len({r['gene_id'] for r in cand})==35
    pmap={r['gene_id']:r for r in prior};bmap={r['gene_id']:r for r in binding}
    rows=[]
    for c in cand:
        gene=c['gene_id'];p=pmap.get(gene);b=bmap.get(gene);ps=[r for r in peaks if r['gene_id']==gene]
        reps=sorted({r['replicate'] for r in ps})
        assert len(reps)==int(c['rin_binding_replicates']) and len(reps)>=2
        overlap=any(a['replicate']!=b['replicate'] and a['CHROMOSOME']==b['CHROMOSOME']
                    and max(int(a['PEAK_START']),int(b['PEAK_START']))<=min(int(a['PEAK_END']),int(b['PEAK_END']))
                    for a,b in itertools.combinations(ps,2))
        inferred='lower_in_rin1' if p and p['original_direction']=='positive' else 'higher_in_rin1'
        rows.append({'gene_id':gene,'display_name':c['display_name'],'candidate_grade':c['candidate_grade'],
            'prior_241_membership':'listed' if p else 'not_listed_exact_ID',
            'prior_1200_membership':'listed' if b else 'not_listed_exact_ID',
            'comparison_class':'prior_241' if p else ('prior_binding_only' if b else 'reprocessed_association_only'),
            'prior_direction':p['original_direction'] if p else 'not_available',
            'prior_fc_wt':p['fc_wt'] if p else '', 'prior_fc_rin':p['fc_rin'] if p else '',
            'prior_ecs':p['ecs'] if p else '', 'prior_dataset2_row':p['source_row'] if p else '',
            'prior_dataset1_row':b['source_row'] if b else '',
            'current_rin1_response':c['rin1_response'],
            'cross_allele_sign_compatibility':('compatible' if inferred==c['rin1_response'] else 'different') if p else 'not_assessed',
            'current_binding_replicates':len(reps),'qualifying_peak_rows':len(ps),
            'any_cross_replicate_interval_overlap':overlap,
            'min_abs_feature_to_peak_distance_bp':min(abs(int(r['FEATURE_TO_PEAK_DISTANCE'])) for r in ps),
            'max_abs_feature_to_peak_distance_bp':max(abs(int(r['FEATURE_TO_PEAK_DISTANCE'])) for r in ps),
            'source_url':'https://doi.org/10.1105/tpc.112.108118',
            'mapping_rule':'exact canonical Solyc ID; no gene-model equivalence claim',
            'interpretation':'Prior list membership only; no novelty or causal inference. Gene-level replicate association is not the original site-overlap/window pipeline.'})
    save(args.out/'candidate_prior_comparison.csv',rows)
    augmented=[]
    for c,r in zip(cand,rows):
        augmented.append({**c,**{k:r[k] for k in ['prior_241_membership','prior_1200_membership',
                                                  'comparison_class','any_cross_replicate_interval_overlap']},
                          'occupancy_definition_note':'Frozen gene-level association in >=2 replicates; not the original site/window pipeline. See S17/S18.'})
    save(args.out/'candidate_evidence_with_prior.csv',augmented)
    save(args.out/'candidate_peak_coordinate_audit.csv',[
        {**p,'source_url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+p['replicate'],
         'coordinate_reference':'Original study: tomato genome v2.31, ITAG2; peak-file assembly label not independently specified',
         'coordinate_convention':'Source start/end retained; zero/one-based convention unverified; auxiliary overlap uses inclusive endpoints',
         'coordinate_conversion':'None; not converted to Ensembl Plants release 63'} for p in peaks])
    grade_counts=[{'candidate_grade':g,'candidate_count':sum(r['candidate_grade']==g for r in rows),
                   **{k:sum(r['candidate_grade']==g and r['comparison_class']==k for r in rows)
                      for k in ['prior_241','prior_binding_only','reprocessed_association_only']}} for g in 'ABC']
    save(args.out/'prior_comparison_by_grade.csv',grade_counts)
    manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(args.inputs.glob('*.csv'))}
    summary={'candidates':35,'prior_dataset2_genes':241,'prior_positive':137,'prior_negative':104,
             'prior_dataset1_genes':1200,'classes':dict(collections.Counter(r['comparison_class'] for r in rows)),
             'overlapping_original_directions':dict(collections.Counter(r['prior_direction'] for r in rows if r['prior_241_membership']=='listed')),
             'cross_allele_signs':dict(collections.Counter(r['cross_allele_sign_compatibility'] for r in rows if r['prior_241_membership']=='listed')),
             'without_cross_replicate_interval_overlap':[r['gene_id'] for r in rows if not r['any_cross_replicate_interval_overlap']],
             'by_grade':grade_counts,'input_sha256':manifest,
             'candidate_sha256':hashlib.sha256(args.candidates.read_bytes()).hexdigest(),
             'claim_boundary':'Not-list membership is not novelty. Original spontaneous rin ECS and current engineered rin-1 contrast are not equivalent. Coordinate overlap is an auxiliary diagnostic, not a reconstruction of original peak selection.',
             'frozen_membership_and_grades_unchanged':True}
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if 'sha256' not in k},indent=2))
if __name__=='__main__':main()
