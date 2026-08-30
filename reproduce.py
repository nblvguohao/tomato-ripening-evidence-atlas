#!/usr/bin/env python3
"""Reproduce frozen tables without raw inputs or manuscript files."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'reproduced'

def run(script, *args):
    subprocess.run([sys.executable, str(ROOT/'src'/script), *map(str,args)], cwd=ROOT, check=True)

def main():
    for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines():
        checksum, relative = line.split('  ', 1)
        path = (ROOT/relative).resolve()
        if not path.is_relative_to(ROOT) or hashlib.sha256(path.read_bytes()).hexdigest() != checksum:
            raise ValueError(f'Hash mismatch: {relative}')
    subprocess.run([sys.executable,'-m','unittest','discover','-s','tests','-p','test_*.py','-q'],cwd=ROOT,check=True)
    OUT.mkdir(exist_ok=True)
    run('export_concordant_signature.py',
        '--study-a-matrix','data/processed/GSE42783_rank01.csv.gz',
        '--study-a-metadata','data/processed/GSE42783_sample_metadata.csv',
        '--study-b-matrix','data/processed/GSE108415_rank01.csv.gz',
        '--study-b-metadata','data/processed/GSE108415_sample_metadata.csv',
        '--top-k',500,'--output',OUT/'signature.csv')
    run('build_expression_consensus.py','--signature',OUT/'signature.csv',
        '--validation-a-matrix','data/processed/GSE267238_log2cpm.csv.gz',
        '--validation-a-metadata','data/processed/GSE267238_sample_metadata.csv',
        '--validation-a-genotype','WT','--validation-a-early',0,'--validation-a-late',1,
        '--validation-b-matrix','data/processed/GSE128739_log2cpm.csv.gz',
        '--validation-b-metadata','data/processed/GSE128739_sample_metadata.csv',
        '--validation-b-genotype','Moneymaker','--validation-b-early',0,'--validation-b-late',3,
        '--gene-metadata','data/processed/GSE267238_gene_metadata.csv.gz',
        '--output',OUT/'consensus.csv')
    run('build_evidence_graph.py','--consensus',OUT/'consensus.csv',
        '--binding-seeds','config/published_binding_edges_v2.csv','--output-dir',OUT)
    run('build_ripening_atlas_resource.py','--consensus',OUT/'consensus.csv',
        '--annotation','results/signature_annotation_ensembl_plants63.csv',
        '--edges',OUT/'evidence_graph_edges.csv',
        '--studyshield','results/v3/studyshield_scores.csv.gz',
        '--studyshield-significance','results/v3/studyshield_permutation_significance.csv',
        '--go-enrichment','results/go_enrichment_propagated_reduced.csv',
        '--go-annotations','results/ITAG4_go_propagated.csv.gz',
        '--protein-support','results/PXD051570_signature_support.csv',
        '--tissue-support','results/SRP109982_signature_tissue_support.csv',
        '--perturbation-support','results/GSE210589/deseq2/GSE210589_RIN_binding_perturbation_evidence.csv',
        '--output',OUT/'atlas.csv')
    sys.path.insert(0,str(ROOT/'src'))
    from build_hr_submission_bundle import build_candidates, read_csv, write_csv
    candidates=build_candidates(read_csv(OUT/'atlas.csv'),
        read_csv(ROOT/'results/GSE210589/deseq2/GSE210589_RIN_binding_perturbation_evidence.csv'),
        read_csv(ROOT/'config/grade_a_literature_screen.csv'))
    write_csv(OUT/'candidates.csv',candidates,list(candidates[0]))
    import pandas as pd
    checks=[]
    pairs=[('signature.csv','concordant_mature_green_to_ripe_signature.csv'),
           ('consensus.csv','expression_consensus_4cohort.csv'),
           ('atlas.csv','ripening_atlas_resource_v3.csv'),
           ('candidates.csv','hr_submission/rin_candidate_evidence_table.csv')]
    for regenerated, frozen in pairs:
        actual=pd.read_csv(OUT/regenerated).sort_values('gene_id').reset_index(drop=True)
        expected=pd.read_csv(ROOT/'results'/frozen).sort_values('gene_id').reset_index(drop=True)
        pd.testing.assert_frame_equal(actual,expected,check_dtype=False,check_exact=False,rtol=1e-10,atol=1e-12)
        checks.append({'table':regenerated,'rows':len(actual),'comparison':'pass','rtol':1e-10,'atol':1e-12})
    report={'scope':'processed-input table reproduction; upstream frozen intermediates retained',
            'comparisons':checks,'candidate_count':len(candidates),
            'grade_a_count':sum(x['candidate_grade']=='A' for x in candidates)}
    (OUT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__ == '__main__':
    main()
