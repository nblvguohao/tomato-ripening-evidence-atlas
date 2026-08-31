from pathlib import Path
import json,subprocess,sys
import pandas as pd
ROOT=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(ROOT/'scripts/compare_fujisawa_candidates.py'),
                '--inputs',str(ROOT/'data/processed/fujisawa2013'),
                '--candidates',str(ROOT/'revision_expected/rin_candidate_evidence_revised.csv'),
                '--out',str(ROOT/'reproduced_prior')],check=True)
verified=[]
for p in sorted((ROOT/'prior_expected').glob('*.csv')):
    pd.testing.assert_frame_equal(pd.read_csv(p),pd.read_csv(ROOT/'reproduced_prior'/p.name),check_dtype=False)
    verified.append(p.name)
assert len(verified)==4
(ROOT/'reproduced_prior/validation.json').write_text(json.dumps({'table_count':4,'verified':verified,'status':'pass'},indent=2)+'\n')
print('Four prior-comparison tables reproduced exactly.')
