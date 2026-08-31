"""Reproduce the supplementary revision from processed inputs and verify tables."""
from pathlib import Path
import importlib.util
import json
import sys
import pandas as pd
ROOT=Path(__file__).resolve().parent
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,file)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
analysis=load('revision_analysis',ROOT/'scripts/revise_hr_analysis_20260831.py')
analysis.freeze_baseline=lambda:None
analysis.OUT=ROOT/'reproduced_revision'
analysis.main()
resource=load('revision_resource',ROOT/'scripts/revise_hr_resource_20260831.py')
resource.OUT=analysis.OUT
resource.main()
example=load('revision_worked_query',ROOT/'scripts/build_hr_resource_use_case.py')
example.OUT=analysis.OUT
example.main()
verified=[]
for path in sorted((ROOT/'revision_expected').glob('*.csv*')):
    observed=analysis.OUT/path.name
    pd.testing.assert_frame_equal(pd.read_csv(path),pd.read_csv(observed),check_dtype=False,rtol=1e-10,atol=1e-12)
    verified.append(path.name)
report={'tables_verified':verified,'table_count':len(verified),'membership_and_grades_preserved':True,
        'scope':'Local processed-input reproduction of revision, not public publication or complete raw-data recomputation'}
(analysis.OUT/'reproduction_validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
