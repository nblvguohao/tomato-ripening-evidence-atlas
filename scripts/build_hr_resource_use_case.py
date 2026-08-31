"""A post-review worked query, not a prospective validation or grade revision."""
from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/hr_revision_20260831"
EXTERNAL=["GSE267238","GSE128739","GSE235023","GSE78733"]

def truth(value):
    return str(value).lower()=="true"

def main():
    candidates=pd.read_csv(OUT/"rin_candidate_evidence_revised.csv")
    contexts=pd.read_csv(OUT/"candidate_external_contexts.csv")
    records=[]
    for _,r in candidates.iterrows():
        e=contexts.loc[contexts.gene_id.eq(r.gene_id)].set_index("context")
        query_entry=r.annotation_theme=="Flavour and specialised metabolism"
        external_measured=sum(truth(e.loc[n,"measured"]) for n in EXTERNAL)
        external_same=sum(truth(e.loc[n,"same_direction"]) for n in EXTERNAL)
        protein_record=truth(r.protein_differential) and truth(r.protein_direction_concordant)
        retained=query_entry and protein_record and external_measured==4 and external_same==4
        if not query_entry:
            disposition="outside_frozen_annotation_theme; not a biological rejection"
        elif not protein_record:
            disposition="protein_measurement_follow_up; not evidence against function"
        elif external_same<4:
            disposition="stage_context_follow_up; point-sign discrepancy is not significant reversal"
        else:
            disposition="retained_for_expression_and_protein_remeasurement; not an efficacy ranking"
        records.append({
            "gene_id":r.gene_id,"display_name":r.display_name,"frozen_grade":r.candidate_grade,
            "biological_anchor":r.biological_anchor,"annotation_theme":r.annotation_theme,
            "query_entry":query_entry,"external_measured_count":external_measured,
            "external_same_direction_count":external_same,
            "concordant_differential_protein_record":protein_record,
            "retained_for_worked_query":retained,"disposition":disposition,
            "GSE78733_effect":e.loc["GSE78733","effect"],
            "boundary_No19_measured":e.loc["GSE285925_No19","measured"],
            "boundary_No19_same_direction":e.loc["GSE285925_No19","same_direction"],
            "boundary_No20_measured":e.loc["GSE285925_No20","measured"],
            "boundary_No20_same_direction":e.loc["GSE285925_No20","same_direction"],
            "development_reuse_same_direction":e.loc["GSE183836_healthy","same_direction"],
            "tissue_measured_n":r.tissue_context_count,
            "tissue_same_direction_n":r.tissue_same_direction_context_count,
            "query_status":"post_review_illustration; not preregistered or a validated selection rule"
        })
    table=pd.DataFrame(records)
    table.to_csv(OUT/"horticultural_use_case.csv",index=False)
    report={"question":"Which annotation-linked flavour/metabolism candidates have a concordant differential-protein record and matching point directions in all four external studies, for targeted expression/protein remeasurement?",
            "entry_rule":"Exact frozen annotation_theme == Flavour and specialised metabolism; all 35 candidates retained in the audit table",
            "retention_rule":"4/4 external studies measured and same direction, plus concordant source differential-protein record",
            "status":"Post-review worked example; not prospective validation, novelty evidence, grade reassignment, or experimental success prediction",
            "entry_n":int(table.query_entry.sum()),"retained_n":int(table.retained_for_worked_query.sum()),
            "retained_ids":table.loc[table.retained_for_worked_query,"gene_id"].tolist()}
    (OUT/"horticultural_use_case.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))

if __name__=="__main__": main()
