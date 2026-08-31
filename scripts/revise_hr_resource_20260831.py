"""Export revision-specific resource semantics without changing frozen grades."""
from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/hr_revision_20260831"

def main():
    atlas=pd.read_csv(ROOT/"results/ripening_atlas_resource_v3.csv")
    atlas["occupancy_and_perturbation_supported"]=atlas.multiomics_high_confidence
    atlas["four_study_consensus_role"]="two_selection_sources_and_two_external_transfer_studies"
    atlas["legacy_high_confidence_note"]="Compatibility alias only; occupancy plus perturbation, not calibrated confidence or mandatory protein support"
    atlas["atlas_claim"]=atlas.all_four_same_direction.map({True:"two_source_two_external_direction_consensus",False:"frozen_signature_with_partial_direction_consensus"})
    phospho_note="Positive differential-site records only. No record is not a tested negative; the full measurable/testable phosphosite universe was not reconstructed."
    atlas["phosphosite_support_status"]=atlas.phosphosite_differential_BR8_vs_MG.eq(True).map({True:"differential_site_record_present",False:"no_record_testability_unknown"})
    atlas["phosphosite_testability_note"]=phospho_note
    atlas.to_csv(OUT/"ripening_atlas_resource_revised.csv",index=False)
    candidates=pd.read_csv(ROOT/"results/hr_submission/rin_candidate_evidence_table.csv")
    candidates["display_order"]=candidates.submission_priority_rank
    candidates["display_order_definition"]="Grade A/B/C; anchor first; limited annotation first; differential protein first; phosphosite first; descending absolute rin-1 effect; gene ID. Display order only, not calibrated functional priority."
    candidates["four_study_consensus_sources"]="GSE42783;GSE108415;GSE267238;GSE128739"
    candidates["selection_source_count"]=2
    candidates["external_grade_source_count"]=2
    candidates["protein_measurement_status"]=candidates.protein_support_status
    candidates["protein_statistical_support_note"]=candidates.protein_support_status.map(
        lambda s:"Source differential-protein record present" if s=="differential_direction_concordant" else
        "No differential-protein record provided; not evidence of a completed nonsignificant test" if s=="measured_direction_concordant" else
        "Measured opposite point-estimate direction" if s=="measured_direction_discordant" else "Not measured")
    candidates["phosphosite_record_status"]=candidates.phosphosite_differential.eq(True).map({True:"differential_site_record_present",False:"no_record_testability_unknown"})
    candidates["phosphosite_testability_note"]=phospho_note
    candidates["broad_tissue_measurement_coverage"]=candidates.tissue_context_count.ge(9)
    candidates["tissue_same_direction_context_count"]=(candidates.tissue_context_count*candidates.tissue_direction_consistency).round().astype(int)
    for table in (candidates,):
        for col in table.select_dtypes(include="object"):
            mask=table.gene_id.eq("Solyc06g051800")
            table.loc[mask,col]=table.loc[mask,col].str.replace("Uluisik et al.","Minoia et al.",regex=False)
    candidates.to_csv(OUT/"rin_candidate_evidence_revised.csv",index=False)
    # The standalone old registry predates the completed submission screen.
    # Use the completed, frozen candidate-table records, as in baseline Table S7.
    literature=candidates.loc[candidates.candidate_grade.eq("A"),[
        "gene_id","display_name","candidate_grade","submission_priority_rank",
        "literature_screen_status","literature_screen_scope","literature_screen_source",
        "literature_screen_url","literature_screen_conclusion","literature_screen_date"
    ]].rename(columns={"literature_screen_date":"screen_date"}).copy()
    assert not literature.literature_screen_status.str.contains("pending",case=False).any()
    for col in literature.select_dtypes(include="object"):
        mask=literature.gene_id.eq("Solyc06g051800")
        literature.loc[mask,col]=literature.loc[mask,col].str.replace("Uluisik et al.","Minoia et al.",regex=False)
    literature.to_csv(OUT/"grade_a_literature_registry_revised.csv",index=False)
    dictionary=[
        ("display_order","Deterministic figure/table ordering, not a validated prioritization score"),
        ("submission_priority_rank","Deprecated compatibility alias of display_order; no probability interpretation"),
        ("occupancy_and_perturbation_supported","Replicated admitted RIN occupancy plus independent rin-1 differential expression"),
        ("multiomics_high_confidence","Deprecated alias of occupancy_and_perturbation_supported; does not require protein evidence"),
        ("four_cohort_direction_consensus","Two selection sources plus GSE267238 and GSE128739; not four external validations"),
        ("bootstrap_same_direction_frequency","Fraction of all ordered endpoint-column resamples with reference sign; not a posterior probability"),
        ("effect_bootstrap_p025/p975","Conditional empirical-bootstrap percentiles on the stated effect scale; small-sample and normalization-fixed"),
        ("rank01","Within-sample percentile rank, late-minus-early group means"),
        ("log2cpm","log2(CPM+1), late-minus-early group means, not DESeq2 log2 fold change"),
        ("protein_measurement_status","Distinguishes source differential record, measured same sign only, measured opposite sign and absent measurement"),
        ("phosphosite_support_status / phosphosite_record_status",phospho_note),
        ("phosphosite_differential_BR8_vs_MG / phosphosite_differential","Legacy Boolean indicating presence in a differential-site record table; False does not establish a measured nonsignificant site"),
        ("broad_tissue_measurement_coverage","True iff measured in at least 9 of the 11 tissue contexts; contains no direction-agreement condition and is not a grade input"),
        ("tissue_same_direction_context_count","Number of measured tissue contexts matching reference direction; denominator is tissue_context_count"),
        ("GSE235023 role mapping","Historical planned_role=sensitivity; independent external-transfer sensitivity; does not select the signature or enter frozen candidate grades"),
    ]
    pd.DataFrame(dictionary,columns=["field","definition"]).to_csv(OUT/"revised_data_dictionary.csv",index=False)
    original=pd.read_csv(ROOT/"results/hr_submission/rin_candidate_evidence_table.csv")
    assert candidates.gene_id.equals(original.gene_id)
    assert candidates.candidate_grade.equals(original.candidate_grade)
    assert candidates.submission_priority_rank.equals(original.submission_priority_rank)
    (OUT/"resource_revision_validation.json").write_text(json.dumps({"signature_rows":len(atlas),"candidate_rows":len(candidates),"grades":candidates.candidate_grade.value_counts().to_dict(),"membership_grades_and_display_order_unchanged":True,"public_archive_rewritten":False},indent=2)+"\n")

if __name__=="__main__": main()
