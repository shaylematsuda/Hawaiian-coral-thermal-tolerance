# Uncommitted-change inventory — 2026-11-15

Last commit: **f4d01df — 2026-11-04 16:40:02 -0600 — Update supplemental physiology and community figures**.

Scope: every tracked modification/deletion and every non-ignored untracked file in this checkout. This includes work already present before the current Figure 1/Figure 3 edits.

**190 paths:** 2 deleted, 12 modified, 176 untracked/new.

**Nothing is currently staged.** An ordinary commit now would contain no file changes. This inventory describes what would be included if all listed paths were staged. No staging or commit was performed.

## Main substantive changes

- Figure 1: satellite collection map; coded trait table and extracted photos; real APEX treatment means ± SE; documented logger-removal exclusions; elapsed-day x-axis with day 0 = March 20; automatic export/import and assembly.
- Figure 3: DIV composition above the smaller ITS2-profile row, shared ALG sample order/labels, profile legend, class/genus composition with original plotting rules, explicit unclassified bacteria, approved widths/outlines, assembled portrait figure.
- Other checkout work: 16S/ITS2 manuscript output mappings, supplementary ASV and ANCOM figures, S17 tables, physiology/community figure updates, and review documents.
- Panel C citation placeholders * and ** remain unchanged pending citation numbers.

## Files to review separately before staging

- `.DS_Store` is local Finder metadata; `ACE21.Rproj` only adds a project ID.
- Two Figure 4 heatmap-candidate files are currently deleted; these deletions predate this request.
- Review folders, earlier previews, draft assets, and manuscript/audit files are included in this complete inventory, but need not be included in a figure-focused commit.


## Complete path list

Status: `M` = tracked modification; `D` = tracked deletion; `??` = untracked/new.

### Local/project metadata (2)

| Status | Path | Size |
|---|---|---:|
| `M` | `.DS_Store` | 14.0 KB |
| `M` | `ACE21.Rproj` | 0.2 KB |

### Analysis and figure code (11)

| Status | Path | Size |
|---|---|---:|
| `M` | `ACE21_16s_analysis.Rmd` | 149.2 KB |
| `M` | `ACE21_Sensor data.Rmd` | 1.7 KB |
| `M` | `Data/16s_analysis/bac_figures/procedure/ambient_figures_ace21.R` | 30.0 KB |
| `M` | `ITS2_analysis_update.Rmd` | 120.0 KB |
| `??` | `scripts/assemble_figure3.R` | 5.9 KB |
| `??` | `scripts/build_figure1.py` | 7.7 KB |
| `??` | `scripts/import_apex_records.py` | 3.7 KB |
| `??` | `scripts/make_its2_figure3a.R` | 5.9 KB |
| `??` | `scripts/plot_T0_class_unclassified.R` | 1.8 KB |
| `??` | `scripts/plot_apex_temperature.R` | 5.5 KB |
| `??` | `scripts/render_figure3.R` | 1.3 KB |

### Other loose project files (5)

| Status | Path | Size |
|---|---|---:|
| `M` | `CombinedPhysioPlots_T0.pdf` | 12.7 KB |
| `??` | `HC_counts_for_FC_20181213.xlsx` | 41.4 KB |
| `??` | `Supplement_16S_ASV_composition.png` | 970.9 KB |
| `??` | `Supplementary_Figure_S5_Physiological_responses_T1_TF.pdf` | 20.7 KB |
| `??` | `University of Hawaii Mail - zoox FCM (dragged).pdf` | 168.0 KB |

### Figure outputs (including older previews and deleted files) (51)

| Status | Path | Size |
|---|---|---:|
| `??` | `Figures/16S_tree_by_clade.pdf` | 21.3 KB |
| `D` | `Figures/ACE21_Figure4_effect_size_heatmap_candidate.pdf` | deleted |
| `D` | `Figures/ACE21_Figure4_effect_size_heatmap_candidate.png` | deleted |
| `??` | `Figures/Fig1.pdf` | 1893.8 KB |
| `??` | `Figures/Fig1.png` | 680.4 KB |
| `??` | `Figures/Fig1A_collection_map.pdf` | 907.8 KB |
| `??` | `Figures/Fig1A_collection_map.png` | 2711.3 KB |
| `??` | `Figures/Fig1B_temperature_available.pdf` | 160.9 KB |
| `??` | `Figures/Fig1B_temperature_available.png` | 267.6 KB |
| `??` | `Figures/Fig1B_temperature_available_plot.rds` | 713.8 KB |
| `??` | `Figures/Fig1B_temperature_mean_SE.pdf` | 158.4 KB |
| `??` | `Figures/Fig1B_temperature_mean_SE.png` | 497.4 KB |
| `??` | `Figures/Fig1B_temperature_mean_SE_plot.rds` | 932.6 KB |
| `??` | `Figures/Fig1C_traits_draft.pdf` | 461.5 KB |
| `??` | `Figures/Fig3.pdf` | 54.0 KB |
| `??` | `Figures/Fig3.png` | 364.2 KB |
| `??` | `Figures/Fig3A_ITS2.pdf` | 51.6 KB |
| `??` | `Figures/Fig3A_ITS2.png` | 178.1 KB |
| `??` | `Figures/Fig3_class_data.csv` | 24.3 KB |
| `??` | `Figures/Fig3_genus_data.csv` | 588.9 KB |
| `??` | `Figures/Fig3_plot.rds` | 3024.7 KB |
| `??` | `Figures/Figure3_T0_class_with_unclassified.pdf` | 6.4 KB |
| `??` | `Figures/Figure3_T0_class_with_unclassified.png` | 115.0 KB |
| `??` | `Figures/Figure3_T0_class_with_unclassified_data.csv` | 24.3 KB |
| `M` | `Figures/Figure_S7_ANCOM_treatment_by_species.pdf` | 87.1 KB |
| `M` | `Figures/Figure_S7_ANCOM_treatment_by_species.png` | 387.7 KB |
| `??` | `Figures/ITS2_Fig3A_plot.rds` | 1727.8 KB |
| `??` | `Figures/ITS2_Fig3A_sample_order.csv` | 11.2 KB |
| `??` | `Figures/ITS2_T0_samp_leg.pdf` | 51.6 KB |
| `??` | `Figures/ITS2_T0_samp_leg.png` | 178.1 KB |
| `M` | `Figures/ITS2_lineages_across_time_FigS15.pdf` | 17.1 KB |
| `M` | `Figures/ITS2_lineages_across_time_FigS15.png` | 566.1 KB |
| `??` | `Figures/New SFig20.png` | 566.1 KB |
| `M` | `Figures/PCA_T0_physio.pdf` | 9.8 KB |
| `??` | `Figures/Supplementary_Figure_S10_pvar_ASV_taxa.pdf` | 31.9 KB |
| `??` | `Figures/Supplementary_Figure_S10_pvar_ASV_taxa.png` | 425.4 KB |
| `??` | `Figures/Supplementary_Figure_S11_pacu_ASV_taxa.pdf` | 29.6 KB |
| `??` | `Figures/Supplementary_Figure_S11_pacu_ASV_taxa.png` | 414.0 KB |
| `??` | `Figures/Supplementary_Figure_S13_16S_T0_by_ALG.pdf` | 6.7 KB |
| `??` | `Figures/Supplementary_Figure_S13_16S_T0_by_ALG.png` | 116.6 KB |
| `??` | `Figures/Supplementary_Figure_S20_ITS2_ALG_stability.pdf` | 17.1 KB |
| `??` | `Figures/Supplementary_Figure_S20_ITS2_ALG_stability.png` | 566.1 KB |
| `??` | `Figures/Supplementary_Figure_S21B_16S_Endozoicomonadaceae_by_ALG.pdf` | 21.3 KB |
| `??` | `Figures/Supplementary_Figure_S21B_16S_Endozoicomonadaceae_by_ALG.png` | 529.9 KB |
| `??` | `Figures/Supplementary_Figure_S8_mcap_ASV_taxa.pdf` | 36.8 KB |
| `??` | `Figures/Supplementary_Figure_S8_mcap_ASV_taxa.png` | 484.7 KB |
| `??` | `Figures/Supplementary_Figure_S9_pcomp_ASV_taxa.pdf` | 50.4 KB |
| `??` | `Figures/Supplementary_Figure_S9_pcomp_ASV_taxa.png` | 620.5 KB |
| `??` | `Figures/_png_test.png` | 12.7 KB |
| `??` | `Figures/rerun_previews/Supplementary_Figure_S8_Mcap_taxa_rerun.pdf` | 37.5 KB |
| `??` | `Figures/rerun_previews/Supplementary_Figure_S8_Mcap_taxa_rerun.png` | 506.8 KB |

### APEX original inputs and archived notebooks (9)

| Status | Path | Size |
|---|---|---:|
| `??` | `Data/APEX/ACE21_Sensor_data_legacy.Rmd` | 5.9 KB |
| `??` | `Data/APEX/Apex_Temp.Rmd` | 1.2 KB |
| `??` | `Data/APEX/Tank1234.csv` | 172.2 KB |
| `??` | `Data/APEX/Tank5678.csv` | 169.0 KB |
| `??` | `Data/APEX/Tank911.csv` | 106.7 KB |
| `??` | `Data/APEX/source_workbooks/14March2017 Tank9-12.xlsx` | 420.3 KB |
| `??` | `Data/APEX/source_workbooks/14march2017 5-8.xlsx` | 473.4 KB |
| `??` | `Data/APEX/source_workbooks/APEX Temp Tanks 1-4.xlsx` | 536.2 KB |
| `??` | `Data/APEX/source_workbooks/Tanks 1-4 June2017.xlsx` | 157.5 KB |

### Derived data and APEX documentation (14)

| Status | Path | Size |
|---|---|---:|
| `??` | `Data/APEX/README.md` | 1.9 KB |
| `??` | `Data/APEX/derived/APEX_sampling_window_treatment_mean_SE.csv` | 817.4 KB |
| `??` | `Data/APEX/derived/agreeing_records.csv` | 3264.2 KB |
| `??` | `Data/APEX/derived/conflict_source_rows.csv` | 457.4 KB |
| `??` | `Data/APEX/derived/conflicting_timestamps.csv` | 186.1 KB |
| `??` | `Data/APEX/derived/gaps_over_24h.csv` | 0.6 KB |
| `??` | `Data/APEX/derived/logger_removal_May1_15_tanks1_4.csv` | 93.5 KB |
| `??` | `Data/APEX/derived/sampling_window_cleaning_exclusions.csv` | 0.2 KB |
| `??` | `Data/APEX/derived/summary.csv` | 0.2 KB |
| `??` | `Data/APEX/temperature_available_long.csv` | 3095.4 KB |
| `??` | `Data/APEX/temperature_coverage.csv` | 0.8 KB |
| `??` | `Data/derived/Table_S17a_16S_Treatment_x_ALG_PERMANOVA.csv` | 2.5 KB |
| `??` | `Data/derived/Table_S17b_16S_Treatment_x_ALG_dispersion.csv` | 1.0 KB |
| `??` | `Data/derived/Table_S17c_16S_Treatment_x_ALG_pairwise_dispersion.csv` | 1.1 KB |

### Figure source assets (10)

| Status | Path | Size |
|---|---|---:|
| `??` | `Figures/assets/Fig1/Fig1B_temperature.pdf` | 151.9 KB |
| `??` | `Figures/assets/Fig1/Fig1B_temperature_plot.rds` | 952.4 KB |
| `??` | `Figures/assets/Fig1/coral_traits_source.csv` | 0.6 KB |
| `??` | `Figures/assets/Fig1/kaneohe_satellite.jpg` | 724.6 KB |
| `??` | `Figures/assets/Fig1/kaneohe_satellite_metadata.json` | 0.8 KB |
| `??` | `Figures/assets/Fig1/reference-000.jpg` | 417.5 KB |
| `??` | `Figures/assets/Fig1/reference-001.jpg` | 111.9 KB |
| `??` | `Figures/assets/Fig1/reference-002.jpg` | 82.9 KB |
| `??` | `Figures/assets/Fig1/reference-003.jpg` | 55.7 KB |
| `??` | `Figures/assets/Fig1/reference-004.jpg` | 44.2 KB |

### Review and diagnostic artifacts (79)

| Status | Path | Size |
|---|---|---:|
| `??` | `output/reviews/ancom_current_2026-10-29/ANCOM_current_T1_TF_all_species.pdf` | 42.1 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/ANCOM_current_T1_TF_all_species.png` | 196.9 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/ancom_current_summary.csv` | 0.4 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/mcap_T1_ancom.pdf` | 12.3 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/mcap_T1_ancom_res.csv` | 30.3 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/mcap_TF_ancom.pdf` | 6.6 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/mcap_TF_ancom_res.csv` | 11.3 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/pacu_T1_ancom.pdf` | 10.6 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/pacu_T1_ancom_res.csv` | 25.2 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/pacu_TF_ancom.pdf` | 10.7 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/pacu_TF_ancom_res.csv` | 28.1 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/pcomp_T1_ancom.pdf` | 7.4 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/pcomp_T1_ancom_res.csv` | 4.6 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/pcomp_TF_ancom.pdf` | 5.7 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/pcomp_TF_ancom_res.csv` | 2.9 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/pvar_T1_ancom.pdf` | 8.4 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/pvar_T1_ancom_res.csv` | 15.4 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/pvar_TF_ancom.pdf` | 8.7 KB |
| `??` | `output/reviews/ancom_current_2026-10-29/pvar_TF_ancom_res.csv` | 11.3 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/ANCOM_current_T1_TF_all_species.pdf` | 86.8 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/ANCOM_current_T1_TF_all_species.png` | 190.6 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/ancom_current_summary.csv` | 0.3 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/mcap_T1_ancom.pdf` | 21.9 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/mcap_T1_ancom_res.csv` | 31.1 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/mcap_TF_ancom.pdf` | 11.1 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/mcap_TF_ancom_res.csv` | 11.6 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/pacu_T1_ancom.pdf` | 19.8 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/pacu_T1_ancom_res.csv` | 25.8 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/pacu_TF_ancom.pdf` | 21.3 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/pacu_TF_ancom_res.csv` | 28.8 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/pcomp_T1_ancom.pdf` | 7.6 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/pcomp_T1_ancom_res.csv` | 4.7 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/pcomp_TF_ancom.pdf` | 6.5 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/pcomp_TF_ancom_res.csv` | 3.0 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/pvar_T1_ancom.pdf` | 14.5 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/pvar_T1_ancom_res.csv` | 15.8 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/pvar_TF_ancom.pdf` | 11.8 KB |
| `??` | `output/reviews/ancom_current_historical_style_2026-10-29/pvar_TF_ancom_res.csv` | 11.5 KB |
| `??` | `output/reviews/ancom_styled_corrected_2026-10-29/ANCOM_styled_all_species.pdf` | 88.4 KB |
| `??` | `output/reviews/ancom_styled_corrected_2026-10-29/ANCOM_styled_all_species.png` | 289.6 KB |
| `??` | `output/reviews/ancom_styled_corrected_2026-10-29/mcap_T1.pdf` | 21.9 KB |
| `??` | `output/reviews/ancom_styled_corrected_2026-10-29/mcap_TF.pdf` | 11.5 KB |
| `??` | `output/reviews/ancom_styled_corrected_2026-10-29/pacu_T1.pdf` | 19.6 KB |
| `??` | `output/reviews/ancom_styled_corrected_2026-10-29/pacu_TF.pdf` | 21.1 KB |
| `??` | `output/reviews/ancom_styled_corrected_2026-10-29/pcomp_T1.pdf` | 8.3 KB |
| `??` | `output/reviews/ancom_styled_corrected_2026-10-29/pcomp_TF.pdf` | 6.5 KB |
| `??` | `output/reviews/ancom_styled_corrected_2026-10-29/pvar_T1.pdf` | 14.4 KB |
| `??` | `output/reviews/ancom_styled_corrected_2026-10-29/pvar_TF.pdf` | 11.7 KB |
| `??` | `output/reviews/ancom_styled_corrected_2026-10-29/summary.csv` | 0.2 KB |
| `??` | `output/reviews/apex_discrepancies/APEX_discrepancies.pdf` | 1049.0 KB |
| `??` | `output/reviews/apex_discrepancies/APEX_discrepancies_May_detail.png` | 329.9 KB |
| `??` | `output/reviews/apex_discrepancies/APEX_discrepancies_overview.png` | 433.0 KB |
| `??` | `output/reviews/apex_discrepancies/APEX_other_tanks_May1_15.png` | 301.1 KB |
| `??` | `output/reviews/apex_discrepancies/APEX_sampling_window.pdf` | 161.7 KB |
| `??` | `output/reviews/apex_discrepancies/APEX_sampling_window.png` | 501.5 KB |
| `??` | `output/reviews/apex_discrepancies/APEX_sampling_window_mean_SE.pdf` | 158.4 KB |
| `??` | `output/reviews/apex_discrepancies/APEX_sampling_window_mean_SE.png` | 350.3 KB |
| `??` | `output/reviews/apex_discrepancies/APEX_sampling_window_treatment_mean_SE.csv` | 663.1 KB |
| `??` | `output/reviews/apex_discrepancies/APEX_two_versions.pdf` | 38.2 KB |
| `??` | `output/reviews/apex_discrepancies/APEX_two_versions.png` | 304.8 KB |
| `??` | `output/reviews/apex_discrepancies/agreeing_records.csv` | 3264.2 KB |
| `??` | `output/reviews/apex_discrepancies/conflict_source_rows.csv` | 457.4 KB |
| `??` | `output/reviews/apex_discrepancies/conflicting_timestamps.csv` | 186.1 KB |
| `??` | `output/reviews/apex_discrepancies/gaps_over_24h.csv` | 0.6 KB |
| `??` | `output/reviews/apex_discrepancies/sampling_window_cleaning_exclusions.csv` | 0.2 KB |
| `??` | `output/reviews/apex_discrepancies/summary.csv` | 0.2 KB |
| `??` | `output/reviews/commit_inventory_2026-11-15.md` | this report |
| `??` | `output/reviews/flow_cytometry_audit/flow_fragids_absent_from_metadata.csv` | 2.0 KB |
| `??` | `output/reviews/flow_cytometry_audit/flow_plate_map_audit.csv` | 12.1 KB |
| `??` | `output/reviews/flow_cytometry_audit/flow_wells_missing_symbiodinium_counts.csv` | 1.1 KB |
| `??` | `output/reviews/pvar_unclassified_bacteria/Figure3_unclassified_bacteria_note.txt` | 0.6 KB |
| `??` | `output/reviews/pvar_unclassified_bacteria/JQ236445_reference.gb` | 2.1 KB |
| `??` | `output/reviews/pvar_unclassified_bacteria/S1982_candidate1_blast.xml` | 43.5 KB |
| `??` | `output/reviews/pvar_unclassified_bacteria/S1982_candidate1_blast_hits.csv` | 3.0 KB |
| `??` | `output/reviews/pvar_unclassified_bacteria/S1982_candidate1_blast_submission.txt` | 41.9 KB |
| `??` | `output/reviews/pvar_unclassified_bacteria/S1982_candidate_sequences.fasta` | 5.3 KB |
| `??` | `output/reviews/pvar_unclassified_bacteria/S1982_candidate_summary.csv` | 0.4 KB |
| `??` | `output/reviews/pvar_unclassified_bacteria/T0_OTU_relative_abundances.csv` | 0.8 KB |
| `??` | `output/reviews/statistical_audit_packet/ACE21_statistical_audit_packet.docx` | 38.6 KB |

### Manuscript and statistical audit documents (9)

| Status | Path | Size |
|---|---|---:|
| `??` | `statistic audit packet/ACE21_Manuscript_20250722.docx` | 1120.2 KB |
| `??` | `statistic audit packet/ACE21_statistical_audit_packet.docx` | 38.6 KB |
| `??` | `statistic audit packet/Codex Image Oct 29, 2026, 10_42_43 AM.png` | 190.6 KB |
| `??` | `statistic audit packet/Supplementary_Tables-ACE21_CORRECTED_2026-10-27.xlsx` | 199.7 KB |
| `??` | `statistic audit packet/Supplementary_Tables-ACE21_CORRECTED_2026-10-27_PACU_REMOVED.xlsx` | 215.2 KB |
| `??` | `statistic audit packet/Supplementary_Tables-ACE21_CORRECTED_2026-10-27_REBUILT.xlsx` | 219.6 KB |
| `??` | `statistic audit packet/T0_top10_taxa_by_host_mean_percent.xlsx` | 13.5 KB |
| `??` | `statistic audit packet/output/ACE21_pvalue_audit_64_expressions.docx` | 40.6 KB |
| `??` | `statistic audit packet/output/ACE21_pvalues_reviewed_highlighted.docx` | 1120.0 KB |
