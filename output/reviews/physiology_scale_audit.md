# Physiology plotting-scale audit

Date: 2026-10-09
Source: ACE21_physiology.Rmd

## Descriptive figures

| Trait | Original-scale input | Audit result |
|---|---|---|
| Biomass | AFDWmg.cm2 | Standalone and combined figure inputs already original-scale. |
| Algal cell density | zoox.cm2 | Standalone figures already original-scale; combined percent-change figures now undo log10(x + 1) on individual values before summaries. |
| Calcification | time1.growth / time2.growth / percent.change | Figure inputs already use growth values before the square-root model transformation. |
| Chlorophyll a | ChlaSA | Figure inputs already original-scale. |
| Fv/Fm | Yield | Retained alongside Yield_cu; percent-difference figures can now select it. Combined dataset now averages original Yield. |
| Protein | mg.mm2 | Standalone figures already original-scale; combined percent-change figures now undo the square root on individual values before summaries. |
| Photosynthesis | Photo / GrossP | Figure inputs use original-scale rates after unit conversions. Net versus gross choices were not changed. |
| Respiration | Resp / RespPos | Figure inputs use original-scale rates; existing sign conventions were retained. |
| P:R | PR / PR_real | Figure inputs are ratios without statistical power/log transformations. Existing ratio definitions were retained. |

Original-scale here includes the existing normalization, unit conversion and quality-control rules, including the 0.1 replacement for low-signal PAM observations.

## Analysis inputs

PAM models continue to use Yield_cu = Yield^3. The combined PhysioData dataset now uses means of Yield for Fv/Fm, affecting the downstream PCA and other multivariate analyses that consume it. Cell-density and protein transformations in that analysis dataset remain unchanged. PCA centering/scaling is retained.

PhysioData_plot is a separate original-scale copy for the two combined percent-change figure workflows. Model diagnostic histograms and residual plots still display the quantities needed for model checks.

## Other repairs needed by retained percent-difference figures

Removed an assignment that replaced actual trait labels with a recycled eight-label vector. Corrected the four-panel T1 plot to use the existing, flipped Mcapt1 object and the TF plot to refer to mcap_TF_plot instead of the undefined TF_TF_plot.

## Scope

This audit checks plotting scales, not the validity of every model or figure. Previously identified Gower metadata inclusion, potential dispersion row-order mismatch, and stale references in the later legacy NMDS section remain separate review items. Existing model warnings also need scientific review. Existing exported figures have not been overwritten by validation.

## Validation

- Parsed and executed active notebook code in a clean R process through the retained T1/TF percent-difference figure section. Figure/file writes were suppressed during validation.
- Confirmed all 655 model-fitting and emmeans code lines match the pre-edit version.
- Matched the combined dataset's Fv/Fm values to group means calculated directly from original Yield.
- Matched original-scale cell-density and protein plotting values to their source records.
- Matched PAM percent differences to direct original-scale calculations.
- Verified all 72 trait/species/time-point rows retain their correct trait labels and built five retained PAM/TF plot objects.
- Existing singular-fit and plotting warnings remain; this is not a claim that the full notebook or statistical design is validated.
