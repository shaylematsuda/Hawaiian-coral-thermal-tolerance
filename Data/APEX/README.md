# APEX temperature exports

Source: user-provided Downloads/APEX folder, copied without modification on 2026-11-13.

- Tank1234.csv: March 7-May 1, 2017; tanks 1-4.
- Tank5678.csv: March 8-May 1, 2017; tanks 5-8.
- Tank911.csv: March 7-April 24, 2017; tanks 9 and 11.
- Apex_Temp.Rmd: original 2017 analysis, retained as provenance (not executed).

These are NOT the later `Tanks*_full.csv` exports referenced by ACE21_Sensor data.Rmd.
All dates in these exports are month/day/year. The older script's day/month parser
for tanks 1-4 would misinterpret these files.

Reproducible available-records plot: `Rscript scripts/plot_apex_temperature.R`.
No values are removed, interpolated, or extrapolated. Lines break across gaps
longer than one hour. Coverage and long-format tables are generated here.

Tank allocation: Ambient 1, 4, 5, 8, 9; Elevated 2, 3, 6, 7, 11. This follows the
experimental allocation supported by the metadata; the isolated conflicting tank-9
metadata record does not redefine the tank's treatment.

Later records are needed to establish complete coverage through the final sampling.

## Approved workflow (2026-11-14)

The supplied XLSX workbooks are preserved in `source_workbooks/`.
Run `scripts/import_apex_records.py` with Python (openpyxl) to rebuild `derived/`.
Run `Rscript scripts/plot_apex_temperature.R` (or the sensor notebook) to print
and save the treatment mean ± SE plot, with the investigator-confirmed
logger-removal exclusions recorded explicitly in code and audit CSVs.
Available tanks are the replicates at each timestamp; no interpolation or
averaging across time is performed. Conflicting copies are not resolved by averaging.

`Figures/assets/Fig1/Fig1B_temperature.pdf` is the vector export imported by
`scripts/build_figure1.py`. Run that script with Python (reportlab, pypdf) to save
`Figures/Fig1.pdf`. Raw CSV and XLSX measurements are never modified.
