# Shared execution helpers for the sensor and Figure 1 notebooks.
# All entry points are run from the project root.
figure1_python <- function() {
  candidates <- unique(c(Sys.getenv("FIGURE1_PYTHON"),
    path.expand("~/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3"),
    unname(Sys.which("python3"))))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  probe <- shQuote("import openpyxl, reportlab, pypdf")
  for (python in candidates) {
    result <- suppressWarnings(system2(python, c("-c", probe), stdout = FALSE, stderr = FALSE))
    if (identical(result, 0L)) return(python)
  }
  stop("Figure 1 needs Python with openpyxl, reportlab and pypdf. Set FIGURE1_PYTHON to that Python executable.")
}

refresh_apex_records <- function(force = FALSE) {
  inputs <- c(list.files("Data/APEX", pattern = "^Tank.*[.]csv$", full.names = TRUE),
    list.files("Data/APEX/source_workbooks", pattern = "[.]xlsx$", full.names = TRUE),
    "scripts/import_apex_records.py")
  outputs <- file.path("Data/APEX/derived", c("agreeing_records.csv",
    "conflicting_timestamps.csv", "conflict_source_rows.csv", "gaps_over_24h.csv", "summary.csv"))
  rebuild <- force || !all(file.exists(outputs))
  if (!rebuild) rebuild <- max(file.mtime(inputs)) > min(file.mtime(outputs))
  if (rebuild) {
    result <- system2(figure1_python(), shQuote("scripts/import_apex_records.py"))
    if (result != 0L) stop("APEX source import failed; no figure was rebuilt.")
  }
  invisible(outputs)
}
