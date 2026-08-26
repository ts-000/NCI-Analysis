# Pipeline 3 — NCI Analysis and Figures

## Overview

Pipeline 3 takes the model and biomotif mainframes from Pipeline 2 and generates the final NCI analysis, figures, tables, and Excel outputs used in this project.

The analysis is handled by:

```text
nci_analysis.py
```

## Input

The script expects these two files in the working directory:

```text
mainframe.csv
mainframe_biomotifs.csv
```

Empty examples are included under `example/` only to show the required filenames.

## Requirements

```bash
pip install numpy pandas matplotlib XlsxWriter
```

## Run

```bash
python nci_analysis.py
```

The generated report is written to:

```text
mega_nci_report/
```

This includes the analysis figures, presentation and paper figures, supplementary figures, tables, raw figure data, and Excel workbooks.

## Notes

- The script uses the fixed substituent/conformer series and project-specific core-interaction definitions used in this work.
- Biomotif data are used for the corresponding geometry comparison plots.
- Running the script replaces an existing `mega_nci_report/` directory, so previous results should be copied elsewhere if they need to be kept.
