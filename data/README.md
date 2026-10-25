# Data

`phages_db.xlsx`: tabular export of the Actinobacteriophage Database, PhagesDB (https://phagesdb.org/phages/), converted to Excel. One sheet (`Sheet1`), one row per phage, 5,580 rows and 41 columns.

Columns used in the analysis:

| Column | Used as |
|---|---|
| `Phage.Name` | phage identifier |
| `Temperate.` | response: `True` = temperate, `False` = virulent, `None` = unlabelled |
| `Cluster` | PhagesDB cluster of related phages (used only for validation, never as a predictor) |
| `Host` | isolation host; the first word gives `host_genus` |
| `Finished.Sequence.` | only finished genomes are kept |
| `Genome.Length.bp.` | genome length (converted to kb) |
| `GC.` | GC content (%) |
| `X..ORFs`, `X..tRNAs`, `X..tmRNAs` | gene counts (also used per kb) |
| `End.Type` | DNA end type (missing values coded as `Unknown`) |
| `Morphotype` | SIPHO, MYO, PODO, TECTI |

The lifestyle labels are the values recorded in the database at the time of export; they were not re-verified for this project.
