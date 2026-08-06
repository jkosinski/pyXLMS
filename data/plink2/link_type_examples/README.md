# pLink link-type examples (cross- / loop- / mono-linked)

Small real-data slices (8 rows each) from **PRIDE PXD051075** (A549 fixed cells, DSS),
`A549_FixedCells_DSS_5FDR.filtered_*_spectra.csv`. They illustrate that pLink deposits
crosslink results **split across three files that share one identical schema** (the same
21 columns), differing only in how the `Peptide` / `Peptide_Type` / `Proteins` fields are
encoded:

| file | `Peptide_Type` | `Peptide` example | `Proteins` example |
|---|---|---|---|
| `*cross-linked_spectra.csv` | `Cross-Linked` | `KALAAAGYDVEK(1)-AVAASKER(6)` | `sp|A|.. (k1)-sp|B|.. (k2)/` |
| `*loop-linked_spectra.csv`  | `Loop-Linked`  | `DQKKTQEQLALEMAELTAR(3)(4)`   | `sp|A|.. (k)/` (one protein, two in-peptide sites) |
| `*mono-linked_spectra.csv`  | `Mono-Linked`  | `APKPDGPGGGPGGSHMGGNYGDDR(3)` | `sp|A|.. (451)/` (one protein, one site) |

`read_plink` currently reads only the **cross-linked** file: it splits `Peptide` and
`Proteins` on `-` (assuming two peptides / two proteins), so a mono-/loop-linked file
raises `IndexError` on `.split("-")[1]`.

Per pyXLMS issue #201, the natural representation is:
- **mono-link** → a crosslink with the `beta_*` fields set to `None` (one linked site);
- **loop-link** → a regular crosslink with both sites on the same peptide/protein.

These files are provided as ready-made fixtures for that discussion / any future support.
