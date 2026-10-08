# Mitraseeds

Data and analysis code for:

> Portela LHX, Souza EB, Nuñez-Florentín M, Nepomuceno A, Santos FAR, Zappi DC.
> *Seed micromorphology diagnoses species but does not support infrageneric groups in South American* Mitracarpus *(Spermacoceae, Rubiaceae).*
> Submitted to **Plant Ecology and Evolution**.

Seeds of 27 South American species of *Mitracarpus* were scored for nine unordered
qualitative characters and analyzed phenetically: Gower dissimilarity, UPGMA
clustering, PCoA with vector fitting, PERMANOVA/PERMDISP, Mantel tests, and an
exhaustive search for the smallest subset of characters that still separates every
species.

Everything here is reproducible from a single script. Nothing in the paper was
computed by hand.

---

## Contents

| File | What it is |
|---|---|
| `Supplementary_File_2_pipeline.R` | The whole analysis, from the raw matrix to every figure and table. Self-contained: the data are embedded in it as plain text. |
| `Mitracarpus_seed_matrix.xlsx` | The same data in a readable workbook — the matrix, the 52 character states, seed dimensions, and biome assignments. Use this if you want the data without reading R code. |
| `FigS1A_effort_distance.pdf` | Mean Gower distance against the number of specimens examined. |
| `FigS1B_effort_variability.pdf` | Within-species variability against sampling effort. |
| `FigS2_silhouette.pdf` | Mean silhouette width across *k*, with the PAM cross-check. |
| `FigS3_MCA.pdf` | Multiple Correspondence Analysis of the complete cases. |
| `FigS4_accumulation.pdf` | Greedy character-accumulation curve for species resolution. |
| `FigS5_permdisp.pdf` | PERMDISP dispersion by grouping. |
| `FigS6_PCoA_biome.pdf` | The PCoA coloured by biome instead of by cluster. |

The four main-text figures (UPGMA phenogram, PCoA, ordination fit vs diagnostic
contribution, cluster × biome composition) are produced by the same script and are
not duplicated here.

---

## The data

**`Matrix`** — 27 species × 9 characters, as numeric state codes.

| Code | Character |
|---|---|
| `VGS` | Ventral groove shape |
| `SSh` | Seed shape |
| `AW` | Anticlinal wall configuration |
| `PW` | Periclinal wall topography |
| `EOrn` | Exotestal ornamentation |
| `CS` | Exotestal cell shape |
| `DCD` | Dorsal cruciform depression |
| `TVGDV` | Termination of the ventral groove extensions in dorsal view |
| `Colour` | Seed color |

Two things matter when reusing this matrix:

- **An empty cell (`NA` in the script) means the state could not be determined**
  from the material available. It is not a zero and not an absence. Only
  *M. carajasensis* has missing states, for four characters.
- **The characters are unordered.** State 3 is not "between" states 2 and 4; the
  numbers are labels, not quantities. Gower dissimilarity is computed on factors
  for exactly this reason. Treating them as numeric would produce a different and
  meaningless distance matrix.

**`Seed_dimensions`** — observed length and width ranges per species, in mm, with
means and the length/width ratio.

**`Biomes`** — one terrestrial biome per species, following Dinerstein et al. (2017),
referring to the **distribution of the species** rather than to the locality of any
one specimen. The `conf` column records how the assignment was made: `H` where the
distribution is stated in a revision or protologue, `L` where it was inferred from
the localities of the specimens examined. Three species carry `L` and the paper
flags them.

---

## Running the analysis

```r
install.packages(c("vegan", "cluster", "ape", "dendextend",
                   "FactoMineR", "ggplot2", "ggrepel", "scales"))

source("Supplementary_File_2_pipeline.R")
```

Tested on **R 4.6.1** with `vegan` 2.7-5 and `cluster` 2.1.8.2. The script calls
`sessionInfo()` at the end, so the exact environment of any run is recorded in its
own output.

It writes, into the working directory:

- 11 figures, each as a 600 dpi TIFF and a vector PDF
- 7 tables as CSV (`Table3_characters.csv`, `TableS2_gower_matrix.csv`,
  `TableS3_pcoa_eigenvalues.csv`, `TableS4_accumulation.csv`,
  `TableS5_revised_subtypes.csv`, `TableS6_biomes.csv`, `TableS7_species_scores.csv`)
- `results.txt`, the complete console output

Runtime is a few minutes, most of it in the permutation tests.

### Reproducibility

Every permutation test runs 9,999 permutations and **its own `set.seed()` call
immediately before it**, rather than one seed at the top of the script. This is
deliberate: with a single seed, each test consumes random numbers left by the
previous one, so adding or reordering a test shifts the results of all the tests
after it. Seeding each test individually makes every *p*-value independent of the
order in which the script is run.

One consequence worth knowing: *p* for seed shape in the vector fitting is 0.0486.
It is genuinely below 0.05, but barely. A different RNG or a change in `vegan`'s
permutation code could move it across the threshold.

---

## License

- **Code** (`Supplementary_File_2_pipeline.R`) — [MIT](LICENSE-MIT)
- **Data and figures** (`.xlsx`, `.pdf`) — [CC BY 4.0](LICENSE-CC-BY-4.0)

Reuse either freely, with attribution to the paper above.

---

## Contact

Luís Henrique Ximenes Portela — Escola Nacional de Botânica Tropical, Instituto de
Pesquisas Jardim Botânico do Rio de Janeiro (JBRJ) (ximenes849@gmail.com)

Please open an issue for anything that does not reproduce.
