# Intermediate analysis workflow

Scripts in `03-scripts/R/` that produce **intermediate artefacts** (Seurat
objects, DEG tables, model fits) rather than manuscript figures. Every figure
script in `03-scripts/figures/` assumes these have already been run.

Run them in the numbered order below. All paths are relative to the repository
root; run everything from the root (the `.Rproj` sets this for you, and
`_quarto.yml` sets `execute-dir: project` so notebooks match).

`seq_functions.R` is a utility library, not a stage — every script `source()`s it.

## Run order

| # | Script | Writes | Requires |
|---|--------|--------|----------|
| 00 | `00_make_seurat_object.qmd` | `04-analysis/Seurats/Dec2024/` (full object) | raw/Trailmaker matrices, MapMyCells reassignments |
| 01 | `01_make_seurat_gigaclass.R` | `04-analysis/Seurats/Dec2024/seurat_{excitatory,inhibitory,glia}.Rds` | 00 |
| 02 | `02_dge_pseudobulk_deseq2.R` | `04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/` | 00 |
| 03 | `03_classify_ERG-LRG_deseq2.R` | ERG/LRG calls in the DESeq2 DEG directory | 01, 02 |
| 04 | `04_cache_subclasses.R` | `04-analysis/dge_method_comparison/subclass_cache/` | 00 |
| 05 | `05_dge_pseudobulk_voom.R` | `04-analysis/dge_method_comparison/voom/`, `0_method_concordance.csv` | 02, 04 |
| 06 | `06_classify_ERG-LRG_voom.R` | ERG/LRG calls in `04-analysis/dge_method_comparison/voom/` | 01, 05 |
| 07 | `07_dge_glmm_animal_random_effect.R` | `04-analysis/dge_method_comparison/glmm/` | 04 |
| 08 | `08_sample_sizes_ess.R` | `04-analysis/dge_method_comparison/0_per_subclass_replicates_and_ESS.csv` | 04 |

Stages 04–08 are the alternative-DGE-framework analyses added during revision.
07 and 08 are independent of each other and of 05/06; 05 must precede 06.

**Cache location.** Stage 04 writes ~2.7 GB of per-subclass Seurat caches. Set
`ACTDEP_CACHE_ROOT` to an existing directory (e.g. an external disk) to keep
them off the boot volume; otherwise they land under `04-analysis/`. See
`SubclassCacheDir()` in `seq_functions.R`.

**macOS note.** Figure 2—figure supplements 4 and 5 use RedRibbon, which needs
the patch in `03-scripts/patches/` applied before RedRibbon is installed, or
every entry point segfaults. See `03-scripts/patches/README.md`.

## Dependency graph

```mermaid
flowchart TD
    subgraph data["02-data"]
        RAW["raw / Trailmaker<br/>count matrices"]
        PUB["published data<br/>(Yao2023, Zhang2023)"]
    end

    S00["00_make_seurat_object.qmd"]
    S01["01_make_seurat_gigaclass.R"]
    S02["02_dge_pseudobulk_deseq2.R"]
    S03["03_classify_ERG-LRG_deseq2.R"]
    S04["04_cache_subclasses.R"]
    S05["05_dge_pseudobulk_voom.R"]
    S06["06_classify_ERG-LRG_voom.R"]
    S07["07_dge_glmm_animal_random_effect.R"]
    S08["08_sample_sizes_ess.R"]

    NUCLEI[("Seurats/Dec2024<br/>full object")]
    GIGA[("seurat_{excitatory,<br/>inhibitory,glia}.Rds")]
    DEG2[("DEGs/…_pseudobulk<br/>DESeq2 tables + ERG/LRG")]
    CACHE[("dge_method_comparison/<br/>subclass_cache")]
    VOOM[("dge_method_comparison/voom<br/>tables + ERG/LRG")]
    GLMM[("dge_method_comparison/glmm")]
    ESS[("0_per_subclass_<br/>replicates_and_ESS.csv")]

    RAW --> S00
    PUB --> S00
    S00 --> NUCLEI

    NUCLEI --> S01 --> GIGA
    NUCLEI --> S02 --> DEG2
    GIGA --> S03
    DEG2 --> S03
    S03 --> DEG2

    NUCLEI --> S04 --> CACHE

    CACHE --> S05
    DEG2 --> S05
    S05 --> VOOM

    GIGA --> S06
    VOOM --> S06
    S06 --> VOOM

    CACHE --> S07 --> GLMM
    CACHE --> S08 --> ESS

    NUCLEI --> FIGS["03-scripts/figures/*.R"]
    GIGA --> FIGS
    DEG2 --> FIGS
    VOOM --> FIGS
    GLMM --> FIGS
    ESS --> FIGS
```

## Which figures consume which artefact

```mermaid
flowchart LR
    NUCLEI[("full Seurat object")]
    GIGA[("gigaclass Seurats")]
    DEG2[("DESeq2 DEGs<br/>+ ERG/LRG")]
    VOOM[("voom DEGs<br/>+ ERG/LRG")]
    GLMM[("NB-GLMM fits")]

    NUCLEI --> F1["Figure1B / 1C-D / 1E-H<br/>Figure1_SuppFig1, 3-4"]
    NUCLEI --> F2["Figure2C-E / 2F-G<br/>Figure2_SuppFig6"]
    NUCLEI --> F4["Figure4A-B / 4C / 4D-E / 4F<br/>Figure4_SuppFig1, 2"]
    NUCLEI --> F5["Figure5A / 5B-E<br/>Figure6A-J / Figure6_SuppFig1-2"]

    DEG2 --> F2AB["Figure2A-B<br/>Figure2_SuppFig7"]
    DEG2 --> F3["Figure3AEI / 3BFJ /<br/>3CGK / 3DHL"]
    GIGA --> F3

    VOOM --> F2S1["Figure2_SuppFig1"]
    VOOM --> F3S1["Figure3_SuppFig1AEI / BFJ / CGK"]
    GIGA --> F3S1

    DEG2 --> FOREST["Figure2_SuppFig2, 3<br/>(forest plots)"]
    VOOM --> FOREST
    GLMM --> FOREST

    DEG2 --> RRHO["Figure2_SuppFig4, 5<br/>(RRHO)"]
    VOOM --> RRHO
    GLMM --> RRHO
```
