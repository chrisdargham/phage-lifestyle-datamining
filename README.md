# Predicting bacteriophage lifestyle from genome data

A data mining project on 4,320 phage genomes from the [Actinobacteriophage Database (PhagesDB)](https://phagesdb.org/phages/). The question is whether a phage's lifestyle (**temperate** or **virulent**) can be predicted from simple genome annotations, and how well such a model would work on new phages.

**Read the report:** [`phage_datamining.md`](phage_datamining.md) (opens on GitHub with all code, tables and figures). An HTML version is in [`phage_datamining.html`](phage_datamining.html).

## Main finding

On a random train/test split, a random forest reaches **99.3% accuracy**. But every PhagesDB cluster of related phages has a single lifestyle, so a random split puts close relatives in both the training and test sets. A baseline that only looks up the cluster name, with no genome data, scores **99.8%**.

When whole clusters are held out, the picture changes:

| Validation | Random forest |
|---|---|
| Random 70/30 split | 99.3% accuracy |
| Held-out clusters, 10 repeated 70/30 splits | 84.9% balanced accuracy (SD 11.5%) |
| Grouped 10-fold cross-validation | 93.7% accuracy, AUC 0.974 |
| Same, each cluster counted once | 80.4% accuracy |
| Same, clusters of 1–3 phages only | 58.8% accuracy (19 clusters) |
| Permutation test (100 shuffles of cluster labels) | p = 0.0099 |

The model is clearly better than chance on new clusters, but much weaker than the random split suggests, especially for small clusters. This matters because the unlabelled phages have no labelled relatives, and most of them are in clusters of 1 to 3 phages.

![Random split vs held-out clusters](figures/cluster-holdout-1.png)

## Workflow

| Step | Methods | Report section |
|---|---|---|
| 1. Data preparation | Filtering, feature construction, missing-value handling, cohort table | 1 |
| 2. Exploratory analysis | Summaries, distributions, group comparisons | 2 |
| 3. Regression | Simple and multiple linear models, interaction, polynomial term, F-tests, AIC/BIC, residual diagnostics, cluster-level t-test | 3 |
| 4. Classical classification | Logistic regression, LDA, QDA, ROC/AUC, 10-fold CV, bootstrap, best subset (BIC), stepwise (AIC), Youden threshold | 4 |
| 5. Trees and ensembles | Decision tree with cost-complexity pruning (1-SE rule), bagging, random forest, boosting, OOB error, variable importance | 5 |
| 6. Validation check | No-genomics baseline, random vs held-out-cluster splits, grouped 10-fold CV, per-cluster accuracy, feature ablation, permutation test | 6 |
| 7. Unsupervised learning | PCA, k-means (elbow, silhouette), hierarchical clustering (4 linkages) | 7 |
| 8. Prediction | Lifestyle predictions for the 237 unlabelled phages | 8 |

## Other results

- **Regression:** genome size depends mainly on morphotype and host (adjusted R² = 0.545). Virulent phages tend to have larger genomes, but the difference is not significant when each cluster counts once (p = 0.11 instead of p < 2e-16).
- **Classical classifiers:** about 85% accuracy on a random split, 78% under grouped cross-validation. High coding density and GC content go with a temperate lifestyle.
- **Feature ablation:** on new clusters, genome size/ORFs, GC content and DNA end type carry the most information that other features cannot replace. Host genus and morphotype add little.
- **Unsupervised learning:** PCA and clustering group phages by genome size, RNA genes and morphotype, not by lifestyle.

## Repository structure

```
phage_datamining.Rmd     the whole analysis (one file)
phage_datamining.md      report rendered for GitHub
phage_datamining.html    report rendered as HTML
render.R                 renders both report formats
install_packages.R       installs the required R packages
data/phages_db.xlsx      PhagesDB export (input)
figures/                 every figure in the report, as PNG
output/                  result tables (CSV) written by the analysis
```

Main result files in `output/`:

| File | Content |
|---|---|
| `model_comparison_random_split.csv` | All 8 classifiers on the random split |
| `holdout_repeats.csv`, `holdout_summary.csv` | Random vs held-out-cluster splits |
| `grouped_cv_predictions.csv` | Out-of-fold prediction for every labelled phage |
| `grouped_cv_per_cluster.csv` | Accuracy for each cluster |
| `feature_ablation.csv` | AUC lost when each feature group is removed |
| `permutation_null.csv` | Balanced accuracy for each label shuffle |
| `regression_models.csv` | Adjusted R², AIC and BIC of the regression models |
| `unlabelled_predictions.csv` | Predicted lifestyle and probability for the 237 unlabelled phages |

## How to reproduce

Requirements: R (4.1 or later) and RStudio.

1. Open `phage-lifestyle-datamining.Rproj` in RStudio.
2. Run `install_packages.R` once.
3. Run `render.R`, or open `phage_datamining.Rmd` and click **Knit**.

From a terminal:

```bash
Rscript install_packages.R
Rscript render.R
```

One render takes about 4 minutes (`render.R` builds both formats, so about 8). All paths are relative and all random steps use `set.seed(123)`. The number of permutations can be lowered with `N_PERM` at the top of the Rmd for a faster run. Small differences can appear with other package versions; the numbers in the report text are computed directly from the results.

## Limitations

Lifestyle labels in PhagesDB may have been assigned per cluster rather than tested for every phage. The data are dominated by siphoviruses and *Mycobacterium* hosts. The features are genome summaries, not the genes that directly control lysogeny (such as integrase or repressor genes).

## Author

Chris Dargham, B.Sc. Bioinformatics, Lebanese American University
