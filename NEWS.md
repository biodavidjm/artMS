# artMS 1.33.3 (2026-11-25)

- Fix the vignette link to the changelog (it pointed to the old `NEWS` file)

# artMS 1.33.2 (2026-11-25)

- The website (<https://biodavidjm.github.io/artMS/>) is now built and published automatically by a GitHub Actions workflow on every push to `master`; the built files are no longer stored in the repository
- `NEWS` converted to `NEWS.md`, so the website has a Changelog page
- Website: new navigation bar (Get started, Reference, Changelog, Bioconductor) and functions grouped by task on the Reference page
- DESCRIPTION `URL` now lists the website and the GitHub repository (replaces http://artms.org)
- Fix broken links to MaxQuant, CORUM and the configuration section of the vignette

# artMS 1.33.1 (2026-11-16)

- Remove unused dependencies: `plotly` and `getopt`
- Replace Travis CI with a GitHub Actions workflow that checks the package on the Bioconductor devel Docker image
- Update documentation (roxygen2 8.1.1), README badges and the pkgdown website

# artMS 1.31.1 (2026-11-16)

- Fix the unit test that made R CMD check fail on Bioconductor: MSstats now returns an extra `LABEL` column in the `ModelQC` output of `artmsQuantification()`. The test now checks that the columns used by artMS are present instead of an exact number of columns.
- Enrichment analysis now uses the `gprofiler2` package instead of `gProfileR`. `gProfileR` is deprecated and its server no longer accepts requests, so the enrichment in `artmsEnrichProfiler()`, `artmsEnrichLog2fc()` and `artmsAnalysisQuantifications()` was failing. The results keep the same column names as before, with these differences:
    - the `domain` column now uses the g:Profiler source names (`GO:BP`, `KEGG`, `REAC`, ...)
    - the gProfileR-only columns `subgraph.number` and `relative.depth` are no longer included
    - `OMIM` is no longer available in g:Profiler and is ignored
    - install `gprofiler2` to use these functions: `install.packages("gprofiler2")`
- Fix typo in the enrichment error message

# artMS 1.10.2 (2021-07-13)

- Bug fix affecting `artmsAnalysisQuantifications()`
- Allow `artmsAnalysisQuantifications()` to process previous versions of artMS
- Change Extension of `artms_sessionInfo_quantification` file from `.txt` to `.log`
- New parameters available in `artmsQuantification()`, including:
    - Parameter `printTables`. Default `TRUE`, prints tables. FALSE otherwise.
    - Parameter `return_results_object`. If TRUE, it returns a list of data frames with MSstats results. Default is FALSE
    - If both `printTables` and `printPDF` are FALSE, then `return_results_object` becames TRUE and a list of data frames is returned
- Change default parameters of the configuration files: less verbose

# artMS 1.10.1 (2021-06-30)

- Addressing major changes in MSstats
    - R version larger than 4.1 is now required
    - Fractions: the option "Fractions" is removed from the configuration file. If fractions are present, the user must include a "Fraction" column in the keys.txt file, which artMS will detect automatically.
    - `artmsQualityControlEvidenceBasic`: fraction parameter no longer required (automatically detected from the keys file)
    - `keys.txt`: use `Fraction` instead of `FractionKey`
- External packages used exclusively by the `artmsAnalysisQuantifications` function are not required. Those packages will have to be installed before running this function.
- `artmsAvgIntensityRT`: argument `species` is not longer required
- Example datasets:
    - `artms_data_ph_evidence`: the size has been significantly reduced (bioconductor requirement). Only two biological replicates and 1/20 of the lines selected randomly. Including only 36 columns from the original evidence file
    - `artms_data_ph_msstats_results`: output from running `artmsQuantification` on the full version of the evidence file, including 4 biological replicates (instead of the reduced version available in the package)
    - `artms_data_ph_msstats_modelqc`: output from running `artmsQuantification` on the full version of the evidence file, including 4 biological replicates (instead of the reduced version available in the package)

# artMS 1.10.0 (2021-05-12)

- Bioconductor release 3.13
- Several adjustments in `artmsAnalysisQuantifications` function:
    - Use NAs instead of 0 when reshaping the data for missing values (`-log2fc-wide.txt` and `-log2fc-long.txt` files)
    - `plotClusteringAnalysis` only available when 3 or more comparisons available in the configuration file
- R > v4.0.0 is now required
- Update documentation & vignette
- Code cleaning

# artMS 1.8.3 (2021-04-05)

- `artmsProtein2SiteConversion` now supports uniprot id isoforms (thanks Emily King)
- Update vignette to make clearer how to provide several protein ids in "normalization_reference" (thanks Olga Schubert)

# artMS 1.8.2 (2021-03-18)

- `artmsProtein2SiteConversion`: New PTM available as argument, a new user defined `PTM:XXX:yy`. Check documentation to find out more
- QC plots: by default, all qc plots are now output to a folder directory, by type (qc-basic, qc-extended, qc_summary)
- artMS working directory: artMS will create all the folders and subfolders relative to the working directory. No need to specify the full path to the working directory, but the user must set the working directory: setwd("/path/to/working/directory/")
- Configuration file data object updates:
    - "output" the user can add a folder where would like to have the output results file. For example, "output: results_202003/example-results.txt" would create the "resutls_202003" folder (if it does not exist) with all the results files available there
    - "LFC" (log2fc) updated to range -0.58 to 0.58, i.e., a fold change larger than a 1.5 (instead of 2 as before)
- Update and improve documentation and vignette
- Several bug fixes

# artMS 1.8.1 (2020-10-27)

- Update "plotPCA" message
- Update documentation and vignettes
- Code cleaning

# artMS 1.8.0 (2020-10-27)

- New bioconductor version

# artMS 1.6.7 (2020-10-22)

- The following plots are now deprecated
    - .clustering.log2fcSign.all-zoom.pdf
    - .clustering.log2fc.all-zoom.pdf
- Fix bugs

# artMS 1.6.6 (2020-10-21)

- Convert new MaxQuant format of PTMs to the old format
- MSstats messages are not displayed by default when using artmsQuantification. The user can enable MSstats messages by selecting "display_msstats = TRUE"
- Prevent artmsWriteConfigYamlFile() from overwriting an existing configuration file unless the user allows it (overwrite = TRUE)
- printPDF now available in all functions printing plots to pdf, which means that notebooks can be used and print all plots. Default is still printPDF = TRUE
- Fix bugs

# artMS 1.6.5 (2020-05-20)

- Fix ggplot warnings (caused by NA values)
- Fix artmsAnalysisQuantification reproducibility plots
- Improves artmsQualityControlEvidenceBasic() correlation matrix clustered plot
- Fix pca01.pdf plot
- New pca04.pdf plot (dot plot)
- artmsAnalysisQuantifications check point: check if sufficient data is available

# artMS 1.6.4 (2020-05-12)

- Fix Quality Control functions to handle a small number of runs (less than 5)
- New argument "printPDF" for artmsQualityControlSummaryExtended, to select whether to print plots to PDFs (default = TRUE)
- Vignette: example plots added

# artMS 1.6.3 (2020-05-06)

- Bug Fixes affecting artmsAnalysisQuantifications()

# artMS 1.6.2 (2020-05-05)

- Fix NEWS formatting
- Update vignette with AC options

# artMS 1.6.1 (2020-04-29)

- Fix NEWS formatting
- Update vignette with AC options

# artMS 1.6.0 (2020-04-27)

- Major package updates. It addresses most of the warnings and messages from required packages
- Analysis of Protein Acetylation now supported.
- General to all Quality control plots: update default names for all pdf outputs
- QC artmsQualityControlEvidenceExtended():
    - New plots: peptide and protein overlap across bioreplicates and conditions (plotPEPTOVERLAP and plotPROTOVERLAP)
    - Updates and improves all the plots, affecting font size and other aesthetic options
    - Change parameter "plotIC" to "plotPCA" in artmsQualityControlEvidenceExtended()
- QC artmsQualityControlEvidenceBasic():
    - Depracated: Plot correlation distributions
    - Updated and improved all the plots, affecting font size and other aesthetic
- More data testing datasets available:
    - artms_data_ph_contrast: contrast object for available the PH example dataset
    - artms_data_ph_config: configuration yaml object for the available PH example dataset
    - artms_data_ph_msstats_results: results data object from running artmsQuantification
    - artms_data_ph_msstats_modelqc: modelqc data object from running artmsQuantification
- Vignette updated
- Other minor fixes / improvements

# artMS 1.4.3 (2020-02-08)

- artmsAnalysisQuantifications(), several updates:
    - change default value for "l2fc_thres" to 1
    - change defult value for "output_dir" to "analysis_quant"
    - NEW option: "outliers". It allows to remove outliers from abundance data
- Fix bug affecting artmsQualityControlSummaryExtended()

# artMS 1.4.2 (2019-12-06)

- Fixed a crash when impute fails due to missing outlier data (part 2)

# artMS 1.4.1 (2019-12-06)

- Fixed a crash when impute fails due to missing outlier data

# artMS 1.4.0 (2019-11-01)

- New bioconductor release

# artMS 1.2.9 (2019-10-28)

- Bug fix on artmsProtein2SiteConversion affecting UB conversion

# artMS 1.2.8 (2019-10-21)

- Enhancement: artmsProtein2SiteConversion can handle now full uniprot ids (as extracted by default in MaxQuant). For example, from sp|Q86U42|PABP2_HUMAN it will extract the Uniprot entry id Q86U42 as the main id and the "protein to site"" conversion is performed.

# artMS 1.2.7 (2019-10-14)

- Remove reshape2 dependency (depracated) and use data.table melt function instead
- Update documentation: Notes on imputation

# artMS 1.2.6 (2019-08-01)

- Optimize Analysis of Quantifications: adding check points when enrichment tool gProfileR is not available; clean up code; fix bugs

# artMS 1.2.5 (2019-07-25)

- Update vignette and readme

# artMS 1.2.4 (2019-07-25)

- artMS requires now R > 3.6
- Check point for empty evidence files

# artMS 1.2.3 (2019-06-11)

- Fix: Refactored code to better handle SILAC data in QC and quantifications
- New function to hack SILAC data
- Keeps adding `na.rm = TRUE` to ggplots
- Improves QC plots
- Overall code formating (going back to my prefered style)

# artMS 1.2.2 (2019-05-11)

- Fix: add PH check for new maxquant versions in basic QC plot function
- New function to add error bars in plots

# artMS 1.2.1 (2019-05-07)

- Fix artmsProtein2SiteConversion package dependency

# artMS 1.2.0 (2019-05-01)

- Bioconductor release: update version (from 1.1 to 1.3)

# artMS 1.0.17 (2019-05-27)

- Minor changes in code

# artMS 1.0.16 (2019-04-27)

- The configuration file also includes the option to activate quality control based on the summary.txt file
- Improve messages to the user

# artMS 1.0.15 (2019-04-24)

- artMS detects RefSeq IDs. It does not support yet enrichments and annotations
- Update artmsAnalysisQuantifications to adjust for previous changes in PTMs: it better handles 'protein groups'
- Doesn't remove outliers

# artMS 1.0.14 (2019-04-22)

- Update artmsProtein2SiteConversion: it runs as default the "Leading razor protein".
- Update vignette

# artMS 1.0.13 (2019-04-20)

- Important artmsProtein2SiteConversion update: the protein column id to be converted to the 'site/peptide specific' notation (i.e. ProteinID_AAnumber) can be now selected by the user. This is important because if the user decides to remove the protein groups, the 'leading razor protein' is chosen as the main protein id column, therefore the user would have to choose that protein when converting the file.
- artmsQuantification is now able to process full uniprot ids (sp|ENTRY|ENTRY_NAME ...): it detects and extracts the ENTRY id
- Depracated: artmsReplicatePlots

# artMS 1.0.12 (2019-04-07)

- Fix: filtering by group uses the Leading Razor Protein when the filter "protein_groups: remove" is selected
- Fix: `plotSP` of `artmsQualityControlEvidenceExtended` fails with old versions of MaxQuant
- artmsProtein2SiteConversion also processes RefSeq annotations

# artMS 1.0.11 (2019-03-25)

- Fix: column namne MS/MS count for old versions of MaxQuant

# artMS 1.0.10 (2019-02-22)

- Fix: final fix of build problems

# artMS 1.0.9 (2019-02-21)

- Fix: build problem

# artMS 1.0.8 (2019-02-18)

- Fix: Remove abundance files to prevent users' confusions
- Fix: use "Leading Razor Protein" as main protein ID

# artMS 1.0.7 (2019-02-18)

- Bug fix: update annotations to process species without support by bioconductor annotation packages.
- Improve messages to the user

# artMS 1.0.6 (2019-01-21)

- Bug fix (checkfile)

# artMS 1.0.5 (2019-01-20)

- Update the vignette

# artMS 1.0.3 (2019-01-20)

- bug fix: pre-formatted mss file now works

# artMS 1.0.2 (2019-01-19)

- Fixed error: the condition has length (from bioc 3.9)
- Fixed bug with artmsVolcanoPlot

# artMS 1.0.1 (2018-10-29)

- Fixed bug on artmsProtein2SiteConvertion

# artMS 1.0.0 (2018-10-29)

- artMS is officially released

# artMS 0.99.102 (2018-10-27)

- Fix bug in artmsConvertMetabolomics()

# artMS 0.99.101

- Update vignettes

# artMS 0.99.100

- This version addresses all the suggestions and comments raised by the bioconductor's reviewer
- Incorporates functions to handle Metabolomics based on the Markview alignment table

# artMS 0.99.52

- Updated Vignette extensively

# artMS 0.99.25

- Update Vignette to just output html

# artMS 0.99.02

- Add the Bioconductor webhook

# artMS 0.99.01

- Submit the package to the Bioconductor project
