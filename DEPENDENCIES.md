# artMS Package Dependencies

**Last Updated:** 2025-10-30  
**Package Version:** 1.10.3  

## Overview

This document provides detailed information about artMS dependencies, their purposes, minimum version requirements, and modernization recommendations.

---

## Critical Dependencies

### MSstats (>= 4.0.0)
- **Purpose:** Core quantification engine for relative protein abundance
- **Minimum Version:** 4.0.0 (significant API changes)
- **Current Usage:** Used in `artmsQuantification()`, `artmsRunMSstats()`
- **Status:** ✅ Recently updated (v1.10.1)
- **Notes:** Major API changes in v4.0 required R >= 4.1.0

### data.table (>= 1.14.0)
- **Purpose:** Fast data manipulation, primary data structure
- **Minimum Version:** 1.14.0 recommended
- **Current Usage:** Throughout package for data processing
- **Status:** ✅ Active and maintained
- **Notes:** Critical for performance with large datasets

---

## Visualization Dependencies

### ggplot2 (>= 3.3.0)
- **Purpose:** Primary plotting system
- **Minimum Version:** 3.3.0
- **Current Usage:** All quality control plots, data plots
- **Status:** ✅ Active and maintained
- **Notes:** Core tidyverse package

### pheatmap (>= 1.0.12)
- **Purpose:** Heatmap generation
- **Current Usage:** Enrichment and clustering heatmaps
- **Status:** ✅ Stable
- **Alternative:** Consider ComplexHeatmap (already in Suggests)

### plotly (>= 4.9.0)
- **Purpose:** Interactive plots
- **Current Usage:** Interactive visualizations in analysis functions
- **Status:** ✅ Active
- **Notes:** Name conflicts managed via rawNamespace import

### ggrepel (>= 0.9.0)
- **Purpose:** Non-overlapping text labels in plots
- **Current Usage:** QC plots, PCA plots
- **Status:** ✅ Active

### ggdendro (>= 0.1.20)
- **Purpose:** Dendrogram support for ggplot2
- **Current Usage:** Clustering visualizations
- **Status:** ✅ Maintained

### VennDiagram (>= 1.6.20)
- **Purpose:** Venn diagram generation
- **Current Usage:** Comparison visualizations
- **Status:** ⚠️ Last updated 2018, consider alternatives
- **Alternative:** ggvenn (more modern, ggplot2-based)

### UpSetR (>= 1.4.0)
- **Purpose:** Set intersection plots
- **Current Usage:** Complex set comparisons
- **Status:** ✅ Stable
- **Alternative:** ggupset (tidyverse integration)

### circlize (>= 0.4.0)
- **Purpose:** Circular visualization
- **Current Usage:** Specific data visualizations
- **Status:** ✅ Active

### gplots (>= 3.1.0)
- **Purpose:** Various plots, specifically heatmap.2
- **Current Usage:** `heatmap.2()` function
- **Status:** ⚠️ Maintenance mode
- **Alternative:** pheatmap or ComplexHeatmap already used

### corrplot (>= 0.84)
- **Purpose:** Correlation matrix visualization
- **Current Usage:** QC correlation matrices
- **Status:** ✅ Active

### RColorBrewer (>= 1.1-2)
- **Purpose:** Color palettes
- **Current Usage:** All plotting functions
- **Status:** ✅ Stable

---

## Data Manipulation Dependencies

### dplyr (>= 1.0.0)
- **Purpose:** Data manipulation (tidyverse)
- **Minimum Version:** 1.0.0 (major API update)
- **Current Usage:** Data filtering, summarization
- **Status:** ✅ Active, core tidyverse
- **Notes:** Preferred over plyr for new code

### tidyr (>= 1.1.0)
- **Purpose:** Data reshaping (pivot operations)
- **Minimum Version:** 1.1.0 (pivot_* functions)
- **Current Usage:** Data transformation in analysis functions
- **Status:** ✅ Active, core tidyverse

### stringr (>= 1.4.0)
- **Purpose:** String manipulation
- **Current Usage:** Text processing throughout
- **Status:** ✅ Active, core tidyverse

### plyr (>= 1.8.6)
- **Purpose:** Split-apply-combine (DEPRECATED)
- **Current Usage:** `ddply()`, `rename()`, `summarise()`
- **Status:** ⚠️ SUPERSEDED by dplyr
- **Recommendation:** **REMOVE in Phase 1**
- **Migration:** Replace with dplyr equivalents
- **Files to update:** `plots.R`, `enrichments.R`, `annotations.R`

---

## Annotation and Enrichment Dependencies

### AnnotationDbi (>= 1.50.0)
- **Purpose:** Bioconductor annotation infrastructure
- **Current Usage:** Annotation functions
- **Status:** ✅ Active, Bioconductor core

### org.Hs.eg.db (>= 3.11.0)
- **Purpose:** Human genome annotation
- **Current Usage:** Human proteome annotations
- **Status:** ✅ Updated regularly (Bioconductor release cycle)

### org.Mm.eg.db (Suggested)
- **Purpose:** Mouse genome annotation
- **Current Usage:** Mouse proteome annotations
- **Status:** ✅ Updated regularly (Bioconductor release cycle)

### limma (>= 3.44.0)
- **Purpose:** Statistical analysis for genomics
- **Current Usage:** Statistical functions
- **Status:** ✅ Active, Bioconductor core

---

## File I/O Dependencies

### openxlsx (>= 4.2.0)
- **Purpose:** Excel file writing (no Java dependency)
- **Current Usage:** Export functions
- **Status:** ✅ Active
- **Notes:** Preferred over xlsx (no rJava requirement)

### yaml (>= 2.2.0)
- **Purpose:** Configuration file parsing
- **Current Usage:** `artmsWriteConfigYamlFile()`, config reading
- **Status:** ✅ Active

---

## Statistical and Analysis Dependencies

### cluster (>= 2.1.0)
- **Purpose:** Clustering algorithms (PAM)
- **Current Usage:** Clustering analysis
- **Status:** ✅ Base R recommended package

### stats (Base R)
- **Purpose:** Statistical functions
- **Current Usage:** Throughout (correlation, PCA, etc.)
- **Status:** ✅ Core R

### scales (>= 1.1.0)
- **Purpose:** Scale functions for visualization
- **Current Usage:** Plot formatting
- **Status:** ✅ Active

---

## Utility Dependencies

### bit64 (>= 4.0.0)
- **Purpose:** 64-bit integer support
- **Current Usage:** Large integer handling in MaxQuant files
- **Status:** ✅ Active
- **Notes:** Essential for accurate data import

### getopt (>= 1.20.3)
- **Purpose:** Command-line argument parsing
- **Current Usage:** Command-line interface support
- **Status:** ✅ Stable

### seqinr (>= 3.6-1)
- **Purpose:** Sequence analysis
- **Current Usage:** Protein sequence processing
- **Status:** ✅ Active

---

## Suggested Dependencies (Optional)

### gProfileR (>= 0.7.0)
- **Purpose:** Gene enrichment analysis (DEPRECATED)
- **Current Usage:** `artmsEnrichLog2fc()`, `artmsEnrichProfiler()`
- **Status:** ⚠️ **DEPRECATED** (archived on CRAN)
- **Replacement:** **gprofiler2** (>= 0.2.0)
- **Recommendation:** **UPDATE in Phase 1**
- **Priority:** HIGH - package is deprecated

### gprofiler2 (>= 0.2.0)
- **Purpose:** Gene enrichment analysis (modern replacement)
- **Status:** ✅ Active replacement for gProfileR
- **Recommendation:** Add to Suggests, migrate enrichment functions

### ComplexHeatmap (>= 2.4.0)
- **Purpose:** Advanced heatmap visualization
- **Current Usage:** Optional for `artmsAnalysisQuantifications()`
- **Status:** ✅ Active, Bioconductor

### factoextra (>= 1.0.7)
- **Purpose:** PCA/clustering visualization
- **Current Usage:** Optional for `artmsAnalysisQuantifications()`
- **Status:** ✅ Active

### FactoMineR (>= 2.3)
- **Purpose:** Multivariate analysis
- **Current Usage:** Optional for `artmsAnalysisQuantifications()`
- **Status:** ✅ Active

### PerformanceAnalytics (>= 2.0.4)
- **Purpose:** Performance and risk analysis
- **Current Usage:** Optional correlation analysis
- **Status:** ✅ Active

### BiocStyle (>= 2.16.0)
- **Purpose:** Vignette styling
- **Current Usage:** Vignette compilation
- **Status:** ✅ Active, Bioconductor

### testthat (>= 3.0.0)
- **Purpose:** Unit testing framework
- **Minimum Version:** 3.0.0 (3rd edition recommended)
- **Current Usage:** Package testing
- **Status:** ✅ Active
- **Recommendation:** Expand test coverage

### knitr (>= 1.30)
- **Purpose:** Vignette compilation
- **Current Usage:** Vignette building
- **Status:** ✅ Active

### rmarkdown (>= 2.5)
- **Purpose:** R Markdown support
- **Current Usage:** Vignette building
- **Status:** ✅ Active

---

## Recommended Additions

### checkmate (>= 2.0.0)
- **Purpose:** Fast and consistent input validation
- **Status:** ✅ Active
- **Recommendation:** Add to Imports (Phase 1)
- **Benefit:** Standardized input validation across all functions

### cli (>= 3.0.0)
- **Purpose:** Modern console output and progress bars
- **Status:** ✅ Active
- **Recommendation:** Add to Imports (Phase 2)
- **Benefit:** Better user feedback and progress indication

### lifecycle (>= 1.0.0)
- **Purpose:** Function lifecycle management
- **Status:** ✅ Active
- **Recommendation:** Add to Imports (Phase 2/3)
- **Benefit:** Clear communication of function maturity and deprecation

### rlang (>= 0.4.0)
- **Purpose:** Tidy evaluation support
- **Status:** ✅ Active (already imported via dplyr)
- **Recommendation:** Use explicitly for NSE if needed (Phase 3)

---

## Dependency Modernization Priorities

### Phase 1: Critical Updates

1. **Remove plyr** → Use dplyr exclusively
   - Effort: Medium
   - Impact: Code simplification, reduced conflicts
   - Files: 3-4 R files

2. **Replace gProfileR → gprofiler2**
   - Effort: Low-Medium
   - Impact: Continued functionality (gProfileR is deprecated)
   - Files: `enrichments.R`

3. **Add checkmate**
   - Effort: Low (addition only)
   - Impact: Better input validation
   - Files: All exported functions

### Phase 2: Enhancements

4. **Add cli**
   - Effort: Low
   - Impact: Better UX
   - Files: Long-running functions

5. **Add lifecycle**
   - Effort: Low
   - Impact: Better deprecation management
   - Files: Documentation updates

6. **Consider ggvenn** (replace VennDiagram)
   - Effort: Low
   - Impact: More modern, consistent plotting
   - Files: `plots.R`

### Phase 3: Optional Optimizations

7. **Evaluate ComplexHeatmap** for all heatmaps
   - Effort: Medium
   - Impact: More powerful heatmaps
   - Files: `plots.R`, analysis functions

8. **Consider ggupset** (replace UpSetR)
   - Effort: Medium
   - Impact: Better tidyverse integration
   - Files: Comparison visualization functions

---

## Version Compatibility Matrix

| Dependency | Min Version | Tested Version | R Version Required |
|------------|-------------|----------------|-------------------|
| R | 4.1.0 | 4.3.0 | - |
| MSstats | 4.0.0 | 4.6.0 | >= 4.1.0 |
| data.table | 1.14.0 | 1.14.8 | >= 3.1.0 |
| ggplot2 | 3.3.0 | 3.4.2 | >= 3.3.0 |
| dplyr | 1.0.0 | 1.1.2 | >= 3.4.0 |
| tidyr | 1.1.0 | 1.3.0 | >= 3.4.0 |

---

## Bioconductor Release Compatibility

- **Current:** Bioconductor 3.13+ (R >= 4.1.0)
- **Tested:** Bioconductor 3.17 (R >= 4.3.0)
- **Target:** Compatible with Bioc release and devel

---

## Dependency Health Check

### ✅ Healthy (Active Development)
- MSstats, ggplot2, dplyr, tidyr, data.table
- Bioconductor annotation packages
- Core visualization packages (pheatmap, plotly)

### ⚠️ Needs Attention
- **plyr** - Superseded, should be removed
- **gProfileR** - Deprecated, must be replaced
- **VennDiagram** - Inactive, consider alternatives
- **gplots** - Maintenance mode, alternatives available

### ❌ Critical Issues
- None currently, but gProfileR deprecation is urgent

---

## Testing Strategy for Dependencies

1. **Before updating any dependency:**
   - Run full test suite
   - Check for deprecation warnings
   - Review NEWS/changelog for breaking changes

2. **After updating:**
   - Run `R CMD check` 
   - Test all examples
   - Run full test suite
   - Test on multiple R versions (oldrel, release, devel)

3. **Continuous monitoring:**
   - Subscribe to dependency release announcements
   - Regular `devtools::install_deps(upgrade = "ask")`
   - Monitor Bioconductor release schedules

---

## License Compatibility

All dependencies are compatible with artMS GPL (>= 3) license:
- GPL-2, GPL-3: Compatible
- MIT, BSD: Compatible
- Artistic-2.0 (Bioconductor): Compatible

---

## Future Considerations

1. **Reducing dependency count:**
   - Current: 39 Imports + 11 Suggests = 50 total
   - Consider consolidating where possible
   - Evaluate rarely-used dependencies

2. **Alternative packages to monitor:**
   - `gprofiler2` for enrichment
   - `ggvenn` for Venn diagrams
   - `patchwork` for plot composition
   - `progressr` for unified progress reporting

3. **Staying current:**
   - Annual dependency audit
   - Follow tidyverse and Bioconductor updates
   - Monitor MSstats API changes carefully
