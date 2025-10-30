# artMS Comprehensive Refactoring Plan
## R Package Modernization Roadmap

**Version:** 1.0  
**Date:** 2025-10-30  
**Current Package Version:** 1.10.3  
**Target R Version:** >= 4.1.0  

---

## Executive Summary

This document provides a comprehensive, prioritized refactoring plan for modernizing the artMS R package. The package currently consists of ~14,000 lines of code across 22 R files with 33 exported functions, providing proteomics data analysis workflows in four core functional areas: Quality Control, MSstats-based Quantification, Downstream Analysis, and Export Utilities.

### Current State Analysis

**Strengths:**
- Well-established Bioconductor package with active maintenance
- Comprehensive functionality covering full proteomics workflow
- Integration with MSstats (recently updated for v4.1.0)
- Good documentation structure (roxygen2, vignettes)
- Active user base and clear use cases

**Areas for Improvement:**
- Mixed coding styles (base R, data.table, some tidyverse)
- Dependency overlap (both plyr and dplyr)
- Limited test coverage (only 2 test files)
- Large monolithic functions (some >500 lines)
- Inconsistent function naming conventions
- Code duplication in plotting functions
- Missing input validation in several functions
- No lifecycle badges for function maturity

---

## Phase 1: Critical Updates and Foundation (Priority: HIGH)
**Estimated Duration:** 4-6 weeks  
**Effort:** High  
**Breaking Changes:** Minimal  

### 1.1 Dependency Audit and Modernization
**Priority:** P0 - Critical  
**Effort:** Medium  
**Breaking Changes:** None (if done carefully)

**Current Issues:**
- Both `plyr` and `dplyr` imported (plyr is superseded)
- `gProfileR` suggested but deprecated (replaced by `gprofiler2`)
- Multiple overlapping plotting packages

**Actions:**
1. **Replace plyr with dplyr** (Files affected: `plots.R`, `enrichments.R`, `annotations.R`)
   - Replace `plyr::ddply()` with `dplyr::group_by() %>% summarise()`
   - Replace `plyr::rename()` with `dplyr::rename()`
   - Remove plyr from DESCRIPTION/NAMESPACE
   - **Estimated effort:** 2-3 days

2. **Update gProfileR to gprofiler2**
   - Update `artmsEnrichLog2fc()` and `artmsEnrichProfiler()`
   - Update documentation and examples
   - **Estimated effort:** 1 day

3. **Document version requirements**
   - Create `DEPENDENCIES.md` documenting each dependency and minimum version
   - Document MSstats API version compatibility
   - **Estimated effort:** 1 day

**Files to modify:**
- `DESCRIPTION` (Imports/Suggests sections)
- `R/enrichments.R`
- `R/plots.R`
- `R/annotations.R`
- New: `DEPENDENCIES.md`

### 1.2 Code Style Standardization
**Priority:** P0 - Critical  
**Effort:** Medium  
**Breaking Changes:** None (internal only)

**Actions:**
1. **Adopt tidyverse style guide consistently**
   - Run `styler::style_pkg()` on entire codebase
   - Configure `.lintr` for ongoing compliance
   - **Estimated effort:** 1 day

2. **Standardize function naming**
   - Current inconsistency: `artmsQualityControlEvidenceBasic` vs `artms_data_ph_evidence`
   - Internal functions: ensure all start with `.artms_`
   - Exported functions: keep `artms` prefix (maintain backward compatibility)
   - **Estimated effort:** 1 day (documentation updates needed)

3. **Setup pre-commit hooks**
   - Use `precommit` package
   - Configure for style checking, documentation checks
   - **Estimated effort:** 0.5 day

**Files to modify:**
- All R files (via styler)
- New: `.lintr`
- New: `.pre-commit-config.yaml`

### 1.3 Input Validation and Type Checking
**Priority:** P1 - High  
**Effort:** High  
**Breaking Changes:** None (additions only)

**Actions:**
1. **Add input validation to all exported functions**
   - Use `checkmate` package for efficient validation
   - Validate file paths, data.frame structures, parameter values
   - **Estimated effort:** 1 week

2. **Improve error messages**
   - Make error messages actionable and user-friendly
   - Include expected vs. actual in error messages
   - **Estimated effort:** 3 days (concurrent with validation)

**Example implementation:**
```r
artmsQualityControlEvidenceBasic <- function(evidence_file, keys_file, ...) {
  # Input validation
  checkmate::assert(
    checkmate::check_string(evidence_file),
    checkmate::check_data_frame(evidence_file)
  )
  checkmate::assert_choice(prot_exp, choices = c('AB', 'PH', 'UB', 'AC', 'APMS'))
  
  # Load and validate
  evidence_file <- .artms_checkIfFile(evidence_file)
  .artms_validateEvidenceStructure(evidence_file)
  ...
}
```

**Files to modify:**
- All files with exported functions (22 files)
- `DESCRIPTION` (add checkmate to Imports)

### 1.4 Testing Infrastructure Expansion
**Priority:** P1 - High  
**Effort:** High  
**Breaking Changes:** None

**Actions:**
1. **Create test structure by module**
   ```
   tests/testthat/
     test-qc-basic.R
     test-qc-extended.R
     test-qc-summary.R
     test-quantification.R
     test-downstream-analysis.R
     test-exports-saint.R
     test-exports-other.R
     test-annotations.R
     test-enrichments.R
     test-plots.R
     test-utilities.R
   ```
   - **Estimated effort:** 2 weeks

2. **Create minimal test datasets**
   - Reduce size of existing test data
   - Create focused test data for each module
   - **Estimated effort:** 3 days

3. **Aim for 60%+ code coverage** (Phase 1 target)
   - Focus on exported functions first
   - Use `covr` package for tracking
   - **Estimated effort:** Ongoing

**New files:**
- `tests/testthat/test-*.R` (10+ new test files)
- `tests/testthat/fixtures/` (test data directory)
- Updated: `tests/testthat/test-evidence.R`

---

## Phase 2: Code Organization and Enhancement (Priority: MEDIUM)
**Estimated Duration:** 6-8 weeks  
**Effort:** High  
**Breaking Changes:** Minimal (with deprecation warnings)

### 2.1 Modular File Structure Reorganization
**Priority:** P2 - Medium  
**Effort:** Medium  
**Breaking Changes:** None (internal only)

**Current Structure Issues:**
- `plots.R` contains diverse plotting functions (>400 lines)
- `qualityControlEvidenceBasic.R` is monolithic (>800 lines)
- Utility functions scattered across files

**Proposed New Structure:**
```
R/
  # Core Modules (matching functional areas)
  qc/
    qc-evidence-basic.R
    qc-evidence-extended.R
    qc-summary.R
    qc-metabolomics.R
    qc-plots.R
  
  quantification/
    msstats-runner.R
    msstats-config.R
    msstats-format.R
    msstats-summary.R
  
  analysis/
    enrichment-complexes.R
    enrichment-profiler.R
    clustering.R
    pca.R
  
  plots/
    plot-abundance.R
    plot-heatmap.R
    plot-qc.R
    plot-enrichment.R
  
  exports/
    export-saint-express.R
    export-saint-q.R
    export-photon.R
    export-phosfate.R
  
  utils/
    annotations.R
    file-io.R
    validators.R
    conversions.R
  
  # Keep at root for backward compatibility
  main.R (import/export statements)
  data.R (data documentation)
```

**Migration Strategy:**
1. Create new directory structure
2. Move functions to new locations
3. Update imports in main.R
4. Test that all exports work
5. Update internal documentation
6. **Estimated effort:** 1 week

**Note:** R package structure doesn't support subdirectories in R/, so we'll use file naming convention instead:
- `qc-evidence-basic.R`, `qc-evidence-extended.R`, etc.
- `export-saint-express.R`, `export-saint-q.R`, etc.

### 2.2 Function Decomposition and DRY Principles
**Priority:** P2 - Medium  
**Effort:** High  
**Breaking Changes:** None (internal refactoring)

**Actions:**
1. **Break down monolithic functions**
   - Target: No function >200 lines
   - Extract plotting logic into separate internal functions
   - Example: `artmsQualityControlEvidenceBasic()` → multiple internal helpers

2. **Identify and eliminate code duplication**
   - Plotting code (themes, color palettes)
   - Data transformation patterns
   - File I/O operations
   - **Estimated effort:** 2 weeks

3. **Create helper function library**
   - Common plotting themes
   - Data validation helpers
   - File path handling
   - **Estimated effort:** 1 week

**Example Refactoring:**
```r
# Before: monolithic function
artmsQualityControlEvidenceBasic <- function(...) {
  # 800+ lines of code
}

# After: decomposed
artmsQualityControlEvidenceBasic <- function(...) {
  .artms_qc_validate_inputs(evidence_file, keys_file, ...)
  evidence_data <- .artms_qc_load_and_prepare(evidence_file, keys_file)
  
  if (plotINTDIST) .artms_qc_plot_intensity_distribution(evidence_data, ...)
  if (plotREPRO) .artms_qc_plot_reproducibility(evidence_data, ...)
  if (plotCORMAT) .artms_qc_plot_correlation_matrix(evidence_data, ...)
  ...
}
```

### 2.3 Plotting System Modernization
**Priority:** P2 - Medium  
**Effort:** Medium  
**Breaking Changes:** None (output identical)

**Actions:**
1. **Standardize on ggplot2**
   - Already primary plotting system
   - Remove or minimize base R plotting
   - Create consistent theme system
   - **Estimated effort:** 1 week

2. **Separate plotting from computation**
   - Functions should return data + optional plotting
   - Enable programmatic use without PDF generation
   - **Estimated effort:** 1 week

3. **Create plot object return option**
   ```r
   artmsQualityControlEvidenceBasic(..., return_plots = FALSE)
   # If TRUE, returns list of ggplot objects instead of printing PDFs
   ```
   - **Estimated effort:** 3 days

**Files affected:**
- `plots.R` → multiple `plot-*.R` files
- All QC functions
- `analysisQuantifications.R`

### 2.4 Progress Indicators and User Experience
**Priority:** P2 - Medium  
**Effort:** Low  
**Breaking Changes:** None

**Actions:**
1. **Add progress bars for long operations**
   - Use `cli` package progress bars
   - Target functions: QC functions, MSstats runs, enrichment analysis
   - **Estimated effort:** 3 days

2. **Improve messaging system**
   - Consistent use of `cli::cli_alert_*()` functions
   - Different levels: info, success, warning, error
   - **Estimated effort:** 2 days

3. **Better verbose control**
   - Respect `verbose` parameter consistently
   - Add quiet mode option
   - **Estimated effort:** 2 days

**Example:**
```r
if (verbose) cli::cli_alert_info("Processing {nrow(data)} features...")
pb <- cli::cli_progress_bar("Quality control analysis", total = n_steps)
for (step in steps) {
  cli::cli_progress_update()
  ...
}
cli::cli_alert_success("QC analysis complete!")
```

---

## Phase 3: Advanced Features and Optimization (Priority: LOW)
**Estimated Duration:** 4-6 weeks  
**Effort:** Medium  
**Breaking Changes:** None

### 3.1 Performance Optimization
**Priority:** P3 - Low  
**Effort:** Medium  
**Breaking Changes:** None

**Actions:**
1. **Profile computational bottlenecks**
   - Use `profvis` to identify slow operations
   - Focus on QC functions and large dataset processing
   - **Estimated effort:** 1 week

2. **Implement parallel processing**
   - Use `future` + `furrr` for parallel operations
   - Target: replicate-level operations, protein-level plots
   - Make it optional with `parallel = FALSE` default
   - **Estimated effort:** 1 week

3. **Optimize memory usage**
   - Use data.table efficiently (in-place operations)
   - Stream large file reading where possible
   - **Estimated effort:** 1 week

**Example:**
```r
artmsDataPlots <- function(..., parallel = FALSE) {
  if (parallel) {
    future::plan("multisession")
    plots <- furrr::future_map(unique_proteins, .artms_plot_protein, ...)
  } else {
    plots <- purrr::map(unique_proteins, .artms_plot_protein, ...)
  }
}
```

### 3.2 Enhanced Return Object Structures
**Priority:** P3 - Low  
**Effort:** Medium  
**Breaking Changes:** Controlled (with deprecation)

**Actions:**
1. **Create S3 classes for results**
   ```r
   artms_qc_result
   artms_quantification_result
   artms_enrichment_result
   ```
   - **Estimated effort:** 1 week

2. **Implement print/summary methods**
   - User-friendly display of results
   - Extract key metrics easily
   - **Estimated effort:** 3 days

3. **Enable method chaining**
   - Return objects that can be piped
   - Support tidyverse workflows
   - **Estimated effort:** 1 week

### 3.3 Modern R Practices Integration
**Priority:** P3 - Low  
**Effort:** Medium  
**Breaking Changes:** None

**Actions:**
1. **Add pipe operator support**
   - Already imports from tidyr, can use magrittr
   - Ensure functions are pipe-friendly
   - **Estimated effort:** 3 days

2. **Implement tidy evaluation where appropriate**
   - For functions accepting column names
   - Use `rlang` for NSE handling
   - **Estimated effort:** 1 week

3. **Add lifecycle badges**
   - Use `lifecycle` package
   - Mark stable/experimental/deprecated functions
   - **Estimated effort:** 2 days

**Example:**
```r
#' @lifecycle stable
artmsQualityControlEvidenceBasic <- function(...) { ... }

#' @lifecycle deprecated
#' @description Use `artmsEnrichProfiler()` with gprofiler2 instead
artmsEnrichLog2fc <- function(...) {
  lifecycle::deprecate_warn("1.12.0", "artmsEnrichLog2fc()", "artmsEnrichProfiler()")
  ...
}
```

---

## Phase 4: Documentation and Vignettes (Priority: MEDIUM)
**Estimated Duration:** 3-4 weeks  
**Effort:** Medium  
**Breaking Changes:** None

### 4.1 Enhanced Function Documentation
**Priority:** P2 - Medium  
**Effort:** Medium  
**Breaking Changes:** None

**Actions:**
1. **Expand roxygen2 documentation**
   - Add more detailed examples
   - Document all parameters thoroughly
   - Add return value structure details
   - **Estimated effort:** 2 weeks

2. **Add @family tags for grouping**
   - Group by functional module
   - Easier navigation in help
   - **Estimated effort:** 1 day

3. **Create comprehensive examples**
   - Real-world use cases
   - Show integration between functions
   - **Estimated effort:** 1 week

### 4.2 Vignette Expansion
**Priority:** P2 - Medium  
**Effort:** Medium  
**Breaking Changes:** None

**Actions:**
1. **Create module-specific vignettes**
   - `vignette-qc.Rmd` - Quality control workflows
   - `vignette-quantification.Rmd` - MSstats integration
   - `vignette-downstream.Rmd` - Enrichment, clustering, PCA
   - `vignette-exports.Rmd` - SAINT/Photon/Phosfate integration
   - **Estimated effort:** 2 weeks

2. **Create quick start guide**
   - 5-minute introduction vignette
   - Common use cases
   - **Estimated effort:** 2 days

3. **Add troubleshooting guide**
   - Common errors and solutions
   - FAQ section
   - **Estimated effort:** 1 week

---

## Backward Compatibility Strategy

### Deprecation Policy

**Phase 1: Warning (v1.12.0 - v1.14.0)**
- Functions marked as deprecated with lifecycle badges
- Warning messages guide users to new alternatives
- Old functions continue to work

**Phase 2: Defunct (v1.16.0)**
- Deprecated functions throw errors with clear migration path
- Documentation updated to show only current functions

**Phase 3: Removal (v2.0.0)**
- Major version bump allows complete removal
- Migration guide provided

### Migration Guide Contents

1. **Function Replacements**
   - plyr → dplyr migration examples
   - gProfileR → gprofiler2 updates
   - Any renamed functions

2. **API Changes**
   - New parameter defaults
   - Changed return structures
   - New required parameters

3. **Code Examples**
   - Before/after for common operations
   - Update scripts for users

---

## File Structure Recommendations

### Proposed Final Structure

```
artMS/
├── .github/
│   └── workflows/
│       ├── R-CMD-check.yaml
│       └── test-coverage.yaml
├── R/
│   ├── analysis-*.R (5 files)
│   ├── export-*.R (4 files)
│   ├── plot-*.R (4 files)
│   ├── qc-*.R (4 files)
│   ├── quantification-*.R (4 files)
│   ├── utils-*.R (5 files)
│   ├── data.R
│   └── main.R
├── data/
│   └── *.rda (existing data objects)
├── inst/
│   ├── CITATION
│   ├── extdata/ (example data)
│   └── scripts/ (helper scripts)
├── man/
│   └── *.Rd (generated documentation)
├── tests/
│   ├── testthat/
│   │   ├── fixtures/ (test data)
│   │   ├── helper-*.R (test utilities)
│   │   └── test-*.R (11+ test files)
│   └── testthat.R
├── vignettes/
│   ├── artMS-quickstart.Rmd
│   ├── artMS-qc.Rmd
│   ├── artMS-quantification.Rmd
│   ├── artMS-downstream.Rmd
│   ├── artMS-exports.Rmd
│   └── artMS-troubleshooting.Rmd
├── .lintr
├── .pre-commit-config.yaml
├── DEPENDENCIES.md
├── DESCRIPTION
├── LICENSE
├── MIGRATION_GUIDE.md
├── NAMESPACE
├── NEWS.md (convert from NEWS)
├── README.md
└── REFACTORING_PLAN.md (this document)
```

---

## Priority Rankings Summary

### P0 - Critical (Start Immediately)
1. Dependency modernization (plyr → dplyr, gProfileR → gprofiler2)
2. Code style standardization
3. Input validation and error messages

### P1 - High (Phase 1)
4. Test coverage expansion (target 60%+)
5. Setup CI/CD improvements

### P2 - Medium (Phase 2)
6. File organization and naming
7. Function decomposition
8. Plotting system improvements
9. Progress indicators
10. Documentation expansion

### P3 - Low (Phase 3)
11. Performance optimization
12. Parallel processing
13. S3 classes and methods
14. Advanced tidyverse integration

---

## Effort Estimates by Task

### High Effort (2-4 weeks)
- Test coverage expansion
- Function decomposition and DRY refactoring
- Input validation across all functions
- Vignette expansion

### Medium Effort (1-2 weeks)
- Dependency updates and testing
- File organization
- Plotting system modernization
- Performance profiling and optimization
- Enhanced documentation

### Low Effort (<1 week)
- Code style standardization (automated)
- Progress indicators
- Lifecycle badges
- Naming conventions
- Migration guide creation

---

## Success Metrics

### Code Quality
- [ ] All functions have input validation
- [ ] Test coverage >60% (Phase 1), >80% (Phase 3)
- [ ] Zero linter warnings
- [ ] All exported functions documented with examples
- [ ] No deprecated dependencies

### User Experience
- [ ] Consistent error messages
- [ ] Progress indicators on operations >30 seconds
- [ ] All long-running functions have verbose control
- [ ] Return objects have print methods

### Performance
- [ ] 20% improvement in QC function runtime (target)
- [ ] Support for parallel processing where beneficial
- [ ] Memory usage profiled and optimized

### Documentation
- [ ] 6 comprehensive vignettes
- [ ] Migration guide for users
- [ ] Troubleshooting guide
- [ ] All functions have @family tags

---

## Risk Assessment

### Low Risk
- Code styling (automated, no logic changes)
- Documentation improvements
- Adding tests
- Input validation additions

### Medium Risk
- Dependency updates (requires thorough testing)
- File reorganization (needs careful NAMESPACE management)
- Function decomposition (must maintain exact behavior)

### High Risk (Requires Extra Caution)
- Changing function signatures
- Modifying return structures
- Removing deprecated functions

### Mitigation Strategies
1. Comprehensive testing before each phase
2. Staged rollout with deprecation warnings
3. Beta testing with known users
4. Clear migration documentation
5. Version number communication (semantic versioning)

---

## Implementation Timeline

### Month 1-2: Phase 1 Foundation
- Week 1-2: Dependency updates and testing
- Week 3-4: Code style and validation
- Week 5-8: Test infrastructure expansion

### Month 3-4: Phase 2 Organization
- Week 9-11: File reorganization and function decomposition
- Week 12-14: Plotting improvements and UX enhancements
- Week 15-16: Documentation expansion

### Month 5-6: Phase 3 Advanced Features
- Week 17-19: Performance optimization
- Week 20-22: S3 classes and advanced features
- Week 23-24: Final testing and documentation

---

## Conclusion

This refactoring plan provides a structured approach to modernizing artMS while maintaining stability and backward compatibility. The phased approach allows for:

1. **Immediate value** through critical updates (Phase 1)
2. **Improved maintainability** through better organization (Phase 2)
3. **Enhanced capabilities** through modern R features (Phase 3)
4. **Better user experience** throughout all phases

The plan prioritizes stability and backward compatibility while systematically improving code quality, test coverage, documentation, and performance. Each phase builds on the previous one, allowing for iterative improvement and validation.

**Recommended Start:** Begin with Phase 1 immediately, focusing on dependency updates and input validation as the foundation for all subsequent improvements.
