# artMS Refactoring Implementation Guide
## Developer Quick Reference

**For:** Developers implementing the refactoring plan  
**Last Updated:** 2025-10-30  

---

## Quick Start

This document provides practical guidance for developers implementing the refactoring plan. Read `REFACTORING_PLAN.md` first for the overall strategy.

---

## Phase 1 Implementation Tasks

### Task 1: Replace plyr with dplyr

**Priority:** P0  
**Estimated Time:** 2-3 days  

**Files to Update:**
- `R/plots.R` (lines with `plyr::ddply`, `plyr::rename`)
- `R/enrichments.R` (lines with `plyr::summarise`)
- `R/annotations.R` (potential plyr usage)

**Common Replacements:**

```r
# OLD: plyr::ddply
result <- plyr::ddply(data, .(group, condition), summarise, 
                     mean_val = mean(value),
                     sd_val = sd(value))

# NEW: dplyr
result <- data %>%
  dplyr::group_by(group, condition) %>%
  dplyr::summarise(mean_val = mean(value),
                   sd_val = sd(value),
                   .groups = 'drop')

# OLD: plyr::rename
data <- plyr::rename(data, c("old_name" = "new_name"))

# NEW: dplyr
data <- dplyr::rename(data, new_name = old_name)  # Note: order reversed!
```

**Testing:**
```r
# Run after each file update
devtools::test()
devtools::check()
```

**DESCRIPTION Update:**
```r
# Remove from Imports:
# plyr,

# Ensure in Imports:
# dplyr (>= 1.0.0),
```

---

### Task 2: Update gProfileR to gprofiler2

**Priority:** P0  
**Estimated Time:** 1 day  

**Files to Update:**
- `R/enrichments.R` - `artmsEnrichLog2fc()`, `artmsEnrichProfiler()`

**API Changes:**

```r
# OLD: gProfileR
library(gProfileR)
result <- gprofiler(query = genes,
                   organism = "hsapiens",
                   src_filter = c("GO:BP", "KEGG", "REAC"))

# NEW: gprofiler2
library(gprofiler2)
result <- gost(query = genes,
              organism = "hsapiens",
              sources = c("GO:BP", "KEGG", "REAC"))
```

**Species Code Mapping:**
```r
# Create internal helper function
.artms_map_species_to_gprofiler2 <- function(species) {
  species_map <- c(
    "human" = "hsapiens",
    "mouse" = "mmusculus"
  )
  
  if (species %in% names(species_map)) {
    return(species_map[species])
  } else if (species %in% species_map) {
    return(species)  # Already in gprofiler2 format
  } else {
    stop("Species '", species, "' not supported")
  }
}
```

**Update Function:**
```r
artmsEnrichLog2fc <- function(dataset, species, background, ...) {
  # Check for gprofiler2
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop("Package 'gprofiler2' required. Install with: install.packages('gprofiler2')")
  }
  
  # Deprecation warning
  lifecycle::deprecate_warn(
    "1.12.0",
    "artmsEnrichLog2fc()",
    "artmsEnrichProfiler()"
  )
  
  # Convert species
  species_code <- .artms_map_species_to_gprofiler2(species)
  
  # Call gprofiler2
  result <- gprofiler2::gost(
    query = dataset$Gene,
    organism = species_code,
    custom_bg = background,
    ...
  )
  
  return(result)
}
```

**DESCRIPTION Update:**
```r
# Remove from Suggests:
# gProfileR,

# Add to Suggests:
# gprofiler2,
```

---

### Task 3: Add Input Validation with checkmate

**Priority:** P1  
**Estimated Time:** 1 week (all functions)  

**Add to DESCRIPTION:**
```r
Imports:
  checkmate,
  # ... other imports
```

**Standard Validation Patterns:**

```r
#' @importFrom checkmate assert_string assert_data_frame assert_choice
#' assert_file_exists assert_number assert_logical

artmsFunctionName <- function(evidence_file, 
                             keys_file,
                             prot_exp = c('AB', 'PH', 'UB', 'AC', 'APMS'),
                             output_dir = "output",
                             verbose = TRUE) {
  
  # Validate inputs
  checkmate::assert(
    checkmate::check_string(evidence_file),
    checkmate::check_data_frame(evidence_file),
    .var.name = "evidence_file",
    add = "Either a file path or data.frame required"
  )
  
  checkmate::assert(
    checkmate::check_string(keys_file),
    checkmate::check_data_frame(keys_file),
    .var.name = "keys_file"
  )
  
  prot_exp <- match.arg(prot_exp)
  checkmate::assert_string(output_dir)
  checkmate::assert_logical(verbose, len = 1)
  
  # If file paths, check they exist
  if (is.character(evidence_file)) {
    checkmate::assert_file_exists(evidence_file, access = "r")
  }
  if (is.character(keys_file)) {
    checkmate::assert_file_exists(keys_file, access = "r")
  }
  
  # Continue with function logic
  ...
}
```

**Validation Helpers to Create:**

```r
# R/utils-validators.R

#' Validate evidence file structure
#' @keywords internal
.artms_validate_evidence_structure <- function(evidence_df) {
  required_cols <- c("Proteins", "Modified.sequence", "Charge", "Intensity")
  
  missing_cols <- setdiff(required_cols, colnames(evidence_df))
  if (length(missing_cols) > 0) {
    stop("Evidence file missing required columns: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  invisible(TRUE)
}

#' Validate keys file structure
#' @keywords internal
.artms_validate_keys_structure <- function(keys_df) {
  required_cols <- c("RawFile", "Condition", "BioReplicate", "Run")
  
  missing_cols <- setdiff(required_cols, colnames(keys_df))
  if (length(missing_cols) > 0) {
    stop("Keys file missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  invisible(TRUE)
}

#' Validate contrast file structure
#' @keywords internal
.artms_validate_contrast_structure <- function(contrast_df) {
  required_cols <- c("Condition1", "Condition2")
  
  missing_cols <- setdiff(required_cols, colnames(contrast_df))
  if (length(missing_cols) > 0) {
    stop("Contrast file missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  invisible(TRUE)
}
```

---

### Task 4: Expand Test Coverage

**Priority:** P1  
**Estimated Time:** 2 weeks  

**Test File Structure:**

Create these test files in `tests/testthat/`:

```
test-qc-basic.R          # Quality control basic functions
test-qc-extended.R       # Quality control extended functions
test-quantification.R    # MSstats integration
test-enrichment.R        # Enrichment analysis
test-plots.R             # Plotting functions
test-export-saint.R      # SAINT export functions
test-annotations.R       # Annotation functions
test-validators.R        # Input validation
test-utils.R             # Utility functions
```

**Test Template:**

```r
# tests/testthat/test-qc-basic.R

context("Quality Control - Basic Functions")

# Setup test data
setup({
  test_evidence <- artms_data_ph_evidence
  test_keys <- artms_data_ph_keys
})

test_that("artmsQualityControlEvidenceBasic validates inputs correctly", {
  # Test with NULL inputs
  expect_error(
    artmsQualityControlEvidenceBasic(
      evidence_file = NULL,
      keys_file = NULL
    ),
    "evidence_file"
  )
  
  # Test with invalid prot_exp
  expect_error(
    artmsQualityControlEvidenceBasic(
      evidence_file = test_evidence,
      keys_file = test_keys,
      prot_exp = "INVALID"
    ),
    "must be element of set"
  )
})

test_that("artmsQualityControlEvidenceBasic processes data correctly", {
  # Test with valid inputs
  result <- artmsQualityControlEvidenceBasic(
    evidence_file = test_evidence,
    keys_file = test_keys,
    prot_exp = "PH",
    printPDF = FALSE,
    verbose = FALSE
  )
  
  # Result checks
  expect_true(is.list(result) || is.null(result))
})

test_that("artmsQualityControlEvidenceBasic handles file inputs", {
  # Create temp files
  temp_evidence <- tempfile(fileext = ".txt")
  temp_keys <- tempfile(fileext = ".txt")
  
  write.table(test_evidence, temp_evidence, sep = "\t", row.names = FALSE)
  write.table(test_keys, temp_keys, sep = "\t", row.names = FALSE)
  
  # Test with file paths
  expect_silent(
    artmsQualityControlEvidenceBasic(
      evidence_file = temp_evidence,
      keys_file = temp_keys,
      prot_exp = "PH",
      printPDF = FALSE,
      verbose = FALSE
    )
  )
  
  # Cleanup
  unlink(c(temp_evidence, temp_keys))
})
```

**Run Tests:**
```r
devtools::test()                    # Run all tests
devtools::test_file("test-qc-basic.R")  # Run specific file
covr::package_coverage()            # Check coverage
```

**Coverage Goals:**
- Phase 1: 60%+ overall
- Focus on exported functions first
- Internal functions: test indirectly through exports

---

## Phase 2 Implementation Tasks

### Task 5: File Reorganization

**Priority:** P2  
**Estimated Time:** 1 week  

**Important:** R packages cannot have subdirectories in `R/`, so use naming conventions:

**Current → New Naming:**
```
R/qualityControlEvidenceBasic.R    → R/qc-evidence-basic.R
R/qualityControlEvidenceExtended.R → R/qc-evidence-extended.R
R/qualityControlSummaryExtended.R  → R/qc-summary-extended.R
R/evidenceToSAINTqFormat.R         → R/export-saint-q.R
R/evidenceToSaintExpressFormat.R   → R/export-saint-express.R
R/plots.R                          → R/plot-general.R
                                     R/plot-qc.R
                                     R/plot-heatmap.R
R/enrichments.R                    → R/analysis-enrichment.R
R/MSstats_functions.R              → R/quantification-helpers.R
R/runMSstats.R                     → R/quantification-runner.R
```

**Migration Process:**

1. Create new file with new name
2. Move functions to new file
3. Update documentation in new file
4. Test that exports work
5. Delete old file
6. Update `R/main.R` if needed
7. Run `devtools::document()` to update NAMESPACE
8. Run `devtools::check()`

**Example:**
```bash
# Step by step
git mv R/qualityControlEvidenceBasic.R R/qc-evidence-basic.R
# Edit R/qc-evidence-basic.R to ensure all functions are properly documented
devtools::document()
devtools::check()
git commit -m "Reorganize: rename qualityControlEvidenceBasic.R to qc-evidence-basic.R"
```

---

### Task 6: Break Down Monolithic Functions

**Priority:** P2  
**Estimated Time:** 2 weeks  

**Target:** Functions >200 lines should be decomposed

**Example Decomposition:**

```r
# Before: One large function (800+ lines)
artmsQualityControlEvidenceBasic <- function(...) {
  # Input validation (50 lines)
  # Data loading (100 lines)
  # Plot 1 generation (100 lines)
  # Plot 2 generation (100 lines)
  # Plot 3 generation (100 lines)
  # etc.
}

# After: Main function + helpers
artmsQualityControlEvidenceBasic <- function(evidence_file,
                                            keys_file,
                                            prot_exp,
                                            output_dir = "qc_basic",
                                            plotINTDIST = FALSE,
                                            plotREPRO = FALSE,
                                            plotCORMAT = TRUE,
                                            plotINTMISC = TRUE,
                                            plotPTMSTATS = TRUE,
                                            printPDF = TRUE,
                                            verbose = TRUE) {
  
  # Validate and prepare
  .artms_qc_validate_inputs(evidence_file, keys_file, prot_exp, output_dir)
  data_prepared <- .artms_qc_prepare_data(evidence_file, keys_file, prot_exp)
  
  # Create output directory
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Generate plots based on parameters
  if (plotINTDIST) {
    .artms_qc_plot_intensity_dist(data_prepared, output_dir, printPDF, verbose)
  }
  
  if (plotREPRO) {
    .artms_qc_plot_reproducibility(data_prepared, output_dir, printPDF, verbose)
  }
  
  if (plotCORMAT) {
    .artms_qc_plot_correlation(data_prepared, output_dir, printPDF, verbose)
  }
  
  if (plotINTMISC) {
    .artms_qc_plot_intensity_misc(data_prepared, output_dir, printPDF, verbose)
  }
  
  if (plotPTMSTATS & data_prepared$has_ptm) {
    .artms_qc_plot_ptm_stats(data_prepared, output_dir, printPDF, verbose)
  }
  
  if (verbose) message("QC analysis complete. Results in: ", output_dir)
  invisible(data_prepared)
}

# Helper functions (internal, not exported)
.artms_qc_validate_inputs <- function(evidence_file, keys_file, 
                                     prot_exp, output_dir) {
  # ~30 lines of validation
  ...
}

.artms_qc_prepare_data <- function(evidence_file, keys_file, prot_exp) {
  # ~80 lines of data preparation
  ...
}

.artms_qc_plot_intensity_dist <- function(data, output_dir, 
                                         printPDF, verbose) {
  # ~60 lines for this specific plot
  ...
}

# etc. for each plot type
```

**Benefits:**
- Easier to test individual components
- Easier to maintain
- Easier to reuse code
- More readable

---

### Task 7: Add Progress Indicators

**Priority:** P2  
**Estimated Time:** 3 days  

**Add to DESCRIPTION:**
```r
Imports:
  cli,
  # ... other imports
```

**Implementation Patterns:**

```r
#' @importFrom cli cli_alert_info cli_alert_success cli_progress_bar 
#' cli_progress_update cli_progress_done

artmsLongRunningFunction <- function(..., verbose = TRUE) {
  
  if (verbose) cli::cli_alert_info("Starting analysis...")
  
  # For operations with known steps
  n_steps <- 5
  if (verbose) {
    pb <- cli::cli_progress_bar(
      format = "Processing {cli::pb_current}/{cli::pb_total} steps",
      total = n_steps
    )
  }
  
  # Step 1
  ...
  if (verbose) cli::cli_progress_update()
  
  # Step 2
  ...
  if (verbose) cli::cli_progress_update()
  
  # ... continue for all steps
  
  if (verbose) {
    cli::cli_progress_done()
    cli::cli_alert_success("Analysis complete!")
  }
}
```

**For loops:**
```r
if (verbose) {
  cli::cli_alert_info("Processing {length(protein_list)} proteins...")
  pb <- cli::cli_progress_bar(total = length(protein_list))
}

for (protein in protein_list) {
  # Process protein
  ...
  
  if (verbose) cli::cli_progress_update()
}

if (verbose) cli::cli_progress_done()
```

---

## Testing Guidelines

### Unit Test Best Practices

```r
# Good test structure
test_that("function does X when Y", {
  # Arrange
  input_data <- setup_test_data()
  expected_result <- create_expected_output()
  
  # Act
  actual_result <- artms_function(input_data)
  
  # Assert
  expect_equal(actual_result, expected_result)
})

# Test edge cases
test_that("function handles empty input", {
  expect_error(artms_function(data.frame()))
})

test_that("function handles NA values", {
  data_with_na <- data.frame(x = c(1, NA, 3))
  result <- artms_function(data_with_na)
  expect_false(anyNA(result))
})
```

### Coverage Targets

```r
# Check coverage
covr::package_coverage()

# View coverage report
covr::report()

# Coverage by file
covr::file_coverage("R/qc-evidence-basic.R", "tests/testthat/test-qc-basic.R")
```

---

## Documentation Standards

### Roxygen2 Template

```r
#' Short one-line description
#'
#' Longer description with more details about what the function does.
#' Can span multiple paragraphs.
#'
#' @param param1 (type) Description of parameter 1
#' @param param2 (type) Description of parameter 2. Default is `default_value`
#' @param verbose (logical) If `TRUE` (default), shows function messages
#'
#' @return Description of return value. Be specific about the type and structure.
#'   If the function returns invisibly, note that.
#'
#' @family qc-functions
#' @seealso [artmsRelatedFunction()] for related functionality
#'
#' @examples
#' # Basic usage
#' result <- artmsFunction(
#'   param1 = "value1",
#'   param2 = 10
#' )
#'
#' # Advanced usage
#' result <- artmsFunction(
#'   param1 = "value1",
#'   param2 = 20,
#'   verbose = FALSE
#' )
#'
#' @export
artmsFunction <- function(param1, param2 = 10, verbose = TRUE) {
  ...
}
```

### Internal Function Documentation

```r
#' Short description of internal function
#'
#' More details about what this internal helper does.
#'
#' @param param1 Description
#' @return Description
#'
#' @keywords internal
#' @noRd
.artms_internal_helper <- function(param1) {
  ...
}
```

---

## Git Workflow

### Branch Strategy

```bash
# Create feature branch
git checkout -b feature/task-name

# Make changes
# ... edit files ...

# Commit frequently
git add R/file1.R tests/testthat/test-file1.R
git commit -m "Add input validation to function X"

# Run checks before pushing
R CMD check .
# or in R: devtools::check()

# Push when ready
git push origin feature/task-name

# Create PR on GitHub
```

### Commit Message Format

```
<type>: <short description>

<longer description if needed>

Fixes #issue_number
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `refactor`: Code refactoring
- `test`: Adding tests
- `docs`: Documentation updates
- `style`: Code style changes
- `perf`: Performance improvements

Examples:
```
feat: Add input validation with checkmate

- Add checkmate to DESCRIPTION
- Implement validation in all QC functions
- Add validator utility functions
- Update tests

Addresses #123
```

---

## Quality Checklist

Before submitting PR:

- [ ] Code follows tidyverse style guide
- [ ] All new functions have roxygen2 documentation
- [ ] All new functions have unit tests
- [ ] `devtools::check()` passes with no errors, warnings, or notes
- [ ] `devtools::test()` passes all tests
- [ ] Coverage maintained or improved (`covr::package_coverage()`)
- [ ] NEWS.md updated with changes
- [ ] No browser() or debug code left in
- [ ] All temporary files cleaned up

---

## Useful Commands

```r
# Development workflow
devtools::load_all()          # Load package for testing
devtools::document()          # Update documentation
devtools::test()              # Run tests
devtools::check()             # Full package check
devtools::install()           # Install package locally

# Styling
styler::style_pkg()           # Auto-format code

# Coverage
covr::package_coverage()      # Check coverage
covr::report()                # View coverage report

# Dependency checks
deps <- devtools::dev_package_deps()
update(deps)                  # Update dependencies
```

---

## Getting Help

- **Tidyverse style guide:** https://style.tidyverse.org/
- **R Packages book:** https://r-pkgs.org/
- **testthat documentation:** https://testthat.r-lib.org/
- **checkmate documentation:** https://mllg.github.io/checkmate/

---

*This guide will be updated as refactoring progresses.*
