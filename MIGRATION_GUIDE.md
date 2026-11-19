# artMS Migration Guide
## Transitioning to Modern artMS (v1.12.0+)

**Target Audience:** Existing artMS users upgrading from v1.10.x  
**Last Updated:** 2025-10-30  

---

## Overview

This guide helps existing artMS users transition to the modernized package versions. The refactoring maintains backward compatibility while introducing improvements in code quality, performance, and user experience.

---

## Version Timeline

### Current Stable: v1.10.3 (2021-07-14)
- Last pre-refactoring release
- All existing code continues to work
- MSstats 4.0+ compatible

### Planned Releases

**v1.12.0 (Phase 1 - Q2 2026)**
- Dependency updates (plyr → dplyr, gProfileR → gprofiler2)
- Enhanced input validation
- Improved error messages
- Expanded test coverage
- No breaking changes

**v1.14.0 (Phase 2 - Q4 2026)**
- Improved code organization
- Enhanced plotting capabilities
- Progress indicators
- Better documentation
- Deprecated functions marked with warnings

**v1.16.0 (Phase 3 - Q2 2027)**
- Performance optimizations
- Parallel processing support
- S3 classes for results
- Deprecated functions become defunct (error if called)

**v2.0.0 (Future - TBD)**
- Major version bump
- Complete removal of deprecated functions
- Potential API refinements

---

## What's Changing (and What's Not)

### ✅ Staying the Same

**Function names and signatures remain unchanged:**
- `artmsQualityControlEvidenceBasic()`
- `artmsQuantification()`
- `artmsAnalysisQuantifications()`
- All other exported functions maintain their names

**File formats remain compatible:**
- evidence.txt (MaxQuant output)
- keys.txt (experimental design)
- contrasts.txt (comparisons)
- config.yaml (configuration)

**Output formats remain consistent:**
- Same PDF plots
- Same result table structures
- Same file naming conventions

### 🔄 What's New/Improved

**Better error messages:**
```r
# Old (v1.10.x)
Error in artmsQualityControlEvidenceBasic(...) : 
  Missed (one or many) required argument(s)

# New (v1.12.0+)
Error in artmsQualityControlEvidenceBasic() :
  Assertion on 'evidence_file' failed: Must be either a string or a data.frame
  You provided: NULL
```

**Input validation:**
```r
# Now catches errors early
artmsQualityControlEvidenceBasic(
  evidence_file = "nonexistent.txt",  # Validates file exists
  prot_exp = "INVALID"                # Validates valid choice
)
# Error: File 'nonexistent.txt' not found
# Error: Must be element of set {'AB','PH','UB','AC','APMS','PTM:XXX:yy'}
```

**Progress indicators:**
```r
# New: Visual feedback for long operations
artmsQualityControlEvidenceBasic(...)
# ℹ Processing 15,234 features...
# ■■■■■■■■■■■■■■■■░░░░ 80% | ETA: 30s
```

---

## Dependency Changes

### Critical: gProfileR Replacement

**gProfileR is deprecated** and removed from CRAN. Enrichment functions now use **gprofiler2**.

#### Before (v1.10.x) - DEPRECATED
```r
# This will show warning in v1.12.0, error in v1.16.0+
library(gProfileR)
enrich <- artmsEnrichLog2fc(
  dataset = data_annotated,
  species = "human",
  background = background_genes
)
```

#### After (v1.12.0+) - RECOMMENDED
```r
# Install gprofiler2 (one time)
if (!requireNamespace("gprofiler2", quietly = TRUE))
    install.packages("gprofiler2")

# Same function call - gprofiler2 used automatically
enrich <- artmsEnrichLog2fc(
  dataset = data_annotated,
  species = "human",
  background = background_genes
)

# Or use new recommended function
enrich <- artmsEnrichProfiler(
  dataset = data_annotated,
  species = "hsapiens",  # Note: different format
  background = background_genes
)
```

**Migration steps:**
1. Install gprofiler2: `install.packages("gprofiler2")`
2. Update species codes if using `artmsEnrichProfiler()`:
   - "human" → "hsapiens"
   - "mouse" → "mmusculus"
3. Test enrichment analyses
4. Results structure unchanged

### Internal: plyr Removal

**You don't need to do anything** - this is internal only. But if you have custom code using artMS internals:

#### Before (v1.10.x)
```r
# If you were calling internal functions (not recommended)
result <- plyr::ddply(data, .(group), summarise, mean = mean(value))
```

#### After (v1.12.0+)
```r
# Use dplyr instead
result <- data %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(mean = mean(value))
```

---

## Function Deprecation Path

### Phase 1 (v1.12.0) - Warnings Only

Functions work normally but show warnings:

```r
artmsEnrichLog2fc(...)
# Warning: `artmsEnrichLog2fc()` is deprecated as of artMS 1.12.0
# ℹ Please use `artmsEnrichProfiler()` instead with gprofiler2
# This warning is displayed once per session.
```

**Action required:** None immediately, but plan to update code

### Phase 2 (v1.14.0) - Continued Warnings

Same behavior as v1.12.0, but more prominent warnings.

**Action required:** Update code to new functions

### Phase 3 (v1.16.0) - Defunct

Deprecated functions throw errors:

```r
artmsEnrichLog2fc(...)
# Error: `artmsEnrichLog2fc()` is defunct as of artMS 1.16.0
# Please use `artmsEnrichProfiler()` instead
# See ?artmsEnrichProfiler for details
```

**Action required:** Must update code before upgrading to v1.16.0

---

## Updating Your Analysis Scripts

### Typical artMS Workflow

#### Old Script (v1.10.x) - Still works in v1.12.0+
```r
library(artMS)

# Quality Control
artmsQualityControlEvidenceBasic(
  evidence_file = "evidence.txt",
  keys_file = "keys.txt",
  prot_exp = "PH"
)

# Quantification
artmsQuantification(
  yaml_config_file = "config.yaml"
)

# Downstream analysis
artmsAnalysisQuantifications(
  log2fc_file = "results-log2fc.txt",
  species = "human",
  output_dir = "downstream"
)
```

#### Updated Script (v1.12.0+) - Recommended
```r
library(artMS)

# Quality Control - now with progress indicators
artmsQualityControlEvidenceBasic(
  evidence_file = "evidence.txt",
  keys_file = "keys.txt",
  prot_exp = "PH",
  verbose = TRUE  # See detailed progress
)

# Quantification - unchanged
artmsQuantification(
  yaml_config_file = "config.yaml"
)

# Downstream analysis - ensure gprofiler2 installed
if (!requireNamespace("gprofiler2", quietly = TRUE)) {
  message("Installing gprofiler2 for enrichment analysis...")
  install.packages("gprofiler2")
}

artmsAnalysisQuantifications(
  log2fc_file = "results-log2fc.txt",
  species = "human",
  output_dir = "downstream"
)
```

### Error Handling

#### Old (v1.10.x)
```r
# Errors were sometimes cryptic
result <- artmsQuantification(yaml_config_file = "config.yaml")
# Error in something complicated...
```

#### New (v1.12.0+)
```r
# Errors are more helpful
result <- artmsQuantification(yaml_config_file = "config.yaml")
# Error in artmsQuantification():
#   Configuration file 'config.yaml' not found
#   Expected location: /path/to/working/dir/config.yaml
#   Create one with: artmsWriteConfigYamlFile()
```

---

## Configuration File Updates

### No Changes Required

Your existing config.yaml files continue to work without modification:

```yaml
# This continues to work in all versions
files:
  evidence: evidence.txt
  keys: keys.txt
  contrasts: contrasts.txt
  output: results.txt

msstats:
  normalization_method: equalizeMedians
  summaryMethod: TMP
  # ... all other parameters unchanged
```

### New Optional Parameters (v1.12.0+)

```yaml
# Optional: new parameters for enhanced features
output_options:
  verbose: yes           # More detailed messages
  show_progress: yes     # Progress bars
  
performance:
  parallel: no           # Enable parallel processing (Phase 3)
  n_cores: 4            # Number of cores if parallel enabled
```

---

## Testing Your Migration

### Step-by-Step Testing

1. **Install the new version**
```r
# From Bioconductor
BiocManager::install("artMS")

# Check version
packageVersion("artMS")
# Should be >= 1.12.0
```

2. **Test with example data**
```r
# Use built-in test data first
artmsQualityControlEvidenceBasic(
  evidence_file = artms_data_ph_evidence,
  keys_file = artms_data_ph_keys,
  prot_exp = "PH",
  printPDF = FALSE,  # Don't create PDFs for testing
  verbose = TRUE
)
```

3. **Test with your data**
```r
# Run your existing scripts
source("my_artms_analysis.R")

# Check for any warnings
# Update code if deprecation warnings appear
```

4. **Verify outputs**
```r
# Compare results between versions
# Results should be identical (within numerical precision)
```

---

## Common Migration Issues

### Issue 1: gProfileR Not Found

**Symptom:**
```r
Error: Package 'gProfileR' required but not installed
```

**Solution:**
```r
# Install gprofiler2 instead
install.packages("gprofiler2")

# Update your enrichment calls
artmsEnrichProfiler(...)  # Use this instead of artmsEnrichLog2fc()
```

### Issue 2: Different Enrichment Results

**Symptom:**
Results from gprofiler2 differ slightly from gProfileR

**Explanation:**
- gprofiler2 uses updated databases
- Slightly different algorithm improvements
- This is expected and results are more current

**Solution:**
No action needed - new results are more accurate

### Issue 3: Progress Messages Are Verbose

**Symptom:**
Too many messages during processing

**Solution:**
```r
# Reduce verbosity
artmsQualityControlEvidenceBasic(
  ...,
  verbose = FALSE
)
```

### Issue 4: Function Not Found After Update

**Symptom:**
```r
Error: could not find function "artms_internal_function"
```

**Cause:**
You're calling internal functions (not recommended)

**Solution:**
- Use only exported functions (those starting with `artms`, not `.artms_`)
- If you need functionality from internal functions, submit an issue requesting export

---

## Performance Improvements

### v1.12.0+: Optimized Data Processing

Your analyses may run faster due to:
- Better data.table usage
- Removed redundant operations
- Optimized plotting code

**Typical improvements:**
- QC functions: 10-15% faster
- Quantification: 5-10% faster
- Downstream analysis: 15-20% faster

### v1.16.0+ (Phase 3): Parallel Processing

```r
# Enable parallel processing for faster analysis
artmsQuantification(
  yaml_config_file = "config.yaml",
  parallel = TRUE,
  n_cores = 4  # Adjust based on your system
)
```

---

## Getting Help

### If Something Doesn't Work

1. **Check the documentation**
```r
?artmsQualityControlEvidenceBasic
?artmsQuantification
```

2. **Review the NEWS file**
```r
news(package = "artMS")
```

3. **Check existing issues**
   - Visit: https://github.com/biodavidjm/artMS/issues
   - Search for your problem

4. **Submit a new issue**
   - Include minimal reproducible example
   - Include sessionInfo() output
   - Specify version numbers

5. **Email support**
   - artms.help@gmail.com
   - Include version info and error messages

---

## Version Comparison Quick Reference

| Feature | v1.10.3 | v1.12.0 | v1.14.0 | v1.16.0 |
|---------|---------|---------|---------|---------|
| Basic functions | ✅ | ✅ | ✅ | ✅ |
| gProfileR | ✅ | ⚠️ | ⚠️ | ❌ |
| gprofiler2 | ➖ | ✅ | ✅ | ✅ |
| Input validation | Basic | Enhanced | Enhanced | Enhanced |
| Error messages | Basic | Improved | Improved | Improved |
| Progress bars | ❌ | ✅ | ✅ | ✅ |
| Test coverage | 30% | 60% | 75% | 80% |
| Parallel processing | ❌ | ❌ | ❌ | ✅ |

**Legend:**
- ✅ Available and supported
- ⚠️ Deprecated (warning shown)
- ❌ Removed/Not available
- ➖ Not yet available

---

## Rollback Instructions

If you need to return to v1.10.3:

```r
# Remove current version
remove.packages("artMS")

# Install specific old version (not recommended)
# Note: May have dependency conflicts
library(devtools)
install_github("biodavidjm/artMS@v1.10.3")
```

**Warning:** Rolling back is not recommended. Instead:
1. Report issues to GitHub
2. Use verbose mode to debug problems
3. Contact support for assistance

---

## Summary Checklist

Before upgrading to v1.12.0+:

- [ ] Backup your analysis scripts
- [ ] Read the NEWS file for your version
- [ ] Install gprofiler2 if using enrichment functions
- [ ] Test with example data first
- [ ] Run your analyses with `verbose = TRUE` to see any warnings
- [ ] Update scripts to remove any deprecation warnings
- [ ] Document any issues and report them

After upgrading:

- [ ] Verify outputs match expectations
- [ ] Update documentation/notes about version used
- [ ] Enjoy improved error messages and progress indicators!
- [ ] Consider contributing feedback or improvements

---

## Additional Resources

- **Package website:** http://artms.org
- **GitHub repository:** https://github.com/biodavidjm/artMS
- **Bioconductor page:** https://bioconductor.org/packages/artMS
- **Issue tracker:** https://github.com/biodavidjm/artMS/issues
- **Email support:** artms.help@gmail.com

---

*This migration guide will be updated as new versions are released. Last update: 2025-10-30*
