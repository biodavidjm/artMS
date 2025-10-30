# artMS Refactoring Overview
## Quick Reference Summary

**Package:** artMS - Analytical R Tools for Mass Spectrometry  
**Current Version:** 1.10.3  
**Refactoring Goal:** Modernize R package while maintaining stability  

---

## 📊 Current State Snapshot

```
Package Metrics (v1.10.3):
├── Lines of Code:      ~14,000
├── R Source Files:     22
├── Exported Functions: 33
├── Test Files:         2
├── Test Coverage:      ~30%
├── Dependencies:       50 (39 Imports + 11 Suggests)
└── Bioconductor:       Compatible (v3.13+)

Functional Modules:
├── Quality Control (QC):           7 functions
├── MSstats Quantification:         8 functions
├── Downstream Analysis:            12 functions
└── Export Utilities:               6 functions
```

---

## 🎯 Refactoring Objectives Summary

| Category | Current Issues | Target Improvements |
|----------|----------------|---------------------|
| **Code Quality** | Mixed styles, some monolithic functions | Consistent tidyverse style, functions <200 lines |
| **Dependencies** | Deprecated packages (plyr, gProfileR) | Modern alternatives (dplyr, gprofiler2) |
| **Testing** | 30% coverage, 2 test files | 80% coverage, 11+ test files |
| **Documentation** | 1 vignette, basic examples | 6 vignettes, comprehensive examples |
| **Performance** | No parallelization | Optional parallel processing |
| **UX** | Basic error messages | Helpful errors, progress bars |

---

## 📅 Implementation Timeline

```
┌─────────────────────────────────────────────────────────────────┐
│ Phase 1: Foundation (Months 1-2) - CRITICAL                    │
├─────────────────────────────────────────────────────────────────┤
│ Week 1-2:   Replace plyr → dplyr                               │
│ Week 3-4:   Update gProfileR → gprofiler2                       │
│ Week 5-6:   Add input validation (checkmate)                   │
│ Week 7-8:   Expand test coverage to 60%                        │
│ Output:     v1.12.0 Release                                     │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│ Phase 2: Organization (Months 3-4) - ENHANCEMENTS              │
├─────────────────────────────────────────────────────────────────┤
│ Week 9-11:  Reorganize files, decompose functions              │
│ Week 12-14: Modernize plotting system                          │
│ Week 15-16: Add progress indicators, expand docs               │
│ Output:     v1.14.0 Release                                     │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│ Phase 3: Advanced (Months 5-6) - OPTIMIZATION                  │
├─────────────────────────────────────────────────────────────────┤
│ Week 17-19: Performance profiling & optimization               │
│ Week 20-22: Add S3 classes, parallel processing                │
│ Week 23-24: Final testing, documentation                       │
│ Output:     v1.16.0 Release                                     │
└─────────────────────────────────────────────────────────────────┘
```

---

## 🔑 Key Changes by Phase

### Phase 1: Critical Updates (v1.12.0)

**Dependencies:**
- ❌ Remove: `plyr` (superseded)
- ❌ Replace: `gProfileR` (deprecated) → ✅ `gprofiler2`
- ✅ Add: `checkmate` (input validation)

**Code Quality:**
- Standardize style with `styler::style_pkg()`
- Add comprehensive input validation
- Improve error messages

**Testing:**
- Create 8+ new test files
- Achieve 60%+ code coverage
- Add test fixtures

**Breaking Changes:** ⚠️ **NONE** (deprecation warnings only)

---

### Phase 2: Code Organization (v1.14.0)

**File Structure:**
```
Before (22 files):              After (26 files):
qualityControlEvidence*.R  →    qc-evidence-*.R (4 files)
evidenceToSAINT*.R        →    export-saint-*.R (2 files)
plots.R                   →    plot-*.R (4 files)
enrichments.R             →    analysis-enrichment.R
MSstats_functions.R       →    quantification-*.R (2 files)
+ utility files organized by purpose
```

**Code Improvements:**
- Break down functions >200 lines
- Extract plotting logic
- Create reusable helper functions
- ✅ Add: `cli` (progress bars)

**Documentation:**
- 6 comprehensive vignettes
- Enhanced function documentation
- Troubleshooting guide

**Breaking Changes:** ⚠️ **NONE** (internal refactoring only)

---

### Phase 3: Advanced Features (v1.16.0)

**Performance:**
- Profile and optimize bottlenecks
- Optional parallel processing
- Memory optimization for large datasets

**Modern R Features:**
- S3 classes for results
- Better pipe support
- Lifecycle badges (stable/experimental/deprecated)

**User Experience:**
- Consistent return structures
- Method chaining support
- Enhanced print/summary methods

**Breaking Changes:** ⚠️ Deprecated functions → **DEFUNCT** (throw errors)

---

## 📈 Success Metrics Dashboard

```
Code Quality Targets:
├── Linter Warnings:         0 (currently: ~50)
├── Function Complexity:     Max 200 lines (currently: max 800+)
├── Code Duplication:        <5% (currently: ~15%)
└── Style Consistency:       100% tidyverse

Testing Targets:
├── Phase 1:  60% coverage (currently: 30%)
├── Phase 2:  75% coverage
└── Phase 3:  80% coverage

Documentation Targets:
├── Vignettes:              6 (currently: 1)
├── Function Examples:      100% (currently: ~80%)
└── Help Pages Complete:    100%

Performance Targets:
├── QC Functions:           10-15% faster
├── Quantification:         5-10% faster
├── Downstream Analysis:    15-20% faster
└── Parallel Option:        2-4x speedup (on 4 cores)
```

---

## 🛡️ Backward Compatibility Strategy

### Deprecation Timeline

```
v1.10.3 (Current)
    ↓
v1.12.0 ──→ Warnings: "function X is deprecated, use Y instead"
    ↓       ├─ Old functions work normally
    ↓       └─ lifecycle::deprecate_warn() messages
    ↓
v1.14.0 ──→ Continued Warnings (same as v1.12.0)
    ↓       └─ More prominent documentation updates
    ↓
v1.16.0 ──→ Defunct: Deprecated functions throw errors
    ↓       ├─ lifecycle::deprecate_stop() errors
    ↓       └─ Clear migration path in error message
    ↓
v2.0.0  ──→ Removal: Deprecated functions completely removed
            └─ Major version bump allows breaking changes
```

### User Impact

| Version | User Action Required | Effort |
|---------|---------------------|---------|
| v1.12.0 | None (warnings only) | None |
| v1.14.0 | None (warnings only) | None |
| v1.16.0 | Update deprecated calls | 1-2 hours |
| v2.0.0 | Complete migration | 2-4 hours |

---

## 📚 Documentation Suite

Four comprehensive guides have been created:

### 1. [REFACTORING_PLAN.md](REFACTORING_PLAN.md) - Master Plan
**For:** Project managers, maintainers  
**Contains:**
- Detailed phase breakdown
- Task priorities and effort estimates
- File structure recommendations
- Risk assessment
- 20+ pages of comprehensive planning

### 2. [DEPENDENCIES.md](DEPENDENCIES.md) - Dependency Audit
**For:** Developers, system administrators  
**Contains:**
- All 50 dependencies analyzed
- Version compatibility matrix
- Migration paths for deprecated packages
- Health status and recommendations
- 12+ pages of dependency intelligence

### 3. [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) - User Guide
**For:** artMS users, analysts  
**Contains:**
- Version-by-version changes
- Before/after code examples
- Troubleshooting common issues
- Testing checklist
- 13+ pages of user-focused guidance

### 4. [IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md) - Developer Reference
**For:** Contributors, developers  
**Contains:**
- Task-by-task implementation steps
- Code templates and patterns
- Testing guidelines
- Git workflow standards
- 18+ pages of practical guidance

---

## 🔧 Quick Start for Implementers

### Phase 1 - First Steps

```r
# 1. Setup development environment
install.packages(c("devtools", "testthat", "checkmate", "styler", "lintr"))

# 2. Clone and setup
git clone https://github.com/biodavidjm/artMS.git
cd artMS
devtools::load_all()

# 3. Run existing tests
devtools::test()

# 4. Check current state
devtools::check()
covr::package_coverage()

# 5. Start with highest priority task
# See IMPLEMENTATION_GUIDE.md for detailed steps
```

### Priority Order

1. **Week 1-2:** Replace plyr with dplyr
   - File: `R/plots.R`, `R/enrichments.R`
   - See: IMPLEMENTATION_GUIDE.md § Task 1

2. **Week 3-4:** Update gProfileR to gprofiler2
   - File: `R/enrichments.R`
   - See: IMPLEMENTATION_GUIDE.md § Task 2

3. **Week 5-6:** Add input validation
   - All exported functions
   - See: IMPLEMENTATION_GUIDE.md § Task 3

4. **Week 7-8:** Expand test coverage
   - Create 8+ new test files
   - See: IMPLEMENTATION_GUIDE.md § Task 4

---

## 💡 Key Design Principles

### Maintained Throughout Refactoring

1. **Backward Compatibility First**
   - No breaking changes without major version bump
   - Clear deprecation warnings
   - Migration guides provided

2. **Incremental Improvement**
   - Small, tested changes
   - Each phase adds value
   - Can pause between phases

3. **User-Centric**
   - Better error messages
   - Progress indicators
   - Clear documentation

4. **Code Quality**
   - Consistent style
   - Well-tested
   - Easy to maintain

5. **Performance**
   - Don't make it slower
   - Optimize where beneficial
   - Optional parallelization

---

## 📊 Dependency Modernization Matrix

| Package | Status | Action | Phase | Priority |
|---------|--------|--------|-------|----------|
| plyr | ⚠️ Superseded | REMOVE → use dplyr | 1 | P0 |
| gProfileR | ❌ Deprecated | REPLACE → gprofiler2 | 1 | P0 |
| checkmate | ➕ New | ADD (validation) | 1 | P1 |
| cli | ➕ New | ADD (progress) | 2 | P2 |
| lifecycle | ➕ New | ADD (deprecation) | 2 | P2 |
| VennDiagram | ⚠️ Inactive | CONSIDER ggvenn | 3 | P3 |
| gplots | ⚠️ Maintenance | CONSIDER alternatives | 3 | P3 |

**Legend:**
- ❌ Broken/Deprecated (urgent)
- ⚠️ Superseded/Inactive (should update)
- ➕ Recommended addition
- ✅ Healthy (keep)

---

## 🎓 Learning Resources

### For Understanding artMS
- Website: http://artms.org
- Vignette: `vignette("artMS_vignette")`
- Paper: [Molecular & Cellular Proteomics](https://doi.org/10.1074/mcp.RA118.001218)

### For Refactoring Best Practices
- R Packages Book: https://r-pkgs.org/
- Tidyverse Style: https://style.tidyverse.org/
- Testing with testthat: https://testthat.r-lib.org/
- MSstats Documentation: http://msstats.org/

### For Specific Tasks
- checkmate: https://mllg.github.io/checkmate/
- cli: https://cli.r-lib.org/
- gprofiler2: https://cran.r-project.org/package=gprofiler2
- lifecycle: https://lifecycle.r-lib.org/

---

## 🚀 Expected Outcomes

### After Phase 1 (v1.12.0)
✅ No deprecated dependencies  
✅ Comprehensive input validation  
✅ 60% test coverage  
✅ Better error messages  
✅ Foundation for future improvements  

### After Phase 2 (v1.14.0)
✅ Well-organized codebase  
✅ Excellent documentation  
✅ 75% test coverage  
✅ Progress indicators  
✅ Easier maintenance  

### After Phase 3 (v1.16.0)
✅ Optimized performance  
✅ Modern R features  
✅ 80% test coverage  
✅ Optional parallelization  
✅ Production-ready for next 5+ years  

---

## 📞 Getting Help

### During Refactoring
- Issues: https://github.com/biodavidjm/artMS/issues
- Email: artms.help@gmail.com
- Discussions: GitHub Discussions (if enabled)

### Documentation
- This overview: Quick orientation
- REFACTORING_PLAN.md: Detailed planning
- IMPLEMENTATION_GUIDE.md: Step-by-step tasks
- MIGRATION_GUIDE.md: User transition help
- DEPENDENCIES.md: Package dependencies

---

## ✅ Pre-Implementation Checklist

Before starting refactoring:

- [ ] Read REFACTORING_PLAN.md fully
- [ ] Review DEPENDENCIES.md for context
- [ ] Scan IMPLEMENTATION_GUIDE.md for task details
- [ ] Setup development environment
- [ ] Run current tests to establish baseline
- [ ] Create feature branch for Phase 1
- [ ] Communicate plan to stakeholders
- [ ] Schedule regular progress reviews

---

## 📝 Notes

**Created:** 2025-10-30  
**Status:** Planning Complete, Ready for Implementation  
**Next Step:** Begin Phase 1, Task 1 (plyr → dplyr migration)  

**Estimated Total Effort:**
- Phase 1: 4-6 weeks (1-2 developers)
- Phase 2: 6-8 weeks (1-2 developers)  
- Phase 3: 4-6 weeks (1-2 developers)
- **Total:** ~6 months with 1-2 developers

**Risk Level:** Low-Medium (well-planned, backward compatible)

---

*For detailed information on any topic, see the respective guide documents.*
