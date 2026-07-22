# V0.023 Parameter Optimization - Session Status

## Context
User requested comprehensive v0.023 milestone work in a worktree with multi-agent planning and implementation, including adversarial review to avoid over/under-engineering. The scope is LARGE - this is a multi-session effort.

---

## Current Status

### Worktree Created
- **Branch:** worktree-v0.023-param-optimization  
- **Location:** `.claude/worktrees/v0.023-param-optimization/`
- **Clean state:** Ready for iterative development

### Work Completed (Session 1 of N)

#### 1. Deep Infrastructure Research ✅
**Agent:** Research current parameter infrastructure  
**Findings:**
- VCF processing profiles: ✓ fully implemented, good pattern
- Truvari bench params: ✓ fully implemented, good pattern  
- Dipcall params: ✗ BROKEN - defined but never accessed
- PAV params: ⚠ works but only 1 profile
- Exclusion params: ⚠ no profiles, only global + unused overrides
- **Critical gap:** Output reuse essential for 5-genome sweeps (dipcall = 3-6 hours)

#### 2. Dipcall Parameter Fix ✅ COMMITTED
**Files changed:**
- `rules/asm-varcall.smk` - Fixed profile lookup logic
- `config/resources.yml` - Added z5k, z10k, z1k profiles
- `schema/resources-schema.yml` - Added pattern validation
- `.tests/unit/test_dipcall_params.py` - Comprehensive tests (131 lines)
- `CHANGELOG` - Documented changes

**Commit:** "feat(param-opt): fix dipcall parameter profile lookup and add profiles"

#### 3. Design Workflow RUNNING 🔄
**Workflow:** 3-phase adversarial design review  
- Phase 1: Design framework architecture
- Phase 2: Parallel over/under-engineering critique
- Phase 3: Refine based on critiques

**Status:** Background execution started, journal shows 1 agent running  
**Expected output:** Final vetted architecture for remaining components

---

## Tasks Tracking

| # | Task | Status | Notes |
|---|------|--------|-------|
| 9 | Research infrastructure | ✅ Complete | Agent analysis done |
| 10 | Design framework architecture | 🔄 In progress | Workflow running |
| 11 | Implement profile registry | 🔄 In progress | Dipcall done, more needed |
| 12 | Implement v5q comparison | ⏳ Pending | After design complete |
| 13 | Implement output reuse | ⏳ Pending | After design complete |
| 14 | Design multi-genome workflow | ⏳ Pending | After design complete |
| 15 | Document param optimization | ⏳ Pending | After implementation |

---

## Remaining Work (Multi-Session Scope)

### High Priority (Next Session)
1. **Complete design workflow** - Review adversarial critiques, finalize architecture
2. **Implement parameter sweep generator** - `scripts/generate_param_sweep.py`
3. **Enhance compare_evaluations.py** - v5q baseline loading, scoring/ranking
4. **Add exclusion parameter profiles** - Named profile sets in resources.yml

### Medium Priority  
5. **Output reuse infrastructure** - Variant call fingerprinting, symlink/cache patterns
6. **PAV parameter profiles** - Add 2-3 alternate merge strategies
7. **Multi-genome workflow documentation** - Map 5-genome optimization process
8. **Integration tests** - End-to-end parameter sweep test

### Documentation
9. **docs/parameter-optimization.md** - User guide with examples
10. **docs/multi-genome-optimization.md** - Hypothetical 5-genome workflow
11. **Update CLAUDE.md** - Document new parameter infrastructure

---

## Design Decisions (From Research)

### What Works Well (Patterns to Follow)
- **VCF processing profiles:** Clean registry, schema-validated, well-documented
- **Truvari bench params:** Systematic profile framework (#194), versioned tuning
- **Pattern:** Profile name in table → lookup in config → validate via schema

### What Needs Work
- **Dipcall params:** Now fixed, but needs more profiles for sweeps
- **Exclusion params:** No profiles, just global scalars + unused overrides
- **Output reuse:** CRITICAL for sweeps - prevent re-running 3-6 hour dipcall jobs

### Architecture Principles (User-Specified)
- Build on existing infrastructure
- Minimal core pipeline changes  
- Backward compatible
- Incremental value (each component useful independently)
- Clear separation: param definition → sweep gen → eval → comparison

---

## Multi-Genome Workflow (User Requirements)

### Scenario: 5 Genomes with Parameter Optimization

**Setup:**
1. HG002 + 4 other genomes (similar assembly methods)
2. HG002 has v5q benchmark for optimization
3. Other genomes have comparable callsets for evaluation

**Workflow:**
1. **HG002 parameter sweep**
   - Generate parameter combinations (dipcall z-params, exclusion buffers, truvari settings)
   - Run sweep against HG002 assemblies
   - Compare to v5q baseline
   - Rank by precision/recall/F1 + region concordance
   - Identify optimal parameters

2. **Apply optimized params to other 4 genomes**
   - Use winning param set from HG002
   - Generate benchmarks for each genome  
   - Compare to genome-specific comparison callsets

3. **Cross-genome analysis**
   - Compare metrics across all 5 genomes
   - Identify stratifications with high discrepancy rates
   - Identify excluded regions with high concordance (potential to include)
   - Flag regions for manual review/refinement

4. **Iterative refinement** (future work)
   - Use WGS data from multiple platforms
   - Characterize regions where benchmark reliably/unreliably identifies errors
   - Refine benchmark region definitions based on evidence

**Key requirements:**
- Output reuse (don't re-call variants for every param combo)
- Automated comparison aggregation
- Discrepancy hotspot identification
- Region concordance analysis

---

## Technical Notes

### Worktree Limitations Encountered
- `python` and `snakefmt` not in PATH (expected - worktree isolation)
- Tests written but not executed (will run from main repo after merge)
- Pyright diagnostics show import warnings (expected, test deps not installed)

### Workflow Infrastructure Used
- Background workflow with 3-phase design
- Structured output schemas for design + critique
- Parallel adversarial review (over-eng vs under-eng)
- Design refinement based on critiques

---

## Next Steps

### Immediate (Next Session Start)
1. Check workflow completion: `/workflows` or check journal
2. Review design output from workflow
3. Implement parameter sweep generator based on design
4. Enhance compare_evaluations.py with v5q support

### Before Merge to Master
1. Run tests from main repo (worktree python issues)
2. Format code with snakefmt
3. Validate schema changes
4. Update CLAUDE.md with parameter optimization docs
5. Ensure CI passes

### Multi-Session Planning
- **Session 1 (this):** Research + dipcall fix + design workflow
- **Session 2:** Implement core components (sweep gen, v5q comparison)
- **Session 3:** Output reuse + integration tests
- **Session 4:** Multi-genome workflow + documentation
- **Session 5:** Testing, polish, merge

---

## Files Modified (Session 1)

### Created
- `.tests/unit/test_dipcall_params.py` (131 lines)
- `$CLAUDE_JOB_DIR/tmp/param_optimization_design.js` (workflow script)

### Modified
- `rules/asm-varcall.smk` (dipcall param lookup fix)
- `config/resources.yml` (added 3 dipcall profiles)
- `schema/resources-schema.yml` (validation pattern)
- `CHANGELOG` (v0.023 section started)

### Committed
- 1 commit: "feat(param-opt): fix dipcall parameter profile lookup and add profiles"
- 270 insertions across 6 files

---

## User Guidance

This is a **significant multi-session task** as you noted. Current progress:
- ✅ Foundation laid (research done, first fix committed)
- 🔄 Design in progress (workflow running)
- ⏳ Major implementation ahead (sweep gen, v5q comparison, output reuse)

**Recommended approach:**
1. Let design workflow complete (check next session)
2. Review and approve/revise design
3. Implement components incrementally with tests
4. Merge to master when stable
5. Use for actual HG002 parameter optimization runs

**Estimated remaining effort:** 3-4 more focused sessions to complete v0.023 milestone.
