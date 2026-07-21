# DeFrABB Release Process

## Overview

This document describes the process for releasing a new version of DeFrABB. The process ensures that releases are stable, well-tested, and documented before being pushed to the public GitHub repository.

## Release Stages

### 1. Code Freeze & Stabilization

**Goal:** Ensure master is stable and all planned features are merged.

**Steps:**
1. Complete all planned features and bug fixes for the milestone
2. Merge all validated experiment branches
3. Update CHANGELOG with all changes in the "Unreleased" section
4. Ensure all tests pass: `pytest .tests`
5. Ensure code formatting is clean: `snakefmt .` and `black scripts/`
6. Verify CI is green on GitLab

**Checklist:**
- [ ] All milestone features merged
- [ ] All experiment branches reviewed (merge or defer)
- [ ] CHANGELOG updated
- [ ] Tests passing locally
- [ ] CI green on GitLab
- [ ] No uncommitted changes

---

### 2. Full Pipeline Validation Test

**Goal:** Run comprehensive whole-genome regression test that CI cannot cover.

**Why:** CI only runs chr21 test data due to time/compute constraints. The full pipeline test exercises:
- All three references (GRCh38, GRCh37, CHM13v2.0)
- Both variant callers (dipcall, PAV)
- All variant types (smvar, stvar)
- All major features (exclusions, stratifications, annotations, evaluations)

**Steps:**

1. **Commit all changes** to master before running the test:
   ```bash
   git status  # verify clean
   git log -1  # verify latest commit
   ```

2. **Launch the release test:**
   ```bash
   conda activate snakemake
   scripts/release_test.sh --cores 88 --detach
   ```
   
   This will:
   - Clone repo to `../runs/YYYYMMDD_v0.0XX_fulltest-HG002v1.1/`
   - Generate analyses table with fulltest config
   - Validate and dry-run
   - Launch in tmux session `defrabb-YYYYMMDD_v0.0XX_fulltest-HG002v1.1`
   - Run pipeline in background (detached)
   - Auto-generate comparison results when complete

3. **Monitor progress:**
   ```bash
   tmux attach -t defrabb-<RUNID>  # attach to session
   # Ctrl-b d to detach
   
   # Or check log directly:
   tail -f ../runs/<RUNID>/run.log
   ```

4. **Wait for completion** (typically 8-14 hours for full genome)

5. **Review results:**
   ```bash
   cd ../runs/<RUNID>/
   cat comparison_<RUNID>.tsv  # automated comparison
   ls -lh results/             # all outputs
   ```

**Pass criteria:**
- Pipeline completes without errors
- All expected outputs generated
- Comparison metrics look reasonable (no major regressions vs v5q baseline)
- Known issues (from docs/issues/) accounted for

**Fail actions:**
- Fix bugs and restart from step 1
- Document new issues in docs/issues/
- DO NOT proceed to release until test passes

---

### 3. Prepare Release

**Goal:** Finalize documentation and tag the release.

**Steps:**

1. **Update CHANGELOG:**
   - Move "Unreleased" section to "Milestone v0.0XX YYYY-MM-DD"
   - Add summary of major changes
   - Reference merged issues/branches
   - Keep release focus description (e.g., "stabilization", "parameter optimization")

2. **Update version references:**
   - Check CLAUDE.md for version mentions
   - Update run_full_pipeline_test.sh default version if needed
   - Update release_test.sh default version

3. **Create release tag:**
   ```bash
   git tag -a v0.0XX -m "Release v0.0XX: <summary>

   Major changes:
   - <feature 1>
   - <feature 2>
   - <bug fix 1>

   See CHANGELOG for full details."
   
   git push origin v0.0XX
   git push origin master
   ```

4. **Create GitHub release** (if pushing to GitHub):
   - Use tag v0.0XX
   - Copy CHANGELOG milestone section as release notes
   - Attach fulltest comparison results (optional)

---

### 4. Public Release (GitHub)

**Goal:** Make release available on public GitHub repository.

**Prerequisites:**
- GitLab tag v0.0XX created and pushed
- Full pipeline test passed
- CHANGELOG finalized
- All documentation up-to-date

**Steps:**

1. **Push to GitHub:**
   ```bash
   # Add GitHub remote if not already added
   git remote add github https://github.com/usnistgov/defrabb.git
   
   # Push master and tags
   git push github master
   git push github --tags
   ```

2. **Create GitHub release:**
   - Go to https://github.com/usnistgov/defrabb/releases
   - Click "Draft a new release"
   - Select tag v0.0XX
   - Title: "v0.0XX - <focus description>"
   - Description: Copy from CHANGELOG milestone section
   - Click "Publish release"

3. **Verify GitHub state:**
   - Check that master branch shows latest commit
   - Check that tags are visible
   - Check that release appears correctly
   - Verify README displays properly

---

## Release Checklist Template

Copy this for each release:

```markdown
## Release v0.0XX Checklist

### Code Freeze
- [ ] All milestone features merged
- [ ] Experiment branches reviewed (merged/deferred)
- [ ] CHANGELOG updated (Unreleased → Milestone)
- [ ] Tests passing: `pytest .tests`
- [ ] Formatting clean: `snakefmt .` && `black scripts/`
- [ ] CI green on GitLab
- [ ] No uncommitted changes

### Full Pipeline Test
- [ ] Changes committed to master
- [ ] Launched: `scripts/release_test.sh --cores 88 --detach`
- [ ] Session: `defrabb-<RUNID>`
- [ ] Completed successfully (no errors)
- [ ] Results reviewed: `comparison_<RUNID>.tsv`
- [ ] Metrics reasonable (no major regressions)

### Prepare Release
- [ ] CHANGELOG finalized (Unreleased → Milestone v0.0XX)
- [ ] Version references updated
- [ ] Tag created: `git tag -a v0.0XX -m "..."`
- [ ] Tag pushed: `git push origin v0.0XX`
- [ ] Master pushed: `git push origin master`

### Public Release (GitHub)
- [ ] Pushed to GitHub: `git push github master`
- [ ] Tags pushed: `git push github --tags`
- [ ] GitHub release created with CHANGELOG notes
- [ ] GitHub state verified (master, tags, release)
```

---

## Troubleshooting

### Full Pipeline Test Fails

**Common issues:**

1. **Download failures:**
   - Retry with: `cd <RUNID> && ./run_defrabb run -r <RUNID> -j 88`
   - Check network connectivity
   - Verify conda mirror availability

2. **OOM errors:**
   - Check memory allocations in config/resources.yml
   - Reduce --cores if running on constrained system
   - Review logs for specific rule failures

3. **Rule failures:**
   - Check rule-specific logs in `<RUNID>/logs/`
   - Known issues documented in `docs/issues/`
   - File new issue if unknown problem

### CI/CD vs Full Test

**Why not run full test in CI?**
- Time: 8-14 hours vs 15-minute CI limit
- Compute: 100+ GB memory vs limited CI resources
- Cost: Full genome datasets vs chr21 test data

**CI covers:**
- Code formatting (snakefmt, black)
- Linting (snakemake --lint)
- Unit tests (pytest)
- Dry-run validation on chr21 test data

**Full test covers:**
- All references (GRCh38, GRCh37, CHM13v2.0)
- Both variant callers (dipcall, PAV)
- Large-scale exclusion/stratification/annotation
- Full evaluation against v5q baselines

Both are necessary: CI for fast feedback on every commit, full test for release validation.

---

## Post-Release

1. **Create next milestone** in GitLab
2. **Update development roadmap** if needed
3. **File issues** for next milestone tasks
4. **Archive release test clone** (optional):
   ```bash
   tar -czf <RUNID>.tar.gz <RUNID>/results/ <RUNID>/comparison_*.tsv
   # Move to archive location
   ```

---

## Release Schedule

**Current practice:** Release when milestone objectives met, typically after:
- Major feature completion
- Bug fix accumulation
- Stabilization period

**Target cadence:** Flexible, based on development activity. Typically 2-4 months between releases.

---

## Version Numbering

**Format:** `v0.0XX`

**Scheme:**
- Major version: Always 0 (pre-1.0 development)
- Minor version: Increments with each milestone (021, 022, 023, ...)
- No patch versions currently (minor version increments cover all changes)

**Example progression:**
- v0.021 - PAV integration, HG008 support
- v0.022 - Stabilization, FIPS compatibility, bug fixes
- v0.023 - Parameter optimization infrastructure (planned)

---

## Contact

For questions about the release process:
- Internal: Contact pipeline maintainer
- External: File issue on GitHub repository
