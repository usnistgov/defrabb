# plan: move `giab-test-data/` to a dedicated public bucket

> **Status: completed and archived (2026-04-16).** The migration described
> here has been executed; test data now lives in the dedicated
> `giab-test-data` bucket and `config/resources.yml` has been updated to
> match. The `manage_giab_test_data_s3.py` helper script referenced
> throughout this document was a one-off admin tool used during the
> migration and is **not** maintained inside this repository — it lives
> outside the codebase as personal tooling. This document is retained
> here for historical reference only.

This document outlines a post-restore migration plan for moving `s3://giab-data/giab-test-data/` into a separate bucket so the data remains directly downloadable over HTTPS and is no longer affected by the bucket-wide `archiving` Intelligent-Tiering rule in `giab-data`.

## goal

Create a new bucket dedicated to `giab-test-data` that:

- keeps objects directly accessible over HTTPS
- does **not** transition data into `ARCHIVE_ACCESS` or `DEEP_ARCHIVE_ACCESS`
- is easier to manage independently from the rest of `giab-data`
- supports current DeFrABB resource URLs after a controlled config update

## why move to a new bucket

Current state:

- `giab-data` has a bucket-wide Intelligent-Tiering configuration:
  - `ARCHIVE_ACCESS` after 90 days
  - `DEEP_ARCHIVE_ACCESS` after 180 days
- `giab-test-data/` objects in `INTELLIGENT_TIERING` are therefore eligible for archive transitions
- once archived, public HTTPS reads fail with `403 InvalidObjectState`

Moving `giab-test-data/` to a separate bucket avoids that coupling and lets the team apply a storage policy specifically for always-downloadable public test data.

## recommended target design

### bucket

Suggested bucket name examples:

- `giab-test-data-public`
- `giab-public-test-data`
- `giab-downloads`

Pick a name that is stable and clearly indicates that the bucket contains public, directly downloadable data.

### storage class policy

Recommended options:

- **preferred:** `STANDARD`
- **acceptable:** `INTELLIGENT_TIERING` **without** optional archive tiers

Avoid:

- bucket rules that enable `ARCHIVE_ACCESS`
- bucket rules that enable `DEEP_ARCHIVE_ACCESS`
- lifecycle transitions to Glacier classes for this bucket

### HTTPS access

Use standard S3 HTTPS URLs:

```text
https://<new-bucket>.s3.amazonaws.com/giab-test-data/<path>
```

That keeps the access pattern similar to the current URLs and works well for `curl`, Snakemake rules, and browser access.

### public access model

If the current intent is public anonymous download, configure:

- bucket policy allowing public `s3:GetObject`
- only for the object paths that should be public
- no public write access

If public access is not intended for every object, restrict policy accordingly.

## migration prerequisites

Before starting:

1. wait for currently archived objects in `giab-data/giab-test-data/` to finish restoring
2. confirm representative files are readable again over HTTPS
3. inventory current consumers of the old URLs
4. confirm who will own the new bucket and apply its policy/configuration

## phased migration plan

## phase 1: create and configure the new bucket

### 1. create the bucket

Example:

```bash
aws s3api create-bucket \
  --bucket giab-test-data-public \
  --region us-east-1
```

If the bucket is in another region, include the appropriate location configuration.

### 2. configure public-read access for objects

Use a bucket policy similar to:

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "PublicReadForGiabTestData",
      "Effect": "Allow",
      "Principal": "*",
      "Action": "s3:GetObject",
      "Resource": "arn:aws:s3:::giab-test-data-public/giab-test-data/*"
    }
  ]
}
```

Apply it with:

```bash
aws s3api put-bucket-policy \
  --bucket giab-test-data-public \
  --policy file://bucket-policy.json
```

### 3. keep storage directly readable

Preferred:

- upload/copy objects as `STANDARD`

Alternative:

- use `INTELLIGENT_TIERING`
- but do **not** create any Intelligent-Tiering archive configuration for this bucket

### 4. explicitly avoid archive lifecycle/config rules

Before copying data, verify:

- no lifecycle transitions to Glacier/Deep Archive
- no Intelligent-Tiering archive/deep-archive configuration

Useful checks:

```bash
aws s3api get-bucket-lifecycle-configuration \
  --bucket giab-test-data-public
```

```bash
aws s3api list-bucket-intelligent-tiering-configurations \
  --bucket giab-test-data-public
```

## phase 2: validate restored source objects

Before copying all data, spot-check a few restored source files:

```bash
python scripts/manage_giab_test_data_s3.py audit \
  --all-prefix-objects \
  --bucket giab-data \
  --prefix giab-test-data/ \
  --check-http
```

Check especially:

- `GRCh38_chr21.fa.gz`
- `CHM13v2.0/reference/chm13v2.0_maskedY_toupper.fa.gz`
- `GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz`

## phase 3: copy data to the new bucket

### 1. copy the prefix

Example:

```bash
aws s3 sync \
  s3://giab-data/giab-test-data/ \
  s3://giab-test-data-public/giab-test-data/
```

If you want all copied objects to land in a specific class:

```bash
aws s3 sync \
  s3://giab-data/giab-test-data/ \
  s3://giab-test-data-public/giab-test-data/ \
  --storage-class STANDARD
```

### 2. preserve metadata where needed

After copying, spot-check:

- content types
- object sizes
- checksums/ETags where meaningful
- compression suffixes and filenames

### 3. optional: perform incremental syncs

If the source prefix is still changing, do:

1. initial full sync
2. validation
3. final incremental sync just before cutover

## phase 4: validate the new bucket

### 1. metadata checks

Spot-check several objects:

```bash
aws s3api head-object \
  --bucket giab-test-data-public \
  --key giab-test-data/GRCh38_chr21.fa.gz
```

Confirm:

- expected object exists
- correct size
- readable storage class
- no archive status

### 2. HTTPS download checks

Spot-check direct reads:

```bash
curl -I https://giab-test-data-public.s3.amazonaws.com/giab-test-data/GRCh38_chr21.fa.gz
```

```bash
curl -L -r 0-1023 \
  https://giab-test-data-public.s3.amazonaws.com/giab-test-data/GRCh38_chr21.fa.gz \
  >/dev/null
```

Expected result:

- HTTP success
- no `InvalidObjectState`

### 3. script-based validation

The current helper script can validate the new bucket too:

```bash
python scripts/manage_giab_test_data_s3.py audit \
  --all-prefix-objects \
  --bucket giab-test-data-public \
  --prefix giab-test-data/ \
  --check-http
```

## phase 5: update DeFrABB configuration

Update `config/resources.yml` to replace URLs of the form:

```text
https://giab-data.s3.amazonaws.com/giab-test-data/...
```

with:

```text
https://giab-test-data-public.s3.amazonaws.com/giab-test-data/...
```

After updating:

1. run a targeted grep to verify replacements
2. run a focused workflow smoke test for rules that download references/resources

Suggested quick checks:

```bash
rg -n "giab-data.s3.amazonaws.com/giab-test-data|giab-test-data-public.s3.amazonaws.com" config/resources.yml
```

Then test at least one workflow path that downloads:

- a reference FASTA
- an exclusion BED
- a comparison VCF

## phase 6: cutover and transition

Options:

### option A: immediate cutover

- update configs in DeFrABB and other known consumers
- publish the new base URL
- treat old URLs as deprecated

### option B: transition period

- keep the old data in `giab-data/giab-test-data/` for a limited time
- update active projects first
- remove or freeze the old prefix later

Because standard S3 bucket endpoints do not provide a simple transparent redirect strategy for all existing object URLs, assume most downstream users will need to update their URLs.

## phase 7: long-term cleanup

After cutover:

1. document the new bucket as the canonical home for public GIAB test data
2. make sure no archive-tier policies are later added to the new bucket
3. optionally remove or stop updating the old prefix
4. add periodic audit checks for accessibility

## validation checklist

Before calling the migration done, confirm:

- [ ] all needed source objects finished restore
- [ ] new bucket exists
- [ ] new bucket policy allows intended HTTPS reads
- [ ] new bucket has no archive/deep-archive Intelligent-Tiering config
- [ ] new bucket has no lifecycle transitions to Glacier classes
- [ ] full copy completed
- [ ] representative HTTPS URLs work
- [ ] `config/resources.yml` updated
- [ ] DeFrABB smoke test passes
- [ ] downstream users notified of URL change

## rollback plan

If the cutover exposes a problem:

1. keep using restored objects from `giab-data/giab-test-data/`
2. revert `config/resources.yml` to the old URLs
3. fix policy/storage issues in the new bucket
4. repeat validation before another cutover

Because the migration is copy-based rather than move/delete-first, rollback should be low risk as long as the source objects remain available during the transition.

## suggested follow-up improvements

After the migration, consider:

- renaming `manage_giab_test_data_s3.py` to a more general S3 archive-management name
- extending the script with a copy-validation mode
- adding a machine-readable audit summary for bulk migration checks
- documenting the new bucket in project README or operational docs

## rough execution order

1. wait for restore completion
2. verify source HTTPS access
3. create/configure new bucket
4. verify no archive-tier rules on the new bucket
5. sync data
6. validate new HTTPS URLs
7. update `config/resources.yml`
8. run focused workflow checks
9. announce cutover
10. deprecate old URLs on a planned timeline
