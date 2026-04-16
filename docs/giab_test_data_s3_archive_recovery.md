# giab-test-data archive recovery

`giab-data/giab-test-data` currently contains objects in `INTELLIGENT_TIERING` with `ArchiveStatus=DEEP_ARCHIVE_ACCESS`. Public HTTPS reads against those objects fail with `403 InvalidObjectState`, which matches the pipeline download failure mode.

Verified examples:

- `giab-test-data/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz`
- `giab-test-data/GRCh38_chr21.fa.gz`
- `giab-test-data/CHM13v2.0/reference/chm13v2.0_maskedY_toupper.fa.gz`

## audit the referenced pipeline resources

```bash
python scripts/manage_giab_test_data_s3.py audit \
  --resources-file config/resources.yml \
  --check-http \
  --fail-if-archived
```

This reads the `giab-test-data` URLs from `config/resources.yml`, checks the public object metadata with unsigned `head-object`, and optionally probes the HTTPS URL with a 1-byte GET. Objects stuck in archive access will report `InvalidObjectState`.

To scan the whole prefix instead of only the resources referenced by the workflow:

```bash
python scripts/manage_giab_test_data_s3.py audit \
  --all-prefix-objects \
  --bucket giab-data \
  --prefix giab-test-data/
```

## request restore for archived objects

```bash
python scripts/manage_giab_test_data_s3.py restore \
  --resources-file config/resources.yml \
  --tier Standard
```

Notes:

- For `INTELLIGENT_TIERING` archive tiers, the restore request moves the object back into a directly readable tier.
- `Standard` restores are free for `INTELLIGENT_TIERING`, but deep archive restores can still take hours. This script is for recovery of the current bad state, not for low-latency access.
- The restore request requires AWS credentials with `s3:RestoreObject`.

## prevent future objects from becoming inaccessible

The bucket-side fix is to stop applying optional S3 Intelligent-Tiering archive tiers to `giab-test-data/`. Keeping the data in standard Intelligent-Tiering low-latency tiers (including Archive Instant Access) preserves public HTTPS access.

Inspect the current Intelligent-Tiering archive configuration:

```bash
python scripts/manage_giab_test_data_s3.py scan-config \
  --bucket giab-data \
  --prefix giab-test-data/
```

If the bucket has a prefix-scoped Intelligent-Tiering configuration that directly targets `giab-test-data/` (or a child prefix under it), the script can back it up and delete it:

```bash
python scripts/manage_giab_test_data_s3.py disable-config \
  --bucket giab-data \
  --prefix giab-test-data/ \
  --backup-dir results/s3_config_backups \
  --apply
```

The script intentionally refuses to auto-edit broader or tag-based configurations. Those need a bucket-admin review because S3 Intelligent-Tiering filters are positive matches only, so a bucket-wide rule cannot be safely auto-excluded for one prefix.

## current permission gap

The current `tibanna` IAM user can inspect public objects, but it cannot read or change the bucket lifecycle configuration:

```text
AccessDenied: s3:GetLifecycleConfiguration on arn:aws:s3:::giab-data
```

Bucket-admin credentials are therefore still required to apply the final policy/configuration fix.
