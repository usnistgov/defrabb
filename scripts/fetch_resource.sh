#!/usr/bin/env bash
#
# Fetch a DeFrABB resource that may be either a remote URL or a local
# filesystem path (issue #186). Centralizes the local-vs-remote logic for all
# download rules so config values can point at local files for offline / external
# use.
#
# Usage: fetch_resource.sh <source> <dest>
#
#   <source>  http(s)://, ftp(s)://, or s3:// URL, a file:// URI, or a local path
#   <dest>    output path. Compression handling is driven by the dest extension:
#               - dest ending in .gz/.tar.gz/.vcf.gz -> keep the bytes as-is
#               - any other dest                     -> decompress if the source
#                                                       is gzipped, else copy
#
set -euo pipefail

src="${1:?source required}"
dest="${2:?dest required}"

# file:// URIs are local
case "$src" in
file://*) src="${src#file://}" ;;
esac

is_remote() { [[ "$1" =~ ^(https?|ftps?|s3):// ]]; }

tmp="${dest}.fetch.tmp"
cleanup() { rm -f "$tmp"; }
trap cleanup EXIT

# Materialize the source bytes into a temp file (download or copy).
if is_remote "$src"; then
	# Retry remote downloads up to 3 times with exponential backoff
	max_retries=3
	retry_count=0
	retry_delay=5

	while [ $retry_count -lt $max_retries ]; do
		if curl -f -L --connect-timeout 120 --max-time 600 -o "$tmp" "$src"; then
			break
		fi

		retry_count=$((retry_count + 1))
		if [ $retry_count -lt $max_retries ]; then
			echo "fetch_resource: download failed (attempt $retry_count/$max_retries), retrying in ${retry_delay}s..." >&2
			sleep $retry_delay
			retry_delay=$((retry_delay * 2))
		else
			echo "fetch_resource: download failed after $max_retries attempts: $src" >&2
			exit 1
		fi
	done
else
	if [[ ! -f "$src" ]]; then
		echo "fetch_resource: local source not found: $src" >&2
		exit 1
	fi
	cp "$src" "$tmp"
fi

# Decide compression of the output from the destination extension.
case "$dest" in
*.gz)
	# Consumer expects gzipped bytes; pass through unchanged.
	mv "$tmp" "$dest"
	;;
*)
	# Consumer expects plain bytes; decompress only if the source is gzipped.
	if gzip -t "$tmp" 2>/dev/null; then
		gunzip -c "$tmp" >"$dest"
	else
		mv "$tmp" "$dest"
	fi
	;;
esac
