# GitLab Runner Setup for defrabb CI

Guide for installing and registering a GitLab Runner on a NIST workstation
to execute the defrabb CI pipeline (`.gitlab-ci.yml`). Use this when:

- The runner host has its OS reinstalled (existing config wiped)
- A new workstation is being added as a defrabb runner
- The current runner is broken and needs to be replaced

The defrabb CI uses a **shell executor** (not Docker) and a **pre-built conda
env** to avoid per-run install overhead and to work around the NIST firewall's
block on `anaconda.org`.

## Quick Reference

| Item | Value |
|---|---|
| GitLab instance | `https://gitlab.nist.gov` |
| Project ID | 6652 |
| Project path | `bbd-human-genomics/defrabb` |
| Runner type | project (locked to project 6652) |
| Executor | shell |
| Job user | `gitlab-runner` system user |
| Conda env name | `defrabb-ci` |
| Conda env prefix | `/home/gitlab-runner/miniforge3/envs/defrabb-ci` |
| Channel mirror | `https://prefix.dev/{conda-forge,bioconda}` |

## Prerequisites

- **Ubuntu** host (other Linux distros will need package-manager substitutions)
- **sudo / root** on the host
- A **GitLab personal access token** with `api` scope, exported as `$GITLAB_PAT`
  or saved via `glab auth login`. Used to create the runner record server-side.
- The `glab` CLI (`https://gitlab.com/gitlab-org/cli`) — optional but
  recommended for the API steps; plain `curl` works too.

## Step 1 — Install gitlab-runner

Add the official GitLab apt repo and install the package. The .deb post-install
creates the `gitlab-runner` system user and installs/enables a systemd unit.

```sh
curl -L "https://packages.gitlab.com/install/repositories/runner/gitlab-runner/script.deb.sh" | sudo bash
sudo apt install -y gitlab-runner
```

Verify:

```sh
gitlab-runner --version    # expect 18.x or newer
id gitlab-runner           # confirms the system user exists
systemctl is-enabled gitlab-runner
```

## Step 2 — Create the runner record on GitLab

GitLab 16+ uses authentication-token-based runner registration. Create the
runner record server-side first; the response contains a `glrt-` token that
the local runner uses to authenticate.

```sh
TOKEN_FILE=$(mktemp -t glrt.XXXXXX) && chmod 600 "$TOKEN_FILE"

glab api "user/runners" -X POST \
  -F runner_type=project_type \
  -F project_id=6652 \
  -F "description=shell executor on <hostname> (Ubuntu)" \
  -F run_untagged=true \
  -F locked=true \
  | python3 -c "
import json, sys, os
d = json.load(sys.stdin)
print(f\"runner id: {d['id']}\")
with open(os.environ['TOKEN_FILE'], 'w') as f:
    f.write(d['token'])
print(f\"token saved to: {os.environ['TOKEN_FILE']}\")"
```

The token is sensitive but only useful to register a runner against this
specific runner record — if it leaks, delete that record and create a new one.
We'll consume it in Step 5 and remove the tmpfile.

## Step 3 — Install miniforge for the `gitlab-runner` user

The runner user needs its own conda/mamba so it can build envs and so
snakemake's `--use-conda` can create per-rule envs at job time. Install into
the runner user's home so file ownership is correct.

```sh
sudo -u gitlab-runner bash -c '
  cd /tmp && \
  curl -L -o miniforge.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" && \
  bash miniforge.sh -b -p /home/gitlab-runner/miniforge3 && \
  rm miniforge.sh && \
  /home/gitlab-runner/miniforge3/bin/conda init bash && \
  /home/gitlab-runner/miniforge3/bin/mamba --version
'
```

## Step 4 — Configure `.condarc` for the runner user

NIST blocks `anaconda.org` and `repo.anaconda.com` at the firewall, so we route
conda-forge and bioconda through `prefix.dev` and use the institutional proxy.

```sh
sudo -u gitlab-runner tee /home/gitlab-runner/.condarc > /dev/null <<'EOF'
custom_multichannels:
  conda-forge:
    - https://prefix.dev/conda-forge
  bioconda:
    - https://prefix.dev/bioconda

channels:
  - conda-forge
  - bioconda

disallow:
  - anaconda

auto_activate: false

proxy_servers:
  http: http://kkjl338c0s.proxy.cloudflare-gateway.com
  https: https://kkjl338c0s.proxy.cloudflare-gateway.com
EOF
```

Confirm conda picks it up:

```sh
sudo -u gitlab-runner /home/gitlab-runner/miniforge3/bin/conda config --show-sources
```

You should see `/home/gitlab-runner/.condarc` listed with the
`custom_multichannels` block and the proxy.

## Step 5 — Pre-create the `defrabb-ci` conda env

Pin the same versions the CI was previously installing each run, plus `black`.
Pre-building this env once means CI activates it instead of reinstalling every
pipeline (saves ~5 min per run).

The `-H` flag plus explicit `cd ~` matter: without them, mamba's worker
processes inherit the parent shell's CWD (typically your home directory) and
crash with `Permission denied` because the `gitlab-runner` user can't traverse
it.

```sh
sudo -u gitlab-runner -H bash -c 'cd ~ && /home/gitlab-runner/miniforge3/bin/mamba create -n defrabb-ci -y \
    -c conda-forge -c bioconda --strict-channel-priority \
    python==3.12.3 snakemake==8.30.0 snakefmt==0.10.2 apptainer==1.3.6 pytest==8.3.3 black'
```

Verify:

```sh
sudo -u gitlab-runner -H bash -lc 'cd ~ && \
    source /home/gitlab-runner/miniforge3/etc/profile.d/conda.sh && \
    conda activate defrabb-ci && \
    which snakemake snakefmt pytest black apptainer && \
    snakemake --version && apptainer --version'
```

All five tool paths should resolve under
`/home/gitlab-runner/miniforge3/envs/defrabb-ci/bin/`.

## Step 6 — Register the runner

Use the token from Step 2. The new auth-token flow needs only `--url` and
`--token`; description/tags/locked were already set when we created the
runner record via API.

```sh
sudo gitlab-runner register --non-interactive \
    --url "https://gitlab.nist.gov" \
    --token "$(cat $TOKEN_FILE)" \
    --executor shell \
    --description "shell executor on $(hostname) (Ubuntu)"
```

Then restart the service so it picks up the new config and starts polling:

```sh
sudo systemctl restart gitlab-runner
sudo systemctl status gitlab-runner --no-pager | head
```

Clean up the token tmpfile — no longer needed:

```sh
shred -u "$TOKEN_FILE" && unset TOKEN_FILE
```

## Step 7 — Verify the runner is online

```sh
glab api "projects/6652/runners" | python3 -c "
import json, sys
for r in json.load(sys.stdin):
    print(f\"id={r['id']} status={r['status']} online={r['online']} desc={r['description']}\")"
```

The new runner should show `status=online online=True`.

Trigger a pipeline either by pushing a commit to master or via API:

```sh
glab api "projects/6652/pipeline?ref=master" -X POST | python3 -c "
import json, sys; p=json.load(sys.stdin); print(f\"pipeline #{p['id']} status={p['status']}\")"
```

Then watch the jobs:

```sh
glab api "projects/6652/pipelines/<pipeline_id>/jobs" | python3 -c "
import json, sys
for j in json.load(sys.stdin):
    print(f\"  {j['name']:12s} [{j['status']}] dur={j.get('duration','-')}\")"
```

`smkformat` and `pytest` should pass within ~30 seconds. `run-test` takes
20-30 minutes and may surface real env issues unrelated to the runner setup
(see `truvari-env-debugging-reference.md`).

## Step 8 — Clean up old runner records

If you replaced an existing runner, delete the stale record so it stops
appearing in the project's runner list.

```sh
# List all runners on the project
glab api "projects/6652/runners"

# Delete a specific stale runner by id
glab api "runners/<old_runner_id>" -X DELETE
```

## Troubleshooting

### `sudo: a terminal is required to read the password`
The Bash tool / non-interactive shells can't prompt for sudo. Run sudo
commands from your interactive shell (in Claude Code, prefix with `! `).

### `Permission denied: '/home/<your-user>'` during mamba operations
`sudo -u gitlab-runner` inherits CWD from your shell. The `gitlab-runner`
user can't traverse your home, so any tool that calls `chdir(getcwd())`
crashes. **Always include `cd ~` or `cd /tmp`** before invoking conda:

```sh
sudo -u gitlab-runner -H bash -c 'cd ~ && <command>'
```

This affects manual sanity checks only — actual CI jobs run from
`/home/gitlab-runner/builds/<token>/0/<group>/<project>/`, which the runner
user owns.

### `filesystem error: cannot set current path` from libmamba
Same root cause as above. Same fix.

### `pytest` fails with `analysis.qmd does not exist`
The `analysis*` glob in `.gitignore` will hide the source `analysis.qmd` if
the `!analysis.qmd` exception line is missing. Confirm `git ls-files | grep
analysis.qmd` returns the file.

### `run-test` fails on bwapy / `pkg_resources`
Pre-existing Truvari env bug, not a runner issue. See
`docs/truvari-env-debugging-reference.md` and memory `truvari_env_outstanding.md`.

### Runner shows `online=False` after registration
Check the systemd unit:

```sh
sudo systemctl status gitlab-runner --no-pager
sudo journalctl -u gitlab-runner -n 50 --no-pager
```

Most common causes: wrong `--url`, expired token, firewall blocking outbound
HTTPS to `gitlab.nist.gov`.

## Optional Polish

These don't block CI but improve the runner experience:

### Clear the systemd unit-changed warning

```sh
sudo systemctl daemon-reload
```

(The .deb post-install can leave this warning on first install.)

### Allow parallel jobs

The default `concurrent = 1` in `/etc/gitlab-runner/config.toml` serializes
all jobs. Bumping to 3 lets `smkformat` / `pytest` / `run-test` run in
parallel — useful because the long `run-test` would otherwise hold up the
quick checks:

```sh
sudo sed -i 's/^concurrent = 1$/concurrent = 3/' /etc/gitlab-runner/config.toml
sudo systemctl restart gitlab-runner
```

### Cache per-rule conda envs across CI runs

The `run-test` job spends most of its 30 min rebuilding per-rule conda envs
from scratch each pipeline. Pointing snakemake's `--conda-prefix` at a stable
location lets it reuse envs:

```yaml
# in .gitlab-ci.yml, run-test script
- snakemake -p --use-conda --use-apptainer --cores 1 --verbose
    --conda-prefix /home/gitlab-runner/.snakemake-cache/conda
    --config "_happy_threads=1"
```

Make sure the directory is writable by `gitlab-runner` and pre-create it
once: `sudo -u gitlab-runner mkdir -p /home/gitlab-runner/.snakemake-cache/conda`.

## Updating the `defrabb-ci` env

When the pinned versions in this guide need to change (or the upstream
`.gitlab-ci.yml` previously installed different versions), recreate the env:

```sh
sudo -u gitlab-runner -H bash -c 'cd ~ && \
    /home/gitlab-runner/miniforge3/bin/mamba env remove -n defrabb-ci -y && \
    /home/gitlab-runner/miniforge3/bin/mamba create -n defrabb-ci -y \
        -c conda-forge -c bioconda --strict-channel-priority \
        python==<X> snakemake==<Y> snakefmt==<Z> apptainer==<W> pytest==<V> black'
```

Bump the version pins in the command above to match.

## Architecture Decisions

**Why shell executor over Docker.** The original Docker-executor setup needed
`privileged=true` to run apptainer-in-Docker for the PAV rules. Shell executor
avoids that security wrinkle and runs natively on the host. Tradeoff: less
isolation between jobs and the host, but acceptable for a single-project
workstation runner.

**Why a dedicated `gitlab-runner` system user.** The .deb default. Isolates CI
job processes from any real user account on the host. Some teams run the
runner as their own user for convenience (reuses an existing conda) but that
gives every CI job access to that user's home, S3 credentials, etc.

**Why a pre-built `defrabb-ci` env over per-run install.** The original
`.gitlab-ci.yml` ran `mamba install ...` in `before_script` on every pipeline
(~5 min). With the env pre-built, `before_script` becomes `conda activate
defrabb-ci` (~1 second). The env spec rarely changes; updating it manually
when versions bump is cheap.

**Why prefix.dev as the channel mirror.** NIST blocks `anaconda.org` and
`repo.anaconda.com` at the institutional firewall. `prefix.dev` mirrors
conda-forge and bioconda from inside the perimeter. The user's own conda
config uses the same redirect; we mirror it for the `gitlab-runner` user.
