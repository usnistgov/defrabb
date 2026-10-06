<!-- markdownlint-disable-file MD041 -->
@AGENTS.md

## Claude Code notes (maintainer workstation)

- Use the `snakemake` conda env (`~/miniforge3/envs/snakemake`); it is not on
  `PATH` by default.
- Session and design notes live in the Obsidian vault at
  `~/Documents/ndo-notes/projects/defrabb/` (`dev-docs/` for roadmap and design
  notes). Write session notes there, not in the repo.
- GitLab API: use the numeric project ID (`glab api projects/6652/...`); the
  instance's `/gitlab/` URL prefix breaks `glab -R`.
