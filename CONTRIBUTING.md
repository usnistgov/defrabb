# Contributing to the GIAB Variant Calling Benchmark Set Generation and Evaluation Framework

This guide is a work in progress and meant to help development team formalize development process.

## Steps for creating new milestones/ framework versions

1. Create new milestone with a brief high-level description for new functionality, improvements, bug fixes
1. Assign relevant issues to the milestone
1. Address issues on separate branches and merge onto the dev branch
1. create merge request onto master branch when ready to release branch, using `milestone v#.###` as merge request title - review, address comments or create new issues when relevant
1. Update CHANGELOG, copy change log text into milestone and merge request description
1. Close merge request
1. Create tag for milestone version (tag name: v#.###) including CHANGELOG text in release notes
1. Clean-up old branches
