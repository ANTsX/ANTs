# Contributing to ANTs

Contributions to ANTs are welcome. The developers are grateful to the many contributors to ANTs, both those who have [contributed to the ANTs source code](https://github.com/ANTsX/ANTs/graphs/contributors) and those who have contributed documentation, testing, feedback and ideas.

## Contributing code

### 1. Plan contribution with ANTs team

Before doing a lot of work, we recommend that you open an issue or start a discussion thread to discuss your ideas. If fixing a bug, please identify the bug in an issue first, following the bug report template.

### 2. Make changes and test locally

Once your code compiles without error and runs satisfactorily, you can open a pull request. You can mention issues by their number, eg `See issue #1234`. If you're confident that you're completing an open issue, you can say `Fixes #1234`, this will close the associated issue when the PR is merged.

### 3. Open a pull request

All contributions must be made via pull request on a fork of the ANTs repository. The developers will review and give feedback.

Write your commit messages using the standard prefixes for ITK commit messages:

* BUG: Fix for runtime crash or incorrect result
* COMP: Compiler error or warning fix
* DOC: Documentation change
* ENH: New functionality
* PERF: Performance improvement
* STYLE: No logic impact (indentation, comments)
* WIP: Work In Progress not ready for merge (used before PR is ready, or if a draft PR needs feedback or assistance)

A Github continuous integration test will be run on the PR, confirming that it compiles on Linux and doesn't break basic functionality. This is a minimal test, please test new functionality and if at all possible do so on open data, so that the developers can reproduce results.


### Licensing of contributions

All accepted contributions will be incorporated into the ANTs source code and distributed according to the terms of the ANTs license.
