## README

This directory contains Hail scripts used for loading and processing data in the CPG Seqr-Loader workflow.

This workflow has been migrated from [production-pipelines](https://github.com/populationgenomics/production-pipelines), and in that repository there
were two submodules ([gnomad_methods](https://github.com/broadinstitute/gnomad_methods), and [seqr-loader-pipelines](https://github.com/broadinstitute/seqr-loading-pipelines)), and the seqr-loader workflow was importing methods from these submodules. This worked fine, but was subject to two key notes:

1. the submodule pins, at a specific commit, have not been updated in years, and updating them could lead to breaking changes in the workflow. As such, these submodules are reflective of a fixed point in time, and not in sync with the latest version of the code.
2. reliance on submodules creates a more complex build and install process, as the entire content both modules is cloned and copied into Docker builds.

To enable complete independence from other codebases, this repository has copied the relevant methods from these submodules into this repository. Instead of importing large chunks of a separate codebase, the exact content of the methods used in the seqr-loader workflow has been copied into this repository, and is now maintained here.
