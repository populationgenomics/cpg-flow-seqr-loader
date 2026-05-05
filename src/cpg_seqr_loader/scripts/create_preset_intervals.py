#!/usr/bin/env python3

"""
Script to read in the gnomAD 4.1 data table, use it to obtain evenly distributed (variant-balanced) intervals

n.b. requires an additional Mito interval... (does it? we should produce with/without, we have a mito calling workflow)

requires some additional re-splitting of the resulting intervals e.g. to have clean breaks at contig ends

target interval counts - 50, 100, 500, 1000, 2000
"""

from cpg_utils import Path, to_path, hail_batch

...
