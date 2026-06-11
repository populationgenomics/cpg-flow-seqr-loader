FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.138.cpg1-1

ENV PYTHONDONTWRITEBYTECODE=1
ENV VERSION=0.2.6

WORKDIR /cpg_seqr_loader

COPY src src/
COPY LICENSE pyproject.toml README.md ./

# pip install but don't retain the cache files
# pysam (used by the dummy-proband gVCF synthesis script) is also installed explicitly,
# mirroring cpg-flow-align-genotype - the wheel bundles htslib, so no system packages are needed
RUN pip install --no-cache-dir . && \
    pip install --no-cache-dir pysam
