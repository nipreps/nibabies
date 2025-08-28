# NiBabies Docker Container Image distribution
#
# MIT License
#
# Copyright (c) 2023 The NiPreps Developers
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Ubuntu 22.04 LTS - Jammy
ARG BASE_IMAGE=ubuntu:jammy-20240405

# NiBabies wheel
FROM ghcr.io/astral-sh/uv:python3.12-alpine AS src
RUN apk add git
COPY . /src
RUN uvx --from build pyproject-build --installer uv -w /src

# Older Python to support legacy MCRIBS
FROM python:3.6.15-slim as pyenv
RUN pip install --no-cache-dir numpy nibabel scipy pandas numexpr contextlib2 \
    && cp /usr/lib/x86_64-linux-gnu/libffi.so.7* /usr/local/lib

# Intermediate step with utilities for downloading packages
FROM ${BASE_IMAGE} as downloader
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    binutils \
                    bzip2 \
                    ca-certificates \
                    curl \
                    unzip && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# AFNI
FROM downloader as afni
# The download link can point to newer releases
# As a safeguard, take advantage of Docker caching, and
# Bump the date to current to update AFNI
RUN echo "2025.06.12"
RUN mkdir -p /opt/afni-latest \
    && curl -fsSL --retry 5 https://afni.nimh.nih.gov/pub/dist/tgz/linux_openmp_64.tgz \
    | tar -xz -C /opt/afni-latest --strip-components 1 \
    --exclude "linux_openmp_64/*.gz" \
    --exclude "linux_openmp_64/funstuff" \
    --exclude "linux_openmp_64/shiny" \
    --exclude "linux_openmp_64/afnipy" \
    --exclude "linux_openmp_64/lib/RetroTS" \
    --exclude "linux_openmp_64/lib_RetroTS" \
    --exclude "linux_openmp_64/meica.libs" \
    # Keep only what we use
    && find /opt/afni-latest -type f -not \( \
        -name "3dTshift" -or \
        -name "3dUnifize" -or \
        -name "3dAutomask" -or \
        -name "3dvolreg" \) -delete

# Micromamba
FROM downloader as micromamba

# Install a C compiler to build extensions when needed.
# traits<6.4 wheels are not available for Python 3.11+, but build easily.
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /
# Bump the date to current to force update micromamba
RUN echo "2025.06.12" && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
ENV MAMBA_ROOT_PREFIX="/opt/conda"
COPY env.yml /tmp/env.yml
COPY requirements.txt /tmp/requirements.txt
WORKDIR /tmp
RUN micromamba create -y -f /tmp/env.yml && \
    micromamba clean -y -a

ENV PATH="/opt/conda/envs/nibabies/bin:$PATH" \
    UV_USE_IO_URING=0
RUN npm install -g svgo@^3.2.0 bids-validator@1.14.10 && \
    rm -r ~/.npm

# Main container
FROM ${BASE_IMAGE} as nibabies
ENV DEBIAN_FRONTEND="noninteractive" \
    LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8"

# Prepare environment
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    apt-utils \
                    bc \
                    build-essential \
                    ca-certificates \
                    curl \
                    git \
                    gnupg \
                    libtool \
                    locales \
                    lsb-release \
                    netbase \
                    unzip \
                    xvfb \
                    # MCRIBS-required
                    libboost-dev \
                    libeigen3-dev \
                    libflann-dev \
                    libgl1-mesa-dev \
                    libglu1-mesa-dev \
                    libssl-dev \
                    libxt-dev \
                    zlib1g-dev && \
    locale-gen en_US.UTF-8 && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Configure PPAs for libpng12 and libxp6
RUN GNUPGHOME=/tmp gpg --keyserver hkps://keyserver.ubuntu.com --no-default-keyring --keyring /usr/share/keyrings/linuxuprising.gpg --recv 0xEA8CACC073C3DB2A \
    && GNUPGHOME=/tmp gpg --keyserver hkps://keyserver.ubuntu.com --no-default-keyring --keyring /usr/share/keyrings/zeehio.gpg --recv 0xA1301338A3A48C4A \
    && echo "deb [signed-by=/usr/share/keyrings/linuxuprising.gpg] https://ppa.launchpadcontent.net/linuxuprising/libpng12/ubuntu jammy main" > /etc/apt/sources.list.d/linuxuprising.list \
    && echo "deb [signed-by=/usr/share/keyrings/zeehio.gpg] https://ppa.launchpadcontent.net/zeehio/libxp/ubuntu jammy main" > /etc/apt/sources.list.d/zeehio.list

# Dependencies for AFNI; requires a discontinued multiarch-support package from bionic (18.04)
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           ed \
           gsl-bin \
           libglib2.0-0 \
           libglu1-mesa-dev \
           libglw1-mesa \
           libgomp1 \
           libjpeg62 \
           libpng12-0 \
           libxm4 \
           libxp6 \
           netpbm \
           tcsh \
           xfonts-base \
           xvfb \
    && curl -sSL --retry 5 -o /tmp/multiarch.deb http://archive.ubuntu.com/ubuntu/pool/main/g/glibc/multiarch-support_2.27-3ubuntu1.5_amd64.deb \
    && dpkg -i /tmp/multiarch.deb \
    && rm /tmp/multiarch.deb \
    && apt-get install -f \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && gsl2_path="$(find / -name 'libgsl.so.19' || printf '')" \
    && if [ -n "$gsl2_path" ]; then \
         ln -sfv "$gsl2_path" "$(dirname $gsl2_path)/libgsl.so.0"; \
    fi \
    && ldconfig

COPY --from=afni /opt/afni-latest /opt/afni-latest

# AFNI config
ENV PATH="/opt/afni-latest:$PATH" \
    AFNI_IMSAVE_WARNINGS="NO" \
    AFNI_PLUGINPATH="/opt/afni-latest"

# Install FreeSurfer (with Infant Module)
COPY --from=nipreps/freesurfer@sha256:3b895fc732a7080374a15c4f976510f39c0c48dc76c030ab27316febd5e419ee /opt/freesurfer /opt/freesurfer
ENV FREESURFER_HOME="/opt/freesurfer"
ENV SUBJECTS_DIR="$FREESURFER_HOME/subjects" \
    FUNCTIONALS_DIR="$FREESURFER_HOME/sessions" \
    MNI_DIR="$FREESURFER_HOME/mni" \
    LOCAL_DIR="$FREESURFER_HOME/local" \
    MINC_BIN_DIR="$FREESURFER_HOME/mni/bin" \
    MINC_LIB_DIR="$FREESURFER_HOME/mni/lib" \
    MNI_DATAPATH="$FREESURFER_HOME/mni/data" \
    FSL_DIR=${FSLDIR} \
    FREESURFER="/opt/freesurfer"
ENV PERL5LIB="$MINC_LIB_DIR/perl5/5.8.5" \
    MNI_PERL5LIB="$MINC_LIB_DIR/perl5/5.8.5" \
    PATH="$FREESURFER_HOME/bin:$FREESURFER_HOME/tktools:$MINC_BIN_DIR:$PATH"

# MCRIBS (required legacy python)
COPY --from=nipreps/mcribs@sha256:363a5c3f25dd96fd1306329659b89b61718b50b8a6fc82a7e7732fc19af0cbc9 /opt/MCRIBS/ /opt/MCRIBS
COPY --from=pyenv /usr/local/lib/ /usr/local/lib/
ENV PATH="/opt/MCRIBS/bin:/opt/MCRIBS/MIRTK/MIRTK-install/bin:/opt/MCRIBS/MIRTK/MIRTK-install/lib/tools:${PATH}" \
    LD_LIBRARY_PATH="/opt/MCRIBS/lib:/opt/MCRIBS/ITK/ITK-install/lib:/opt/MCRIBS/VTK/VTK-install/lib:/opt/MCRIBS/MIRTK/MIRTK-install/lib:/usr/local/lib:${LD_LIBRARY_PATH}" \
    MCRIBS_HOME="/opt/MCRIBS" \
    PYTHONPATH="/opt/MCRIBS/lib/python:$PYTHONPATH"

# Create a shared $HOME directory
RUN useradd -m -s /bin/bash -G users nibabies && chmod -R 777 /home/nibabies
WORKDIR /home/nibabies
ENV HOME="/home/nibabies" \
    LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH}"

COPY --from=micromamba /bin/micromamba /bin/micromamba
COPY --from=micromamba /opt/conda/envs/nibabies /opt/conda/envs/nibabies
ENV MAMBA_ROOT_PREFIX="/opt/conda"
RUN micromamba shell init -s bash && \
    echo "micromamba activate nibabies" >> $HOME/.bashrc
ENV PATH="/opt/conda/envs/nibabies/bin:$PATH" \
    CPATH="/opt/conda/envs/nibabies/include:$CPATH" \
    LD_LIBRARY_PATH="/opt/conda/envs/nibabies/lib:$LD_LIBRARY_PATH" \
    CONDA_PYTHON="/opt/conda/envs/nibabies/bin/python"

# FSL environment
ENV LANG="C.UTF-8" \
    LC_ALL="C.UTF-8" \
    PYTHONNOUSERSITE=1 \
    FSLDIR="/opt/conda/envs/nibabies" \
    FSLOUTPUTTYPE="NIFTI_GZ" \
    FSLMULTIFILEQUIT="TRUE" \
    FSLLOCKDIR="" \
    FSLMACHINELIST="" \
    FSLREMOTECALL="" \
    FSLGECUDAQ="cuda.q"

# Unless otherwise specified each process should only use one thread - nipype
# will handle parallelization
ENV MKL_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    IS_DOCKER_8395080871=1

# Precaching atlases
COPY scripts/fetch_templates.py fetch_templates.py
RUN ${CONDA_PYTHON} -m pip install --no-cache-dir --upgrade templateflow && \
    ${CONDA_PYTHON} fetch_templates.py && \
    rm fetch_templates.py && \
    find $HOME/.cache/templateflow -type d -exec chmod go=u {} + && \
    find $HOME/.cache/templateflow -type f -exec chmod go=u {} +

# Install pre-built wheel
COPY --from=src /src/dist/*.whl .
RUN ${CONDA_PYTHON} -m pip install --no-cache-dir $( ls *.whl )[telemetry,test]

# Facilitate Apptainer use
RUN find $HOME -type d -exec chmod go=u {} + && \
    find $HOME -type f -exec chmod go=u {} + && \
    rm -rf $HOME/.npm $HOME/.conda $HOME/.empty && \
    ldconfig

WORKDIR /tmp
ARG BUILD_DATE
ARG VCS_REF
ARG VERSION
LABEL org.label-schema.build-date=$BUILD_DATE \
      org.label-schema.name="fMRIPrep Lifespan" \
      org.label-schema.description="fMRIPrep Lifespan - fMRI processing tool from birth and on" \
      org.label-schema.url="https://github.com/nipreps/nibabies" \
      org.label-schema.vcs-ref=$VCS_REF \
      org.label-schema.vcs-url="https://github.com/nipreps/nibabies" \
      org.label-schema.version=$VERSION \
      org.label-schema.schema-version="1.0"

ENTRYPOINT ["/opt/conda/envs/nibabies/bin/nibabies"]
