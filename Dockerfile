# Build wheel separately
FROM python:slim AS src
RUN pip install build
RUN apt-get update && \
    apt-get install -y --no-install-recommends git
COPY . /src/nibabies
RUN python -m build /src/nibabies

# Python to support legacy MCRIBS
FROM python:3.6.15-slim as pyenv
RUN pip install --no-cache-dir numpy nibabel scipy pandas numexpr contextlib2 \
    && cp /usr/lib/x86_64-linux-gnu/libffi.so.7* /usr/local/lib

# Ubuntu 22.04 LTS
FROM ubuntu:jammy-20221130
ENV DEBIAN_FRONTEND="noninteractive" \
    LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8"

# Prepare environment
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    apt-utils \
                    autoconf \
                    build-essential \
                    bzip2 \
                    ca-certificates \
                    curl \
                    git \
                    gnupg \
                    libtool \
                    locales \
                    lsb-release \
                    netbase \
                    pkg-config \
                    unzip \
                    xvfb && \
    curl -sSL https://deb.nodesource.com/setup_14.x | bash - && \
    apt-get install -y --no-install-recommends \
                    nodejs && \
    locale-gen en_US.UTF-8 && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Installing ANTs 2.3.4 (NeuroDocker build)
ENV ANTSPATH="/usr/lib/ants" \
    PATH="/usr/lib/ants:$PATH"
WORKDIR $ANTSPATH
RUN curl -sSL "https://dl.dropbox.com/s/gwf51ykkk5bifyj/ants-Linux-centos6_x86_64-v2.3.4.tar.gz" \
    | tar -xzC $ANTSPATH --strip-components 1

RUN GNUPGHOME=/tmp gpg --keyserver hkps://keyserver.ubuntu.com --no-default-keyring --keyring /usr/share/keyrings/linuxuprising.gpg --recv 0xEA8CACC073C3DB2A \
    && echo "deb [signed-by=/usr/share/keyrings/linuxuprising.gpg] https://ppa.launchpadcontent.net/linuxuprising/libpng12/ubuntu jammy main" > /etc/apt/sources.list.d/linuxuprising.list

# AFNI 2023.04.04
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
ENV PATH="/opt/afni-latest:$PATH" \
    AFNI_IMSAVE_WARNINGS="NO" \
    AFNI_PLUGINPATH="/opt/afni-latest"

# Install AFNI latest (neurodocker build)
ENV AFNI_DIR="/opt/afni"
RUN echo "Downloading AFNI ..." \
    && mkdir -p ${AFNI_DIR} \
    && curl -fsSL --retry 5 https://afni.nimh.nih.gov/pub/dist/tgz/linux_openmp_64.tgz \
       | tar -xz -C ${AFNI_DIR} --strip-components 1 \
    # Keep only what we use
    && find ${AFNI_DIR} -type f -not \( \
        -name "3dTshift" -or \
        -name "3dUnifize" -or \
        -name "3dAutomask" -or \
        -name "3dvolreg" \) -delete


# FSL 6.0.5.1
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           bc \
           dc \
           file \
           libfontconfig1 \
           libfreetype6 \
           libgl1-mesa-dev \
           libgl1-mesa-dri \
           libglu1-mesa-dev \
           libgomp1 \
           libice6 \
           libxcursor1 \
           libxft2 \
           libxinerama1 \
           libxrandr2 \
           libxrender1 \
           libxt6 \
           sudo \
           wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && echo "Downloading FSL ..." \
    && mkdir -p /opt/fsl \
    && curl -fsSL --retry 5 https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-6.0.5.1-centos7_64.tar.gz \
    | tar -xz -C /opt/fsl --strip-components 1 \
    --exclude "fsl/config" \
    --exclude "fsl/data/atlases" \
    --exclude "fsl/data/first" \
    --exclude "fsl/data/mist" \
    --exclude "fsl/data/possum" \
    --exclude "fsl/data/standard/bianca" \
    --exclude "fsl/data/standard/tissuepriors" \
    --exclude "fsl/doc" \
    --exclude "fsl/etc/default_flobs.flobs" \
    --exclude "fsl/etc/fslconf" \
    --exclude "fsl/etc/js" \
    --exclude "fsl/etc/luts" \
    --exclude "fsl/etc/matlab" \
    --exclude "fsl/extras" \
    --exclude "fsl/include" \
    --exclude "fsl/python" \
    --exclude "fsl/refdoc" \
    --exclude "fsl/src" \
    --exclude "fsl/tcl" \
    --exclude "fsl/bin/FSLeyes" \
    && find /opt/fsl/bin -type f -not \( \
        -name "applywarp" -or \
        -name "bet" -or \
        -name "bet2" -or \
        -name "convert_xfm" -or \
        -name "fast" -or \
        -name "flirt" -or \
        -name "fsl_regfilt" -or \
        -name "fslhd" -or \
        -name "fslinfo" -or \
        -name "fslmaths" -or \
        -name "fslmerge" -or \
        -name "fslroi" -or \
        -name "fslsplit" -or \
        -name "fslstats" -or \
        -name "imtest" -or \
        -name "mcflirt" -or \
        -name "melodic" -or \
        -name "prelude" -or \
        -name "remove_ext" -or \
        -name "susan" -or \
        -name "topup" -or \
        -name "zeropad" \) -delete \
    && find /opt/fsl/data/standard -type f -not -name "MNI152_T1_2mm_brain.nii.gz" -delete
ENV FSLDIR="/opt/fsl" \
    PATH="/opt/fsl/bin:$PATH" \
    FSLOUTPUTTYPE="NIFTI_GZ" \
    FSLMULTIFILEQUIT="TRUE" \
    FSLLOCKDIR="" \
    FSLMACHINELIST="" \
    FSLREMOTECALL="" \
    FSLGECUDAQ="cuda.q" \
    LD_LIBRARY_PATH="/opt/fsl/lib:$LD_LIBRARY_PATH"

# Install FreeSurfer
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

# Workbench
WORKDIR /opt
RUN curl -sSLO https://www.humanconnectome.org/storage/app/media/workbench/workbench-linux64-v1.5.0.zip && \
    unzip workbench-linux64-v1.5.0.zip && \
    rm workbench-linux64-v1.5.0.zip && \
    rm -rf /opt/workbench/libs_linux64_software_opengl /opt/workbench/plugins_linux64 && \
    strip --remove-section=.note.ABI-tag /opt/workbench/libs_linux64/libQt5Core.so.5
    # ABI tags can interfere when running on Singularity/Apptainer
ENV PATH="/opt/workbench/bin_linux64:$PATH" \
    LD_LIBRARY_PATH="/opt/workbench/lib_linux64:$LD_LIBRARY_PATH"

# Installing SVGO and bids-validator
RUN npm install -g svgo@^2.3 bids-validator@1.9.9 \
  && rm -rf ~/.npm ~/.empty /root/.npm

# ICA AROMA
WORKDIR /opt/ICA-AROMA
RUN curl -sSL "https://github.com/oesteban/ICA-AROMA/archive/v0.4.5.tar.gz" \
  | tar -xzC /opt/ICA-AROMA --strip-components 1 && \
  chmod +x /opt/ICA-AROMA/ICA_AROMA.py
ENV PATH="/opt/ICA-AROMA:$PATH" \
    AROMA_VERSION="0.4.5"

# Create a shared $HOME directory
RUN useradd -m -s /bin/bash -G users nibabies && chmod -R 777 /home/nibabies
WORKDIR /home/nibabies
ENV HOME="/home/nibabies" \
    LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH}"

# py39_2209.01
COPY --from=nipreps/miniconda@sha256:8894ca17e3c8ba963812a6876093463eab6b88871bcfe23f71ebc84cf38451db /opt/conda /opt/conda

RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ${HOME}/.bashrc && \
    echo "conda activate base" >> ${HOME}/.bashrc

# Set CPATH for packages relying on compiled libs (e.g. indexed_gzip)
ENV PATH="/opt/conda/bin:$PATH" \
    CPATH="/opt/conda/include:$CPATH" \
    PYTHONNOUSERSITE=1 \
    MKL_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    IS_DOCKER_8395080871=1 \
    CONDA_PYTHON="/opt/conda/bin/python"

# Convert3d
RUN conda install -y -n base \
    -c anaconda \
    -c conda-forge \
    convert3d=1.3.0 \
    && sync \
    && conda clean -afy; sync \
    && rm -rf ~/.conda ~/.cache/pip/*; sync \
    && ldconfig

# MCRIBS
COPY --from=nipreps/mcribs@sha256:6c7a8dedd61d0ead8c7c4a57ab158928c1c1d787d87dae33ab7ee43226fb1e0f /opt/MCRIBS/ /opt/MCRIBS
RUN apt-get update && apt-get install -y --no-install-recommends \
        libboost-dev \
        libeigen3-dev \
        libflann-dev \
        libgl1-mesa-dev \
        libglu1-mesa-dev \
        libssl-dev \
        libxt-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
ENV PATH="/opt/MCRIBS/bin:/opt/MCRIBS/MIRTK/MIRTK-install/bin:/opt/MCRIBS/MIRTK/MIRTK-install/lib/tools:${PATH}" \
    LD_LIBRARY_PATH="/opt/MCRIBS/lib:/opt/MCRIBS/ITK/ITK-install/lib:/opt/MCRIBS/VTK/VTK-install/lib:/opt/MCRIBS/MIRTK/MIRTK-install/lib:${LD_LIBRARY_PATH}" \
    MCRIBS_HOME="/opt/MCRIBS" \
    PYTHONPATH="/opt/MCRIBS/lib/python:$PYTHONPATH"

# Precaching atlases
COPY scripts/fetch_templates.py fetch_templates.py
RUN ${CONDA_PYTHON} -m pip install --no-cache-dir --upgrade templateflow && \
    ${CONDA_PYTHON} fetch_templates.py && \
    rm fetch_templates.py && \
    find $HOME/.cache/templateflow -type d -exec chmod go=u {} + && \
    find $HOME/.cache/templateflow -type f -exec chmod go=u {} +

# Install pre-built wheel
COPY --from=src /src/nibabies/dist/*.whl .
RUN ${CONDA_PYTHON} -m pip install --no-cache-dir $( ls *.whl )[all]

# Final settings
RUN ldconfig
WORKDIR /tmp
ARG BUILD_DATE
ARG VCS_REF
ARG VERSION
LABEL org.label-schema.build-date=$BUILD_DATE \
      org.label-schema.name="nibabies" \
      org.label-schema.description="NiBabies - NeuroImaging tools for babies" \
      org.label-schema.url="https://github.com/nipreps/nibabies" \
      org.label-schema.vcs-ref=$VCS_REF \
      org.label-schema.vcs-url="https://github.com/nipreps/nibabies" \
      org.label-schema.version=$VERSION \
      org.label-schema.schema-version="1.0"

COPY --from=pyenv /usr/local/lib/ /usr/local/lib/
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
ENTRYPOINT ["/opt/conda/bin/nibabies"]
