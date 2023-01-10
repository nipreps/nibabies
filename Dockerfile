# Ubuntu 20.04 LTS
FROM python:slim AS src
RUN pip install build
RUN apt-get update && \
    apt-get install -y --no-install-recommends git
COPY . /src/nibabies
RUN python -m build /src/nibabies

FROM ubuntu:focal-20221130
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
                    graphviz \
                    libtool \
                    locales \
                    lsb-release \
                    pandoc \
                    pandoc-citeproc \
                    pkg-config \
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

# # AFNI latest (neurodocker build)
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           apt-utils \
           ed \
           gsl-bin \
           libglib2.0-0 \
           libglu1-mesa-dev \
           libglw1-mesa \
           libgomp1 \
           libjpeg62 \
           libxm4 \
           netpbm \
           tcsh \
           xfonts-base \
           xvfb \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && curl -sSL --retry 5 -o /tmp/multiarch.deb http://archive.ubuntu.com/ubuntu/pool/main/g/glibc/multiarch-support_2.27-3ubuntu1.5_amd64.deb \
    && dpkg -i /tmp/multiarch.deb \
    && rm /tmp/multiarch.deb \
    && curl -sSL --retry 5 -o /tmp/libxp6.deb http://mirrors.kernel.org/debian/pool/main/libx/libxp/libxp6_1.0.2-2_amd64.deb \
    && dpkg -i /tmp/libxp6.deb \
    && rm /tmp/libxp6.deb \
    && curl -sSL --retry 5 -o /tmp/libpng.deb http://snapshot.debian.org/archive/debian-security/20160113T213056Z/pool/updates/main/libp/libpng/libpng12-0_1.2.49-1%2Bdeb7u2_amd64.deb \
    && dpkg -i /tmp/libpng.deb \
    && rm /tmp/libpng.deb \
    && apt-get install -f \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && gsl2_path="$(find / -name 'libgsl.so.19' || printf '')" \
    && if [ -n "$gsl2_path" ]; then \
         ln -sfv "$gsl2_path" "$(dirname $gsl2_path)/libgsl.so.0"; \
    fi \
    && ldconfig \
    && echo "Downloading AFNI ..." \
    && mkdir -p /opt/afni-latest \
    && curl -fsSL --retry 5 https://afni.nimh.nih.gov/pub/dist/tgz/linux_openmp_64.tgz \
    | tar -xz -C /opt/afni-latest --strip-components 1 \
    --exclude "linux_openmp_64/*.gz" \
    --exclude "linux_openmp_64/funstuff" \
    --exclude "linux_openmp_64/shiny" \
    --exclude "linux_openmp_64/afnipy" \
    --exclude "linux_openmp_64/lib/RetroTS" \
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

# Convert3D (neurodocker build)
RUN echo "Downloading Convert3D ..." \
    && mkdir -p /opt/convert3d-1.0.0 \
    && curl -fsSL --retry 5 https://sourceforge.net/projects/c3d/files/c3d/1.0.0/c3d-1.0.0-Linux-x86_64.tar.gz/download \
    | tar -xz -C /opt/convert3d-1.0.0 --strip-components 1 \
    --exclude "c3d-1.0.0-Linux-x86_64/lib" \
    --exclude "c3d-1.0.0-Linux-x86_64/share" \
    --exclude "c3d-1.0.0-Linux-x86_64/bin/c3d_gui"
ENV C3DPATH="/opt/convert3d-1.0.0" \
    PATH="/opt/convert3d-1.0.0/bin:$PATH"

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
    FSL_DIR=${FSLDIR}
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
RUN useradd -m -s /bin/bash -G users nibabies
WORKDIR /home/nibabies
ENV HOME="/home/nibabies"

# py39_2209.01
COPY --from=nipreps/miniconda@sha256:8894ca17e3c8ba963812a6876093463eab6b88871bcfe23f71ebc84cf38451db /opt/conda /opt/conda

RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Set CPATH for packages relying on compiled libs (e.g. indexed_gzip)
ENV PATH="/opt/conda/bin:$PATH" \
    CPATH="/opt/conda/include:$CPATH" \
    PYTHONNOUSERSITE=1 \
    MKL_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    IS_DOCKER_8395080871=1 \
    CONDA_PYTHON="/opt/conda/bin/python"

# Precaching atlases
COPY scripts/fetch_templates.py fetch_templates.py
RUN ${CONDA_PYTHON} -m pip install --no-cache-dir templateflow && \
    ${CONDA_PYTHON} fetch_templates.py && \
    rm fetch_templates.py && \
    find $HOME/.cache/templateflow -type d -exec chmod go=u {} + && \
    find $HOME/.cache/templateflow -type f -exec chmod go=u {} +

COPY --from=src /src/nibabies/dist/*.whl .
RUN /opt/conda/bin/python -m pip install --no-cache-dir $( ls *.whl )[all]

# ABI tags can interfere when running on Singularity/Apptainer
RUN strip --remove-section=.note.ABI-tag /usr/lib/x86_64-linux-gnu/libQt5Core.so.5

# Final settings
RUN ldconfig
WORKDIR /tmp
ARG BUILD_DATE
ARG VCS_REF
ARG VERSION
LABEL org.label-schema.build-date=$BUILD_DATE \
      org.label-schema.name="nibabies" \
      org.label-schema.description="nibabies - NeuroImaging tools for babies" \
      org.label-schema.url="https://github.com/nipreps/nibabies" \
      org.label-schema.vcs-ref=$VCS_REF \
      org.label-schema.vcs-url="https://github.com/nipreps/nibabies" \
      org.label-schema.version=$VERSION \
      org.label-schema.schema-version="1.0"

ENTRYPOINT ["/opt/conda/bin/nibabies"]
