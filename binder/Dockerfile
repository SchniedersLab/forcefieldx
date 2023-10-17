FROM ubuntu:22.04

RUN apt-get update
RUN apt-get install -y wget unzip
RUN apt-get install -y libfreetype6 fontconfig fonts-dejavu
RUN apt-get install -y python3-pip

# Download GraalVM JDK 21 for x64 or Arm64
RUN set -eux; \
    ARCH="$(dpkg --print-architecture)"; \
    case "${ARCH}" in \
       aarch64|arm64) \
         JDK=graalvm-jdk-21_linux-aarch64_bin.tar.gz; \
         ;; \
       amd64|i386:x86-64) \
         JDK=graalvm-jdk-21_linux-x64_bin.tar.gz; \
         ;; \
       *) \
         echo "Unsupported arch: ${ARCH}"; \
         exit 1; \
         ;; \
    esac; \
    echo $JDK; \
    wget https://download.oracle.com/graalvm/21/archive/${JDK}; \
    tar xzf ${JDK};

ENV JAVA_HOME=/graalvm-jdk-21+35.1
ENV PATH="${JAVA_HOME}/bin:${PATH}"

# add requirements.txt, written this way to gracefully ignore a missing file
COPY . .

SHELL ["/usr/bin/bash", "-c"]

RUN ([ -f binder/requirements.txt ] \
    && pip3 install --no-cache-dir -r binder/requirements.txt) \
        || pip3 install --no-cache-dir -r requirements.txt

USER root

# Version of FFX to Download
ENV FFX_VERSION 1.0.0-beta

# Download Force Field X
ENV FFX_TAR ffx-$FFX_VERSION-bin.tar.gz
RUN wget https://ffx.biochem.uiowa.edu/$FFX_TAR; \
    tar xzf $FFX_TAR; \
    rm $FFX_TAR; \
    mv ffx-$FFX_VERSION ffx;

# Set FFX_HOME and the CLASSPATH
ENV FFX_HOME /ffx
ENV FFX_BIN $FFX_HOME/bin
ENV FFX_LIB $FFX_HOME/lib
ENV CLASSPATH $FFX_LIB/
ENV PATH="${FFX_BIN}:${PATH}"

# Set IJava kernel environment variables.
ENV IJAVA_CLASSPATH $CLASSPATH
ENV IJAVA_STARTUP_SCRIPTS_PATH $FFX_BIN/startup.jshell

# Download, Unpack and install the IJava kernel
RUN set -eux; \
  wget https://github.com/SpencerPark/IJava/releases/download/v1.3.0/ijava-1.3.0.zip; \
  unzip ijava-1.3.0.zip -d ijava-kernel; \
  cd ijava-kernel; \
  CLASSPATH=$(echo /ffx/lib/* | tr ' ' ':'); \
  python3 install.py --sys-prefix --startup-scripts-path $IJAVA_STARTUP_SCRIPTS_PATH --classpath $CLASSPATH;

# Set up the user environment
ENV NB_USER ffx
ENV SHELL /usr/bin/bash
ENV NB_UID 1000
ENV HOME /home/$NB_USER

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid $NB_UID \
    --no-create-home \
    --shell $SHELL \
    $NB_USER

RUN mkdir $HOME
# COPY . $HOME

# Set up the Java and Kotlin Notebook Directories
RUN cp -R /ffx/ipynb-kotlin $HOME/.
RUN cp -R /ffx/ipynb-java $HOME/.
RUN cp -R /ffx/lib $HOME/.

RUN chown -R $NB_UID $HOME

USER $NB_USER

# Launch the notebook server
WORKDIR $HOME
CMD ["jupyter", "notebook", "--ip", "0.0.0.0"]
