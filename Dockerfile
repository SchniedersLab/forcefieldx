FROM ubuntu:22.04

RUN apt-get update
RUN apt-get install -y wget unzip curl zip
RUN apt-get install -y libfreetype6 fontconfig fonts-dejavu
RUN apt-get install -y python3-pip

# Download GraalVM JDK 22 for x64 or Arm64
RUN set -eux; \
    ARCH="$(dpkg --print-architecture)"; \
    case "${ARCH}" in \
       aarch64|arm64) \
         JDK=jdk-22_linux-aarch64_bin.tar.gz; \
         ;; \
       amd64|i386:x86-64) \
         JDK=jdk-22_linux-x64_bin.tar.gz; \
         ;; \
       *) \
         echo "Unsupported arch: ${ARCH}"; \
         exit 1; \
         ;; \
    esac; \
    echo ${JDK}; \
    wget https://download.oracle.com/java/22/archive/${JDK}; \
    tar xzf ${JDK};

ENV JAVA_HOME=/jdk-22
ENV PATH="${JAVA_HOME}/bin:${PATH}"

# Add requirements.txt and install them using pip3
COPY requirements.txt .
SHELL ["/usr/bin/bash", "-c"]
RUN pip3 install --no-cache-dir -r requirements.txt

# Become root and set up the user environment
USER root
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

# Create the User home directory and copy in the FFX distro
RUN mkdir $HOME
COPY . $HOME

# Set FFX_HOME and the CLASSPATH
ENV FFX_HOME $HOME
ENV FFX_BIN $FFX_HOME/bin
ENV FFX_LIB $FFX_HOME/lib
ENV CLASSPATH $FFX_LIB/
ENV PATH="${FFX_BIN}:${PATH}"

# Build Force Field X
RUN set -eux; \
  cd /home/ffx; \
  ./mvnw; 

# Download Java Jupyter Kernel
RUN set -eux; \
  wget https://github.com/padreati/rapaio-jupyter-kernel/releases/download/2.1.0/rapaio-jupyter-kernel-2.1.0.jar; 

# Install the rapaio-jupyter-kernel as $NB_USER
RUN chown -R $NB_UID $HOME
USER $NB_USER
RUN set -eux; \
  java -jar rapaio-jupyter-kernel-2.1.0.jar -i -auto; \
  ls /home/ffx/.local/share/jupyter/kernels;

# Set up the FFX Kotlin library
# The allows the fillowing "magic" in Kotlin notebooks.
# %use ffx
RUN mkdir $HOME/.jupyter_kotlin
RUN mkdir $HOME/.jupyter_kotlin/libraries
RUN cp $HOME/ipynb-kotlin/ffx.json $HOME/.jupyter_kotlin/libraries/.
RUN cp -R $HOME/lib $HOME/.jupyter_kotlin/.

# Launch the notebook server
WORKDIR $HOME
CMD ["jupyter", "notebook", "--ip", "0.0.0.0"]

