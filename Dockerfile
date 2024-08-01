FROM ubuntu:22.04

RUN apt-get update
RUN apt-get install -y wget unzip curl
RUN apt-get install -y libfreetype6 fontconfig fonts-dejavu
RUN apt-get install -y python3-pip

# Download GraalVM JDK 21 for x64 or Arm64
RUN set -eux; \
    ARCH="$(dpkg --print-architecture)"; \
    case "${ARCH}" in \
       aarch64|arm64) \
         JDK=openjdk-21.0.2_linux-aarch64_bin.tar.gz; \
         ;; \
       amd64|i386:x86-64) \
         JDK=openjdk-21.0.2_linux-x64_bin.tar.gz; \
         ;; \
       *) \
         echo "Unsupported arch: ${ARCH}"; \
         exit 1; \
         ;; \
    esac; \
    echo ${JDK}; \
    wget https://download.java.net/java/GA/jdk21.0.2/f2283984656d49d69e91c558476027ac/13/GPL/${JDK}; \
    tar xzf ${JDK};

ENV JAVA_HOME=/jdk-21.0.2
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

# Set IJava kernel environment variables.
ENV IJAVA_CLASSPATH $CLASSPATH
ENV IJAVA_STARTUP_SCRIPTS_PATH $FFX_BIN/startup.jshell

# Download Maven, unpack and build Force Field X
RUN set -eux; \
  cd /home/ffx; \
  echo $JAVA_HOME; \
  echo $PATH; \
  wget https://dlcdn.apache.org/maven/maven-3/3.9.8/binaries/apache-maven-3.9.8-bin.tar.gz; \
  tar xvzf apache-maven-3.9.8-bin.tar.gz; \
  apache-maven-3.9.8/bin/mvn; \
  rm apache-maven-3.9.8-bin.tar.gz; \
  rm -rf apache-maven-3.9.8;

# Download, Unpack and install the IJava kernel
# RUN set -eux; \
#  wget https://github.com/SpencerPark/IJava/releases/download/v1.3.0/ijava-1.3.0.zip; \
#  unzip ijava-1.3.0.zip -d ijava-kernel; \
#  cd ijava-kernel; \
#  CLASSPATH=$(echo /home/ffx/lib/* | tr ' ' ':'); \
#  python3 install.py --sys-prefix --startup-scripts-path $IJAVA_STARTUP_SCRIPTS_PATH --classpath $CLASSPATH;


# Download the kernel release
ENV PATH="${HOME}/.jbang/bin:${HOME}/.local/bin:${PATH}"
RUN set -eux; \
  curl -Ls https://sh.jbang.dev | bash -s - app setup; \
  jbang version; \
  jbang trust add https://github.com/jupyter-java/; \
  jbang install-kernel@jupyter-java rapaio --enable-preview --java=21 --jupyter-kernel-dir=/usr/local/share/jupyter/kernels --name="Java";

# Set up the FFX Kotlin library
# The allows the fillowing "magic" in Kotlin notebooks.
# %use ffx
RUN mkdir $HOME/.jupyter_kotlin
RUN mkdir $HOME/.jupyter_kotlin/libraries
RUN cp $HOME/ipynb-kotlin/ffx.json $HOME/.jupyter_kotlin/libraries/.
RUN cp -R $HOME/lib $HOME/.jupyter_kotlin/.

RUN chown -R $NB_UID $HOME
USER $NB_USER

# Launch the notebook server
WORKDIR $HOME
CMD ["jupyter", "notebook", "--ip", "0.0.0.0"]

