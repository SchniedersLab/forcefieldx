FROM openjdk:11.0.7

RUN apt-get update
RUN apt-get install -y python3-pip

# add requirements.txt, written this way to gracefully ignore a missing file
COPY . .
RUN ([ -f requirements.txt ] \
    && pip3 install --no-cache-dir -r requirements.txt) \
        || pip3 install --no-cache-dir jupyter jupyterlab

USER root

# Download the kernel release
RUN wget https://github.com/SpencerPark/IJava/releases/download/v1.3.0/ijava-1.3.0.zip

# Unpack and install the kernel
RUN unzip ijava-1.3.0.zip -d ijava-kernel \
  && cd ijava-kernel \
  && python3 install.py --sys-prefix

# Download Force Field X
RUN wget https://ffx.biochem.uiowa.edu/ffx-1.0.0-beta-bin.tar.gz

# Unpack Force Field X
RUN tar xzf ffx-1.0.0-beta-bin.tar.gz 

# Set FFX_HOME and the CLASSPATH
ENV FFX_HOME /ffx-1.0.0-beta
ENV CLASSPATH /ffx-1.0.0-beta/bin/ffx-all.jar

# Set up the user environment

ENV NB_USER runffx
ENV NB_UID 1000
ENV HOME /home/$NB_USER

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid $NB_UID \
    $NB_USER

COPY . $HOME
RUN chown -R $NB_UID $HOME

USER $NB_USER

# Launch the notebook server
WORKDIR $HOME
CMD ["jupyter", "notebook", "--ip", "0.0.0.0"]
