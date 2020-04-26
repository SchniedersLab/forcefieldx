#! /bin/bash

cd ~
wget https://download.java.net/openjdk/jdk11/ri/openjdk-11+28_linux-x64_bin.tar.gz
tar xvf openjdk-11+28_linux-x64_bin.tar.gz 
export JAVA_HOME="~/jdk-11"
export PATH="$JAVA_HOME/bin:$PATH"

wget https://github.com/SpencerPark/IJava/releases/download/v1.3.0/ijava-1.3.0.zip
unzip ijava-1.3.0.zip -d ijava-kernel   
cd ijava-kernel 
python3 install.py --sys-prefix 

ln -s ~/jdk-11/bin/java ~/bin/java
ln -s ~/jdk-11/bin/jshell ~/bin/jshell

