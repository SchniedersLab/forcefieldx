### Build the Docker image
  docker builder prune 
  
  docker build --no-cache -t ffxdocker -f Dockerfile .

### Test the Docker image in Bash
  docker run -it --rm ffxdocker bash

### Test the Docker image locally by running jupyter notebook
  docker run -it --rm -p 8888:8888 ffxdocker

### Test the Docker image locally by running jupyter lab
  docker run -it --rm -p 8888:8888 ffxdocker jupyter lab --ip=0.0.0.0
