### Build the Docker image
  docker builder prune 
  
  docker build --no-cache -t forcefieldx -f Dockerfile .

### Test the Docker image in Bash
  docker run -it --rm forcefieldx bash

### Test the Docker image locally by running jupyter notebook
  docker run -it --rm -p 8888:8888 forcefieldx

### Test the Docker image locally by running jupyter lab
  docker run -it --rm -p 8888:8888 forcefieldx jupyter lab --ip=0.0.0.0
