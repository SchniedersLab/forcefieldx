### Remove the Docker build cache
  docker builder prune 
  
### Build an image for the architecture of the current platform 
  docker build --no-cache -t forcefieldx -f Dockerfile .

### Build an image for Linux/amd64 and Linux/arm64 and push it to Docker Hub
   #### 1. Create a new multi-platform builder
   docker buildx create --name ffxbuilder --driver docker-container --use  
   #### 2. Boot up the new builder instance
   docker buildx inspect --bootstrap
   #### 3. Build the image and push it to Docker Hub
   docker buildx build --push --platform linux/amd64,linux/arm64 -t mjschnie/forcefieldx:26.0.4 -f Dockerfile .

### Test the Docker image in Bash
  docker run -it --rm forcefieldx bash

### Test the Docker image locally by running jupyter notebook
  docker run -it --rm -p 8888:8888 forcefieldx

### Test the Docker image locally by running jupyter lab
  docker run -it --rm -p 8888:8888 forcefieldx jupyter lab --ip=0.0.0.0
