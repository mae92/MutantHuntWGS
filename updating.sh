# Uploading a container to docker hub
# 
# Log in on https://hub.docker.com/
# 
# Click on Create Repository
# 
# Choose a name (e.g. mutant_hunter) and a description for your repository and click Create.
# 
# Open Terminal

#build container
docker build -t mutant_hunt_wgs - < Dockerfile

#look up images
docker images

#tag image
docker tag 5e574e9ad429 mellison/mutant_hunt_wgs:version1.1

#push image
docker push mellison/mutant_hunt_wgs:version1.1

# Thats it!!

# Pull command
docker pull mellison/mutant_hunt_wgs:version1.1