# Install Stacks software (version 2.59)
The following script describes how I installed this software on our Klone node on Hyak

``` bash 

# Create a conda environment for stacks version 2.59
# environment location: /gscratch/merlab/software/miniconda3/envs/stacks_env
conda create -n stacks_env

# Activate the stacks environment
conda activate stacks_env
 

# Download stacks into the appropriate folder
cd /gscratch/merlab/software/miniconda3/envs/stacks_env
wget https://catchenlab.life.illinois.edu/stacks/source/stacks-2.59.tar.gz

# Install stacks
tar xfvz stacks-2.59.tar.gz
cd stacks-2.59
./configure --prefix=/gscratch/merlab/software/miniconda3/envs/stacks_env 
make
make install

# stacks is ready to use! Test it out by typing:
gstacks
populations

```
