



# Knowledge-based potentials for new coarse-grained protein model

## Goals & Theory

The main goals of the project was to
- determine the initial parameters for the potentials
- examine some possible ways to reduce the number of these parameters

### The model 


### The potentials
- Ramachandran Map potential
- Side-Chain united atom interaction potentials

### Way to reduce # of parameters
- Reduced SVD
- GMM (Gaussian Mixture Model)

## File Structure

Directories:
- my
- **my_utilities** - files containing classes & functions used in scripts and notebooks in top level directory
- **setup**

Also
- `0_S_[name].py` - scripts
- `1_D_[name].py` - jupyter notebooks containing actual plots and analysis

## Setup

Python vesion: `python 3.7`

```pip install -r ./setup/requirements.txt```

## Credits

- Multiple files and statistics used in this project have been obtained using **Bioshell** package
- Excluded volume masks have been created with help of **Pyrosetta**


## Author
Piotr Śmieja, 2023