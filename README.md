# hephaestusFOAM

HephaestusFOAM is an OpenFOAM (v12) addition that allows users superior aerothermochemistry solving along with improved multiresolution capability.

## Getting started

1. If not already done so, [install OpenFOAM v12](https://openfoam.org/download/12-ubuntu/)

2. Pull the git repository into the user's home directory. For fresh installations use the following. For stable versions, select the latest release branch. Generally it is recommended that <install location> is ~.

```
cd <install location>
git clone <repository-url>
git branch -M <branch name>

```

3. Install HephaestusFOAM from scratch via 

```
cd <install location>/hephaestusfoam
./clean_install

```

4. Use the available tutorials to understand how HephaestusFOAM is integrated into OpenFOAM.

## To Do List

- [X] Create Hephaestus-specific application to run the modules in
- [ ] Add higher order schemes
- [X] Get psiThermo Cantera derivative working
    - [X] Get Cantera working in the OpenFOAM environment
    - [X] Create copy of psiThermo
    - [X] Create Cantera pointer in 
- [ ] Port v2412's AMR to v12
- [ ] Add psiThermo changes to rhoThermo
- [ ] Add multiple criteria for AMR
- [ ] 
