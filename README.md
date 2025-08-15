# hephaestusFOAM

HephaestusFOAM is an OpenFOAM (v12) addition that allows users superior aerothermochemistry solving along with improved multiresolution capability.

## Getting started

- Install hephaestusFOAM via:

```
cd <install location>
git remote add origin http://prometheus.utsarr.net/mholland/hephaestusfoam.git
git branch -M main
git push -uf origin main>
```

## To Do List

- [X] Create Hephaestus-specific application to run the modules in
- [ ] Add higher order schemes
- [ ] Get psiThermo Cantera derivative working
    - [X] Get Cantera working in the OpenFOAM environment
    - [X] Create copy of psiThermo
    - [ ] Create Cantera pointer in 
- [ ] Port v2412's AMR to v12
- [ ] Add psiThermo changes to rhoThermo
- [ ] Add multiple criteria for AMR
- [ ] 
