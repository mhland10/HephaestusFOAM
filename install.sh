#!/bin/bash
set -e

echo "Environment ready"

##=================================================================================================
##
##  This script is for a clean installation of HephaestusFOAM onto OpenFOAM v12. If HephaestusFOAM
##  needs to be partially installed, one must currently install those specific modules in the 
##  installation order that appears in this script.
##
##  Note that OpenFOAM v12 must be completely installed.
##
##  Questions should be directed to matthew.holland@utsa.edu
##
##=================================================================================================

## Append hephaestusFOAM RC to bashrc
# Temporarily disable exit-on-error
set +e
source ~/hephaestusfoam/hephaestusrc
rc=$?
set -e

if [ $rc -ne 0 ]; then
    echo "Warning: hephaestusrc returned non-zero ($rc)"
fi

echo "sourced the hephaestusrc"

chmod -R 777 *

echo "Reset permissions"

##=================================================================================================
##
##  Install Cantera
##  
##=================================================================================================

echo "**Installing Cantera**"
echo " "

## Run Cantera install script
cd "$HEPH_CANT"
./cantera_install

echo " "
echo "**Cantera Install Complete**"
echo " "

##=================================================================================================
##
##  Install Utilities
##
##=================================================================================================

echo "**Installing Utilities**"
echo " "

## Run utilities install script
cd "$HEPH_UTIL"
./utils_install

echo " "
echo "**Utilities Installed**"
echo " "

##=================================================================================================
##
##  Install Data
##
##=================================================================================================

echo "**Installing Data**"
echo " "

## Run utilities install script
cd "$HEPH_DATA"
./data_install

echo " "
echo "**Data Installed**"
echo " "

##=================================================================================================
##
##  Install Applications
##
##=================================================================================================

echo "**Installing Applications**"
echo " "

## Run utilities install script
cd "$HEPH_APP"
./app_install

echo " "
echo "**Applications Installed**"
echo " "
