/*---------------------------------------------------------------------------*\
                        |
        //=\\===\\      |    HephaestusFOAM
       //   \\   \\     |      OpenFOAM addition for aerothermal chemistry
      //     \\   \\    |
      ||      \\   \\   |    From the University of Texas at San Antonio
      \\       \\===\\  |         the Laboratory for Turbulence, Sensing, 
       \\      ||   ||  |              and Intelligent Systems
       /\\    //   //   |
      /  \\  //   //    |    Author:  Matthew Holland
     /   /\\//   //     |              matthew.holland@my.utsa.edu
    /   /  ^-----^      |
    \__/                |
                        |
-------------------------------------------------------------------------------
License
    This file is from HephaestusFOAM.

    HephaestusFOAM is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HephaestusFOAM is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HephaestusFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef makeFluidMulticomponentCanteraThermo_H
#define makeFluidMulticomponentCanteraThermo_H

#include "makeExternalThermo.H"
#include "fluidMulticomponentThermo.H"
//#include "../lnInclude/makeExternalThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFluidMulticomponentCanteraThermo(BaseThermo, Mixture, ThermoPhysics)               \
                                                                                               \
    defineExternalThermo(BaseThermo, Mixture, ThermoPhysics);                                  \
                                                                                               \
    addExternalThermo(BaseThermo, BaseThermo, Mixture, ThermoPhysics)
    
//addExternalThermo(fluidMulticomponentThermo, BaseThermo, Mixture);  
/*
#define makeFluidMulticomponentCanteraThermos(                                 \
    BaseBaseThermo, BaseThermo, Mixture)    \
                                                                               \
    defineExternalThermo(BaseThermo, Mixture);                                 \
                                                                               \
    addExternalThermo(basicThermo, BaseThermo, Mixture);                       \
    addExternalThermo(fluidThermo, BaseThermo, Mixture);                       \
    addExternalThermo(BaseBaseThermo, BaseThermo, Mixture);                    \
    addExternalThermo(fluidMulticomponentThermo, BaseThermo, Mixture);         \
    addExternalThermo(BaseThermo, BaseThermo, Mixture)
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
