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
Control
    Aspects of this software may be subject to CUI//EXPT restrictions. 
    Ensure you are using the public version, or confirm you are a U.S. 
    citizen or otherwise authorized under 22 CFR 120.62.

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

#include "psiMulticomponentCanteraThermo.H"

#include "coefficientMulticomponentMixture.H"
#include "coefficientWilkeMulticomponentMixture.H"
#include "singleComponentMixture.H"
#include "multicomponentMixture.H"
#include "homogeneousMixture.H"

#include "forGases.H"

#include "makeFluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{  

    forCoeffGases
    (
        makeFluidMulticomponentThermos,
        psiCanteraThermo,
        psiMulticomponentCanteraThermo,
        coefficientMulticomponentMixture
    );
    forCoeffGases
    (
        makeFluidMulticomponentThermos,
        psiCanteraThermo,
        psiMulticomponentCanteraThermo,
        coefficientWilkeMulticomponentMixture
    );
    forGases
    (
        makeFluidMulticomponentThermo,
        psiMulticomponentCanteraThermo,
        singleComponentMixture
    );
    
}

// ************************************************************************* //
