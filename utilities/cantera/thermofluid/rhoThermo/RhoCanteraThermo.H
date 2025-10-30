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

Class
    Foam::RhoCanteraThermo

Description
    Thermo implementation based on density

SourceFiles
    RhoCanteraThermo.C

\*---------------------------------------------------------------------------*/

#ifndef RhoCanteraThermo_H
#define RhoCanteraThermo_H

#include "BasicThermo.H"

// Cantera headers
#include "cantera/thermo/IdealGasPhase.h" // defines class IdealGasPhase
#include "cantera/base/Solution.h"
//#include "cantera/kinetics/GasKinetics.h"
//#include "cantera/transport.h" // transport properties
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RhoCanteraThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class RhoCanteraThermo
:
    public BaseThermo
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();
        
        bool useEnthalpy_;
        
protected:

    // Create pointer for Cantera solution
    Cantera::shared_ptr<Cantera::Solution> solution_;

public:

    // Constructors

        //- Construct from mesh and phase name
        RhoCanteraThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        RhoCanteraThermo(const RhoCanteraThermo<BaseThermo>&) = delete;


    //- Destructor
    virtual ~RhoCanteraThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const RhoCanteraThermo<BaseThermo>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "RhoCanteraThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
