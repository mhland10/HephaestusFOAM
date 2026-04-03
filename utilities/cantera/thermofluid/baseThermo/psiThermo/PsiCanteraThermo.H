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

Class
    Foam::PsiCanteraThermo

Description
    Thermo implementation based on compressibility, but specifically using 
    Cantera, a free thermochemistry module

SourceFiles
    PsiCanteraThermo.C

\*---------------------------------------------------------------------------*/

#ifndef PsiCanteraThermo_H
#define PsiCanteraThermo_H

#include "BasicThermo.H"

// Cantera headers
#include "cantera/thermo/IdealGasPhase.h" // defines class IdealGasPhase
#include "cantera/base/Solution.h"
#include "cantera/transport.h" // transport properties
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class PsiCanteraThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class PsiCanteraThermo
:
    public Thermo
{
    // Private Member Functions

        //- Calculate the thermo variables
        virtual void calculate();
        
        bool useEnthalpy_;
        
protected:

    // Create pointer for Cantera solution
    Cantera::shared_ptr<Cantera::Solution> solution_;

public:

    // Constructors

        //- Construct from mesh and phase name
        PsiCanteraThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        PsiCanteraThermo(const PsiCanteraThermo<Thermo>&) = delete;


    //- Destructor
    virtual ~PsiCanteraThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PsiCanteraThermo<Thermo>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PsiCanteraThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
