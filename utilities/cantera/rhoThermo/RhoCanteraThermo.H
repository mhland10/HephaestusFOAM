/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

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
