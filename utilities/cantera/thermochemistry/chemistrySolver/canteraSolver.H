/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::canteraSolver

Description
    A dummy solver for chemistry via Cantera

SourceFiles
    canteraSolver.C

\*---------------------------------------------------------------------------*/

#ifndef canteraSolver_H
#define canteraSolver_H

#include "chemistrySolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class canteraSolver Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class canteraSolver
:
    public chemistrySolver<ChemistryModel>
{
    // Private Data

        dictionary coeffsDict_;

        // Solver data
        mutable scalarField cTp_;


public:

    //- Runtime type information
    TypeName("canteraSolver");


    // Constructors

        //- Construct from thermo
        canteraSolver(const fluidMulticomponentThermo& thermo);


    //- Destructor
    virtual ~canteraSolver();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalar& p,
            scalar& T,
            scalarField& c,
            const label li,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "canteraSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
