/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hephaestusFOAM.

    hephaestusFOAM is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or (at your 
    option) any later version.
    
    NOTE: hephaestusFOAM is limited by ITAR/EAR to US Persons.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

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
        void calculate();

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
