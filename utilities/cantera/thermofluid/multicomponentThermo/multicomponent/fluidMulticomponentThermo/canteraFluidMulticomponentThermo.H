/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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
    Foam::canteraFluidMulticomponentThermo

Description
    Base-class for multi-component fluid thermodynamic properties that 
    specifically operates for Cantera implementations.

See also
    Foam::basicThermo

SourceFiles
    canteraFluidMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef canteraFluidMulticomponentThermo_H
#define canteraFluidMulticomponentThermo_H

#include "fluidThermo.H"
#include "multicomponentThermo.H"
#include "CanteraFluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class fluidMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class canteraFluidMulticomponentThermo
:
    virtual public fluidThermo,
    virtual public multicomponentThermo
{
public:

    //- Runtime type information
    TypeName("CanteraFluidMulticomponentThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        canteraFluidMulticomponentThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<canteraFluidMulticomponentThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~canteraFluidMulticomponentThermo();


    // Member Functions

        // Specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mui
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Dynamic viscosity [kg/m/s]
            virtual tmp<volScalarField> mui
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
