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
    Foam::fluidMulticomponentHephaestusThermo

Description
    Base-class for multi-component fluid thermodynamic properties in the
    Hephaestus architecture.

See also
    Foam::basicThermo

SourceFiles
    fluidMulticomponentHephaestusThermo.C

\*---------------------------------------------------------------------------*/

#ifndef fluidMulticomponentHephaestusThermo_H
#define fluidMulticomponentHephaestusThermo_H

#include "fluidThermo.H"
#include "multicomponentThermo.H"
#include "fluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fluidMulticomponentHephaestusThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidMulticomponentHephaestusThermo
:
    virtual public fluidMulticomponentThermo
{
public:

    //- Runtime type information
    TypeName("fluidMulticomponentHephaestusThermo");

    //- Declare run-time constructor selection tables
    /*
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidMulticomponentHephaestusThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );*/

    // Selectors
        /*
        //- Standard selection based on fvMesh
        static autoPtr<fluidMulticomponentHephaestusThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );
        */


    //- Destructor
    virtual ~fluidMulticomponentHephaestusThermo();


    // Member Functions

        //- Function for face evaluation at the solver level
        virtual surfaceScalarField faceEvaluate
        (
            const word& property,
            const surfaceScalarField& pFace,
            const surfaceScalarField& TFace,
            const PtrList<surfaceScalarField>& YFace
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
