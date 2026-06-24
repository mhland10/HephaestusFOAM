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

\*---------------------------------------------------------------------------*/

#include "fluidMulticomponentHephaestusThermo.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    //defineTypeNameAndDebug(fluidMulticomponentHephaestusThermo, 0);
    //defineRunTimeSelectionTable(fluidMulticomponentHephaestusThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

/*
Foam::autoPtr<Foam::fluidMulticomponentHephaestusThermo>
Foam::fluidMulticomponentHephaestusThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    // This assumes you will register derived models
    auto cstrIter =
        fvMeshConstructorTablePtr_->cfind(phaseName);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown fluidMulticomponentHephaestusThermo type: "
            << phaseName << nl
            << exit(FatalError);
    }

    return autoPtr<fluidMulticomponentHephaestusThermo>
    (
        cstrIter()(mesh, phaseName)
    );
}
    */


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidMulticomponentHephaestusThermo::~fluidMulticomponentHephaestusThermo()
{}


// ************************************************************************* //
