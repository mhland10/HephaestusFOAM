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
    Foam::MulticomponentHephaestusThermo

Description
    Multi-component HephaestusFoam thermo implementation

SourceFiles
    MulticomponentHephaestusThermo.C

\*---------------------------------------------------------------------------*/

#ifndef MulticomponentHephaestusThermo_H
#define MulticomponentHephaestusThermo_H

#include "MulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class MulticomponentHephaestusThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class MulticomponentHephaestusThermo
:
    public MulticomponentThermo<BaseThermo>
{
public:

    //- Function for face evaluation at the solver level
    virtual surfaceScalarField faceEvaluate
    (
        const word& property,
        const surfaceScalarField& pFace,
        const surfaceScalarField& TFace,
        const PtrList<surfaceScalarField>& YFace
    ) const = 0;

    virtual ~MulticomponentHephaestusThermo() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "MulticomponentHephaestusThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
