/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::multicomponentHephaestusThermo

Description
    Base-class for multi-component thermodynamic properties in HephaestusFoam

See also
    Foam::basicThermo

SourceFiles
    multicomponentHephaestusThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentHephaestusThermo_H
#define multicomponentHephaestusThermo_H

#include "basicThermo.H"
#include "MulticomponentHephaestusThermo.H"
#include "multicomponentThermo.H"
#include "speciesTable.H"
#include "GeometricFieldListSlicer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class multicomponentHephaestsThermo Declaration
\*---------------------------------------------------------------------------*/

class multicomponentHephaestusThermo
:
    virtual public multicomponentThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Destructor
    virtual ~multicomponentHephaestusThermo();

};


/*---------------------------------------------------------------------------*\
                Class multicomponentHephaestusThermo Declaration
\*---------------------------------------------------------------------------*/

class multicomponentHephaestusThermo::implementation
:
    virtual public multicomponentHephaestusThermo
{

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentHephaestusThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
