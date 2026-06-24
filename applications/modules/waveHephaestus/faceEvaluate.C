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

\*---------------------------------------------------------------------------*/

#include "waveHephaestus.H"

// * * * * * * * Helpful Functions in an Anonymous Namespace * * * * * * * * //

#include <string>

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::waveFluid::faceEvaluate
(
    const word& property,
    const surfaceScalarField& pFaces,
    const surfaceScalarField& TFaces,
    const PtrList<surfaceScalarField>& YFaces,
    surfaceScalarField& field
)   const
{   
    // Pull the species names
    const wordList& speciesNames = thermo_.species();

    // Pull the cantera thermo
    auto canteraThermo = canteraSolution_->thermo();

    // The outer loop is to look through all the faces
    forAll(pFaces, faceI)
    {
        //
        //  Set up the p, T, and Y state
        //
        const scalar p = pFaces[faceI];
        const scalar T = TFaces[faceI];

        std::stringstream composition;
        forAll(speciesNames, specieI)
        {
            if (specieI > 0){ composition << ","; }
            composition << speciesNames[specieI] << ":" << 
                YFaces[specieI][faceI];
        }
        std::string compositionString = composition.str();
        
        // Set the thermodynamic state/mixture
        canteraThermo->setState_TPY( T, p, compositionString );

        if (property=="c")
        {
            field[faceI] = canteraThermo->cp_mass();
        }
        else if (property=="o")
        {
            field[faceI] = canteraThermo->cv_mass();
        }
        else if (property=="a")
        {
            const scalar gamma = canteraThermo->cp_mass() / canteraThermo->cv_mass();
            field[faceI] = sqrt(gamma * pFaces[faceI] / canteraThermo->density());
        }
        else if (property=="h")
        {
            const scalar hT = canteraThermo->enthalpy_mass();
            canteraThermo->setState_TPY(298.15, p, compositionString);
            field[faceI] = hT - canteraThermo->enthalpy_mass();
        }
        else if (property=="e")
        {
            const scalar eT = canteraThermo->intEnergy_mass();
            canteraThermo->setState_TPY(298.15, p, compositionString);
            field[faceI] = eT - canteraThermo->intEnergy_mass();
        }
        else if (property=="d")
        {
            field[faceI] = canteraThermo->density();
        }

    }

    /*
    // Correct dimensions
    if (property=="a")
    {
        field.dimensions().reset(dimVelocity);
    }*/

}