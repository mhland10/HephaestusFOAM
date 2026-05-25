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

// Cantera headers
#include "cantera/thermo/IdealGasPhase.h" // defines class IdealGasPhase
#include "cantera/base/Solution.h"
#include "cantera/transport.h" // transport properties
#include <iostream>

// 
#include "PtrList.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * Helpful Functions in an Anonymous Namespace * * * * * * * * //

using namespace Foam;

std::string printCanteraMixture
(
    const PtrList<surfaceScalarField>& Y,
    const label facei,
    const wordList& speciesNames
)
{
    std::ostringstream oss;
    bool first = true;

    forAll(Y, i)
    {
        const scalar Yi = Y[i]()[facei]; 

        if (Yi <= 1e-14) continue;

        if (!first) oss << ",";
        oss << speciesNames[i] << ":" << Yi;
        first = false;
    }

    return oss.str();
}

using namespace Cantera;

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::waveFluid::facePropertyCalculate
(
    const word& propertyName,
    const surfaceScalarField& pf,
    const surfaceScalarField& Tf,
    const PtrList<surfaceScalarField>& Yif,
    surfaceScalarField& result
)   const
{
    //===========================================
    //
    //  Initialize the result
    //
    //===========================================
    
    // Print the composition
    //Info << "Cell " << celli << " composition: " << nl;
    const wordList& speciesNames = thermo_.species();
    auto canteraThermo = canteraSolution_->thermo();
    //Info << "    " << canteraComposition << nl;
    
    //============================================
    //
    //  Place data in the result field
    //
    //============================================
    
    forAll(pf, facei)
    {
        const scalar p = pf[facei];
        const scalar T = Tf[facei];

        //const scalar Yi = Yif[facei];
        
        // -----------------------------------------------------------
        //  Set up Cantera
        // -----------------------------------------------------------
        
        // Define the composition string for the current face
        std::string canteraComposition = printCanteraMixture(Yif, facei, speciesNames);
        
        // Set the Cantera state
        canteraThermo->setState_TPY(T, p, canteraComposition); 

        // --------------------------------------------------------
        // EOS evaluation (dispatch point)
        // --------------------------------------------------------

        scalar value = 0.0;

        if (propertyName == "rho")
        {
            value = canteraThermo->density();
        }
        else if (propertyName == "he")
        {
            if (useEnthalpy_)
            {
                value = canteraThermo->enthalpy_mass();
                canteraThermo->setState_TPY(298.15, p, canteraComposition);
                value -= canteraThermo->enthalpy_mass();
            }
            else
            {
                value = canteraThermo->intEnergy_mass();
                canteraThermo->setState_TPY(298.15, p, canteraComposition);
                value -= canteraThermo->intEnergy_mass();
            }
            // Restore original state
            canteraThermo->setState_TPY(T, p, canteraComposition);
        }
        else if (propertyName == "a")
        {
            value = sqrt((canteraThermo->cp_mass()/canteraThermo->cv_mass())*(p/canteraThermo->density()));
        }
        else if (propertyName == "Cp")
        {
            value = canteraThermo->cp_mass();
        }
        else if (propertyName == "Cv")
        {
            value = canteraThermo->cv_mass();
        }
        else if (propertyName == "psi")
        {
            value = canteraThermo->density()/p;
        }
        else if (propertyName == "rPsi")
        {
            value = p/canteraThermo->density();
        }
        else if (propertyName == "gamma")
        {
            value = (canteraThermo->cp_mass()/canteraThermo ->cv_mass());
        }
        else
        {
            FatalErrorInFunction
                << "Unknown thermo property: " << propertyName
                << exit(FatalError);
        }

        result[facei] = value;
    }

}