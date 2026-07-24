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

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::waveFluid::prePredictor()
{
    //
    //  Troubleshooting Readout
    //
    /*
    Info<< "min(T)    = " << min(thermo.T()).value() << nl;
    Info<< "max(T)    = " << max(thermo.T()).value() << nl;

    Info<< "min(p)    = " << min(p).value() << nl;
    Info<< "max(p)    = " << max(p).value() << nl;

    Info<< "min(rho)  = " << min(rho).value() << nl;
    Info<< "max(rho)  = " << max(rho).value() << nl;

    forAll(Y, i)
    {
        Info<< Y[i].name()
            << " min=" << min(Y[i]).value()
            << " max=" << max(Y[i]).value()
            << nl;
    }
    */

    fluxScheme_.predictor();

    fluxPredictor();

    reaction->correct();

    correctDensity();

    if (!inviscid && pimple.predictTransport())
    {
        momentumTransport->predict();
        thermophysicalTransport->predict();
    }
}


// ************************************************************************* //
