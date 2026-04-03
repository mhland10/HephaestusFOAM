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

#include "waveFluid.H"
#include "fvmDdt.H"
#include "fvcDiv.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

/*
  TODO: Implement the multicomponent thermophysicalPredictor() method and
        modify to include the flux splitting methods.
*/

void Foam::solvers::waveFluid::thermophysicalPredictor()
{
    speciesPredictor()

    volScalarField& he = thermo_.he();
    
    // Define the total energy flux term
    surfaceScalarField phiHE
    (
        "phiHE",
        aphiv_pos()*(rho_pos()*(HE_pos + 0.5*magSqr(U_pos())))
      + aphiv_neg()*(rho_neg()*(HE_neg + 0.5*magSqr(U_neg())))
    );
    
    // Pressure jump term
    surfaceScalarField phiP_jump("phiP_jump", aSf()*(p_pos() - p_neg()));
    
    // Pressure upwind term
    surfaceScalarField phiP_upwind("phiP_upwind", aphiv_pos()*p_pos() + aphiv_neg()*p_neg());
    
    // Combine
    surfaceScalarField phiP = phiP_jump + phiP_upwind;
    
    // Apply moving mesh correction
    if (mesh.moving())
    {
        phiP += mesh.phi()*(a_pos()*p_pos() + a_neg()*p_neg());
    }
    
    //
    //  Energy equation definition
    //
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvc::div(phiHE)
      + fvc::ddt(rho, K)
      + pressureWork(he.name() == "e" ? fvc::div(phiP) : -dpdt)
     ==
        fvModels().source(rho, he)
      + reaction->Qdot()
    );

    if (!inviscid)
    {
        const surfaceScalarField devTauDotU
        (
            "devTauDotU",
            devTau() & (a_pos()*U_pos() + a_neg()*U_neg())
        );

        EEqn += thermophysicalTransport->divq(he) + fvc::div(devTauDotU);
    }
    
    //
    //  Energy equation solve
    //
    EEqn.relax();

    fvConstraints().constrain(EEqn);

    EEqn.solve();

    fvConstraints().constrain(he);

    thermo_.correct();
}

void Foam::solvers::waveFluid::speciesPredictor()
{

    reaction->correct();

    forAll(Y, i)
    {
        // Define the species volume scalar field
        volScalarField& Yi = Y_[i];
        
        // Define the species concentration fluxes
        surfaceScalarField phiYi
        ( 
          "phiYi", 
          aphiv_pos()*(rho_pos()*Yi_pos[i])
        + aphiv_neg()*(rho_neg()*Yi_neg[i]) 
        );
        
        
        if (thermo_.solveSpecie(i))
        {
            //
            //  Species equation
            //
            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + fvc::div(phiYi)
             ==
                reaction->R(Yi)
              + fvModels().source(rho, Yi)
            );
            
            if (!inviscid)
            {
                YiEqn -= thermophysicalTransport->divj(Yi); 
            }
            
            //
            //  Species equation solve
            //
            YiEqn.relax();

            fvConstraints().constrain(YiEqn);

            YiEqn.solve("Yi");

            fvConstraints().constrain(Yi);
        }
        else
        {
            Yi.correctBoundaryConditions();
        }
    }

    thermo_.normaliseY();
}


// ************************************************************************* //
