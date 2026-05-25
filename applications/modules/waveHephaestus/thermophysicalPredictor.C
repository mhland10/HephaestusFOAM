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
#include "fvmDdt.H"
#include "fvcDiv.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::solvers::waveFluid::pressureWork
(
    const tmp<volScalarField::Internal>& work
) const
{
    if (mesh.moving())
    {
    
        surfaceScalarField workPhi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());
        
        return work + fvc::div(workPhi, p/rho, "div(workPhi,(p|rho))")();
    }
    else
    {
        return move(work);
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::waveFluid::speciesPredictor()
{

    reaction->correct();
    
    //
    //  Species equation fluxes
    //
    List<tmp<surfaceScalarField>> phiYi;
    phiYi.setSize(Y_.size());
    forAll(Y_, i)
    {
        phiYi[i] =
        aphiv_pos()*rhoYi_pos[i]
      + aphiv_neg()*rhoYi_neg[i];
    }

    forAll(Y_, i)
    {
        // Define the species volume scalar field
        //volScalarField Yi = rhoYi_[i] / rho_;
        //volScalarField Yi = Y_[i];
        volScalarField& Yi = Y_[i];
        
        if (thermo_.solveSpecie(i))
        {
            //
            //  Species equation
            //
            fvScalarMatrix rhoYiEqn
            (
                //fvm::ddt(rhoYi_[i])
                fvm::ddt(rho, Yi)
              + fvc::div(phiYi[i])
             ==
                //fvModels().source(rhoYi_[i])
                fvModels().source(rho, Yi)
              + reaction->R(Yi)
            );
            
            if (!inviscid)
            {
                rhoYiEqn -= thermophysicalTransport->divj(Yi);
            }
            
            //
            //  Species equation solve
            //
            rhoYiEqn.relax();

            fvConstraints().constrain(rhoYiEqn);
            
            rhoYiEqn.solve("Yi");
            //rhoYiEqn.solve();
            
            //fvConstraints().constrain(rhoYi_[i]);   // apply BCs / limits FIRST
            fvConstraints().constrain(Yi); 
            
            //Y_[i] = rhoYi_[i] / rho_;                  // THEN recover primitive
            //Y_[i].correctBoundaryConditions();
            Yi.correctBoundaryConditions();
        }
        else
        {
            //Y_[i].correctBoundaryConditions();
            Yi.correctBoundaryConditions();
        }
    }

    thermo_.normaliseY();
}

void Foam::solvers::waveFluid::thermophysicalPredictor()
{
    speciesPredictor();

    volScalarField& he = thermo_.he();
    
    //volScalarField HE("HE", he + K);  
    
    //
    //  Energy equation fluxes
    //
    surfaceScalarField phiHEp
    (
        "phiHEp",
        aphiv_pos()*rhoHE_pos
      + aphiv_neg()*rhoHE_neg
    );
    
    if (he.name() == "e")
    {
    
      // Pressure jump term
      surfaceScalarField phiP_jump("phiP_jump", aSf()*(p_pos() - p_neg()));
      
      // Pressure upwind term
      surfaceScalarField phiP_upwind("phiP_upwind", aphiv_pos()*p_pos() + aphiv_neg()*p_neg());
      
      // Combine
      surfaceScalarField phiP = phiP_upwind + phiP_jump ;
      
      // Apply moving mesh correction
      if (mesh.moving())
      {
          phiP += mesh.phi()*(a_pos()*p_pos() + a_neg()*p_neg());
      }

      phiHEp += phiP;
        
    }   
    
    //
    //  Energy equation definition
    //
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvc::div(phiHEp)
      + fvc::ddt(rho, K)
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
    
    // Pull out kinetic energy
    //he = HE - K;

    fvConstraints().constrain(he);

    thermo_.correct();
}




// ************************************************************************* //
