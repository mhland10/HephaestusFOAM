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

void Foam::solvers::waveFluid::fluxPredictor()
{

    //================================================================
    //
    //  Initialize references
    //
    //================================================================
    
    // Thermo fields
    volScalarField& he = thermo_.he();
    const volScalarField& T = thermo.T();

    //================================================================
    //
    //  Find the Direction of the Faces
    //
    //================================================================
    if (!pos.valid())
    {
        pos = surfaceScalarField::New
        (
            "pos",
            mesh,
            dimensionedScalar(dimless, 1.0)
        );

        neg = surfaceScalarField::New
        (
            "neg",
            mesh,
            dimensionedScalar(dimless, -1.0)
        );
    }
    
    //================================================================
    //
    //  Find the primitive construction
    //
    //================================================================
    
    // Set up pressure
    p_pos   = interpolate(p, pos(), p.name());
    p_neg   = interpolate(p, neg(), p.name());

    surfaceScalarField p_posF(interpolate(p, pos(), p.name()));
    surfaceScalarField p_negF(interpolate(p, neg(), p.name()));

    // Set up temperature
    T_pos   = interpolate(T, pos(), T.name());
    T_neg   = interpolate(T, neg(), T.name());

    surfaceScalarField T_posF(interpolate(T, pos(), T.name()));
    surfaceScalarField T_negF(interpolate(T, neg(), T.name()));
    
    // Set up species
    Yi_pos.setSize(Y_.size());
    Yi_neg.setSize(Y_.size());
    forAll(Y_, i)
    {    
        Yi_pos[i] = interpolate(Y_[i], pos(), Y_[i].name());
        Yi_neg[i] = interpolate(Y_[i], neg(), Y_[i].name());
    }  

    PtrList<surfaceScalarField> Yi_posF(Y_.size());
    PtrList<surfaceScalarField> Yi_negF(Y_.size());

    forAll(Y_, i)
    {
        Yi_posF.set( i, new surfaceScalarField(interpolate(Y_[i], pos(), Y_[i].name()) ) );
        Yi_negF.set( i, new surfaceScalarField(interpolate(Y_[i], neg(), Y_[i].name()) ) );
    }

    // Set up velocity
    U_pos   = interpolate(U, pos(), U.name());
    U_neg   = interpolate(U, neg(), U.name());
    
    //================================================================
    //
    //  Produce the initial flux values
    //
    //================================================================

    rho_pos = interpolate(rho, pos(), rho.name());
    rho_neg = interpolate(rho, neg(), rho.name());
    
    const dimensionedScalar rhoSmall("rhoSmall", rho_.dimensions(), SMALL);
    
    //
    //  Calculate momentum fluxes
    //
    const volVectorField rhoU(rho*U);
    rhoU_pos = interpolate(rhoU, pos(), U.name());
    rhoU_neg = interpolate(rhoU, neg(), U.name());
    
    //
    //  Calculate the energy fluxes
    //
    K = 0.5*magSqr(U);
    const volScalarField rhohe(rho*he);
    const volScalarField rhoHE(rhohe+rho*K);
    const volScalarField HE(he+K);  

    he_pos = surfaceScalarField::New
    (
        "he_pos", mesh,
        dimensionedScalar(thermo_.he().dimensions(), 0.0)
    );
    he_neg = surfaceScalarField::New
    (
        "he_neg", mesh,
        dimensionedScalar(thermo_.he().dimensions(), 0.0)
    );
    
    rhoHE_pos = interpolate(rhoHE, pos(), U.name());
    rhoHE_neg = interpolate(rhoHE, neg(), U.name());      

    //
    //  Calculate the species fluxes
    //
    rhoYi_pos.setSize(Y_.size());
    rhoYi_neg.setSize(Y_.size());
    
    forAll(Y_, i)
    {
        rhoYi_[i] = rho_ * Y_[i];
    
        rhoYi_pos[i] = interpolate(rhoYi_[i], pos(), Y_[i].name());
        rhoYi_neg[i] = interpolate(rhoYi_[i], neg(), Y_[i].name());
    }

    //
    //  Calculate the flux values
    //
    const volScalarField rPsi("rPsi", 1.0/thermo.psi());
    const surfaceScalarField rPsi_pos(interpolate(rPsi, pos(), T.name()));
    const surfaceScalarField rPsi_neg(interpolate(rPsi, neg(), T.name()));

    //================================================================
    //
    //  Calculate the flux values
    //
    //================================================================

    surfaceScalarField phiv_pos("phiv_pos", 1.0* ( U_pos() & mesh.Sf()) );
    surfaceScalarField phiv_neg("phiv_neg", 1.0* ( U_neg() & mesh.Sf()) );
    
    //================================================================
    //
    //  Correct for the mesh motion
    //
    //================================================================

    // Make fluxes relative to mesh-motion
    if (mesh.moving())
    {
        phiv_pos -= mesh.phi();
        phiv_neg -= mesh.phi();
    }
    
    //================================================================
    //
    //  Calculate the wave corrections
    //
    //================================================================

    volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));
        
    surfaceScalarField cSf_pos("cSf_pos",interpolate(c, pos, T.name())*mesh.magSf());
    surfaceScalarField cSf_neg("cSf_neg",interpolate(c, neg, T.name())*mesh.magSf());

    const dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0);
    
    const surfaceScalarField ap
    (
        "ap",
        max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
    );
    const surfaceScalarField am
    (
        "am",
        min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
    );    

    surfaceScalarField denom
    (
        "denom",
        ap - am
    );

    const dimensionedScalar denomSmall
    (
        "denomSmall",
        denom.dimensions(),
        SMALL
    );

    // Preserve sign while enforcing minimum magnitude
    denom = sign(denom)*max(mag(denom), denomSmall);

    a_pos = surfaceScalarField::New
    (
        "a_pos",
        fluxScheme == "Tadmor"
        ? surfaceScalarField::New("a_pos", mesh, 0.5)
        : ap/denom
    );
    a_neg = surfaceScalarField::New("a_neg", 1.0 - a_pos());

    phiv_pos *= a_pos();
    phiv_neg *= a_neg();
    
    aSf = surfaceScalarField::New
    (
        "aSf",
        fluxScheme == "Tadmor"
          ? -0.5*max(mag(am), mag(ap))
          : am*a_pos()
    );
      
    //aSf = surfaceScalarField::New("aSf",0.0*max(mag(am), mag(ap)));

    aphiv_pos = surfaceScalarField::New("aphiv_pos", phiv_pos - aSf());
    aphiv_neg = surfaceScalarField::New("aphiv_neg", phiv_neg + aSf());

    //================================================================
    //
    //  Calculate the flux terms for equations
    //
    //================================================================

    //
    //  Calculate the mass fluxes
    //
    phi_ = aphiv_pos()*rho_pos + aphiv_neg()*rho_neg;

    //
    //  Calculate the momentum fluxes
    //
    phiUp = aphiv_pos()*rhoU_pos() + aphiv_neg()*rhoU_neg() 
            + (a_pos()*p_pos() + a_neg()*p_neg())*mesh.Sf();

    //
    //  Calculate the energy fluxes
    //
    phiHEp = aphiv_pos()*rhoHE_pos + aphiv_neg()*rhoHE_neg;
    
    if (he.name() == "e")
    {
    
      // Pressure jump term
      tmp<surfaceScalarField> phiP_jump = aSf()*(p_pos() - p_neg());
      
      // Pressure upwind term
      tmp<surfaceScalarField> phiP_upwind = aphiv_pos()*p_pos() + aphiv_neg()*p_neg();
      
      // Combine
      phiP = phiP_upwind + phiP_jump ;
      
      // Apply moving mesh correction
      if (mesh.moving())
      {
          phiP.ref() += mesh.phi()*(a_pos()*p_pos() + a_neg()*p_neg());
      }

      phiHEp.ref() += phiP;
        
    }

    
    if (he.name() == "h")
    {
    
      // Pressure jump term
      tmp<surfaceScalarField> phiP = aSf()*(p_pos() - p_neg());

      // Apply moving mesh correction
      if (mesh.moving())
      {
          phiP.ref() += mesh.phi()*(a_pos()*p_pos() + a_neg()*p_neg());
      }

      phiHEp.ref() += phiP;
        
    }
    

    //
    //  Calculate the species fluxes
    //
    phiYi.setSize(Y_.size());
    forAll(Y_, i)
    {
        phiYi[i] = aphiv_pos()*rhoYi_pos[i] + aphiv_neg()*rhoYi_neg[i];
    }
    
}


// ************************************************************************* //
