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
    
    volScalarField& he = thermo_.he();

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
            dimensionedScalar(dimless, 1)
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
    //  Find the fluxes in the equations
    //
    //================================================================

    //rho_pos = interpolate(rho, pos());
    rho_pos = interpolate(rho, pos(), rho.name());
    //rho_neg = interpolate(rho, neg());
    rho_neg = interpolate(rho, neg(), rho.name());
    
    const dimensionedScalar rhoSmall("rhoSmall", rho_.dimensions(), SMALL);
    
    //
    //  Calculate momentum fluxes
    ///
    const volVectorField rhoU(rho*U);
    rhoU_pos = interpolate(rhoU, pos(), U.name());
    rhoU_neg = interpolate(rhoU, neg(), U.name());

    //U_pos = surfaceVectorField::New("U_pos", rhoU_pos()/rho_pos());
    U_pos   = interpolate(U, pos(), U.name());
    //U_neg = surfaceVectorField::New("U_neg", rhoU_neg()/rho_neg());
    U_neg   = interpolate(U, neg(), U.name());
    
    //
    //  Calculate the energy fluxes
    //
    const volScalarField rhohe(rho*he);
    const volScalarField rhoHE(rhohe+rho*K);
    const volScalarField HE(he+K);    
    
    rhoHE_pos = interpolate(rhoHE, pos(), U.name());
    rhoHE_neg = interpolate(rhoHE, neg(), U.name());  
    
    const volScalarField& T = thermo.T();

    const volScalarField rPsi("rPsi", 1.0/thermo.psi());
    const surfaceScalarField rPsi_pos(interpolate(rPsi, pos(), T.name()));
    const surfaceScalarField rPsi_neg(interpolate(rPsi, neg(), T.name()));

    //p_pos = surfaceScalarField::New("p_pos", rho_pos()*rPsi_pos);
    p_pos   = interpolate(p, pos(), p.name());
    //p_neg = surfaceScalarField::New("p_neg", rho_neg()*rPsi_neg);
    p_neg   = interpolate(p, neg(), p.name());
    
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
    
    Yi_pos.setSize(Y_.size());
    Yi_neg.setSize(Y_.size());
    
    forAll(Y_, i)
    {    
        Yi_pos[i] = interpolate(Y_[i], pos(), Y_[i].name());
        Yi_neg[i] = interpolate(Y_[i], neg(), Y_[i].name());
    }        
    
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
    
    
    const volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));
    
    const surfaceScalarField cSf_pos
    (
        "cSf_pos",
        interpolate(c, pos(), T.name())*mesh.magSf()
    );
    const surfaceScalarField cSf_neg
    (
        "cSf_neg",
        interpolate(c, neg(), T.name())*mesh.magSf()
    );
    
    /*
    const dimensionedScalar c0("c0", U.dimensions(), 343.0); // or whatever scale is reasonable
    const surfaceScalarField cSf_pos = c0 * mesh.magSf();
    const surfaceScalarField cSf_neg = c0 * mesh.magSf();
    */

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
    
    //const surfaceScalarField ap( "ap", c0 * mesh.magSf() );
    //const surfaceScalarField am( "am", c0 * mesh.magSf() );
    
    
    a_pos = surfaceScalarField::New
    (
        "a_pos",
        fluxScheme == "Tadmor"
          ? surfaceScalarField::New("a_pos", mesh, 0.5)
          : ap/(ap - am)
    );
    a_neg = surfaceScalarField::New("a_neg", 1.0 - a_pos());
    
    //a_pos = surfaceScalarField::New("a_pos", surfaceScalarField::New("a_pos", mesh, 0.5));
    //a_neg = surfaceScalarField::New("a_neg", surfaceScalarField::New("a_neg", mesh, 0.5));

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

    phi_ = aphiv_pos()*rho_pos() + aphiv_neg()*rho_neg();
    
}


// ************************************************************************* //
