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

    surfaceVectorField U_posF   = interpolate(U, pos(), U.name());
    surfaceVectorField U_negF   = interpolate(U, neg(), U.name());
    
    //================================================================
    //
    //  Produce the initial flux values
    //
    //================================================================

    rho_pos = interpolate(rho, pos(), rho.name());
    rho_neg = interpolate(rho, neg(), rho.name());

    surfaceScalarField rho_posF = interpolate(rho, pos(), rho.name());
    surfaceScalarField rho_negF = interpolate(rho, neg(), rho.name());
    
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
        "he_pos", mesh, dimensionedScalar(thermo_.he().dimensions(), 0.0)
    );
    he_neg = surfaceScalarField::New
    (
        "he_neg", mesh, dimensionedScalar(thermo_.he().dimensions(), 0.0)
    );

    surfaceScalarField he_posF = interpolate(he, pos(), thermo_.he().name());
    surfaceScalarField he_negF = interpolate(he, neg(), thermo_.he().name());

    surfaceScalarField HE_posF = interpolate(HE, pos(), thermo_.he().name());
    surfaceScalarField HE_negF = interpolate(HE, neg(), thermo_.he().name());

    rhoHE_pos = interpolate(rhoHE, pos(), U.name());
    rhoHE_neg = interpolate(rhoHE, neg(), U.name());      

    surfaceScalarField rhoHE_posF = interpolate(rhoHE, pos(), U.name());
    surfaceScalarField rhoHE_negF = interpolate(rhoHE, neg(), U.name());   

    //
    //  Calculate the species fluxes
    //
    rhoYi_pos.setSize(Y_.size());
    rhoYi_neg.setSize(Y_.size());

    PtrList<volScalarField> rhoYi_(Y_.size());
    
    forAll(Y_, i)
    {
        rhoYi_.set( i, new volScalarField( "rhoYi_", rho_* Y_[i] ) );
    
        rhoYi_pos[i] = interpolate(rhoYi_[i], pos(), Y_[i].name());
        rhoYi_neg[i] = interpolate(rhoYi_[i], neg(), Y_[i].name());
    }

    //
    //  Calculate the flux values
    //
    const volScalarField rPsi("rPsi", 1.0/thermo.psi());
    const surfaceScalarField rPsi_pos(interpolate(rPsi, pos(), T.name()));
    const surfaceScalarField rPsi_neg(interpolate(rPsi, neg(), T.name()));
    surfaceScalarField phiv_pos("phiv_pos", 1.0* ( U_posF & mesh.Sf()) );
    surfaceScalarField phiv_neg("phiv_neg", 1.0* ( U_negF & mesh.Sf()) );

    //================================================================
    //
    //  Calculate the flux values
    //
    //================================================================

    //
    //  Correct the mass conservation values
    //
    faceEvaluate( "d", p_posF, T_posF, Yi_posF, rho_posF );
    faceEvaluate( "d", p_negF, T_negF, Yi_negF, rho_negF );

    //
    //  Correct the momentum conservation values
    //
    surfaceVectorField rhoU_posF(rho_posF*U_posF);
    surfaceVectorField rhoU_negF(rho_negF*U_negF);

    //
    //  Correct the energy conservation values
    //
    faceEvaluate( he.name(), p_posF, T_posF, Yi_posF, he_posF );
    faceEvaluate( he.name(), p_negF, T_negF, Yi_negF, he_negF );

    surfaceScalarField K_posF( "K_posF", 0.5 * magSqr(U_posF) );
    surfaceScalarField K_negF( "K_negF", 0.5 * magSqr(U_negF) );

    HE_posF = K_posF + he_posF;
    HE_negF = K_negF + he_negF;

    rhoHE_posF = HE_posF*rho_posF;
    rhoHE_negF = HE_negF*rho_negF;

    //
    //  Correct the species conersvation values
    //
    PtrList<surfaceScalarField> rhoYi_posF(Y_.size());
    PtrList<surfaceScalarField> rhoYi_negF(Y_.size());

    forAll(Y_, i)
    {
        rhoYi_posF.set( i, new surfaceScalarField( "rhoYi_posF",
                Yi_posF[i] * rho_posF
            )
        );
        rhoYi_negF.set( i, new surfaceScalarField( "rhoYi_negF",
                Yi_negF[i] * rho_negF
            )
        );
    }
    
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

    // Correct the sound speed
    surfaceScalarField c_posF("c_posF",interpolate(c, pos, T.name()));
    surfaceScalarField c_negF("c_negF",interpolate(c, neg, T.name()));
    faceEvaluate( "a", p_posF, T_posF, Yi_posF, c_posF );
    faceEvaluate( "a", p_negF, T_negF, Yi_negF, c_negF );
    cSf_pos = c_posF * mesh.magSf();
    cSf_neg = c_negF * mesh.magSf();

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

    aphiv_pos = surfaceScalarField::New("aphiv_pos", phiv_pos - aSf());
    aphiv_neg = surfaceScalarField::New("aphiv_neg", phiv_neg + aSf());

    if (fluxScheme == "HLLC"){

        /*
            Set all the necessary local scalar/vector fields
        */
        surfaceScalarField HE_pos = interpolate(HE, pos(), U.name());
        surfaceScalarField HE_neg = interpolate(HE, neg(), U.name());

        const surfaceScalarField rhoPos = interpolate(rho, pos(), rho.name());
        const surfaceScalarField rhoNeg = interpolate(rho, neg(), rho.name());

        const surfaceScalarField pPos = p_pos();
        const surfaceScalarField pNeg = p_neg();

        const surfaceVectorField rhoUPos = rhoU_pos();
        const surfaceVectorField rhoUNeg = rhoU_neg();

        const surfaceScalarField rhoHEPos = interpolate(rhoHE, pos(), U.name());
        const surfaceScalarField rhoHENeg = interpolate(rhoHE, neg(), U.name());

        surfaceVectorField nHat(mesh.Sf()/mesh.magSf());

        surfaceScalarField cn_pos(interpolate(c, pos, T.name()));
        surfaceScalarField cn_neg(interpolate(c, neg, T.name()));

        /*
            Process fields from the tmp fields
        */
        const surfaceVectorField Up = U_pos();
        const surfaceVectorField Un = U_neg();
        surfaceScalarField un_pos = (U_posF & nHat); 
        surfaceScalarField un_neg = (U_negF & nHat);  

        PtrList<surfaceScalarField> YiPos(Y_.size());
        PtrList<surfaceScalarField> YiNeg(Y_.size());

        forAll(Y_, i)
        {
            YiPos.set(i, new surfaceScalarField(Yi_pos[i]()));
            YiNeg.set(i, new surfaceScalarField(Yi_neg[i]()));
        }
        
        /*
            Set the wavespeed fields
        */
        surfaceScalarField SL_left = un_neg - c_negF;
        surfaceScalarField SL_right = un_pos - c_posF;
        surfaceScalarField SR_left = un_neg + c_negF;
        surfaceScalarField SR_right = un_pos + c_posF;

        surfaceScalarField SL("SL", min(SL_left, SL_right));
        surfaceScalarField SR("SR", min(SR_left, SR_right));

        surfaceScalarField SLmU_pos = SL - un_pos;
        surfaceScalarField SRmU_neg = SR - un_neg;
        surfaceScalarField termL = rho_posF * SLmU_pos;
        surfaceScalarField termR = rho_negF * SRmU_neg;

        surfaceScalarField denomStar = termL - termR;
        
        const dimensionedScalar denomStarSmall( "denomStarSmall", denomStar.dimensions(), SMALL );

        denomStar = ( 2.0*pos0(denomStar) - 1.0 ) * max( mag(denomStar), denomStarSmall );
        
        surfaceScalarField Sstar( "Sstar",
            ( p_negF - p_posF + rho_posF*un_pos*(SL - un_pos) - rho_negF*un_neg*(SR - un_neg) ) / denomStar
        );

        surfaceScalarField SLmSstar("SLmSstar", SL - Sstar);
        surfaceScalarField SRmSstar("SRmSstar", SR - Sstar);

        const dimensionedScalar starSmall("starSmall", dimVelocity, SMALL);

        SLmSstar = (2.0*pos0(SLmSstar)-1.0)*max(mag(SLmSstar), starSmall);
        SRmSstar = (2.0*pos0(SRmSstar)-1.0)*max(mag(SRmSstar), starSmall);
        
        surfaceScalarField rhoStarL("rhoStarL", max(rho_posF*(SL-un_pos)/SLmSstar, rhoSmall));
        surfaceScalarField rhoStarR("rhoStarR", max(rho_negF*(SR-un_neg)/SRmSstar, rhoSmall));
        
        surfaceVectorField rhoUStarL("rhoUStarL", rhoStarL*(U_posF+(Sstar-un_pos)*nHat));
        surfaceVectorField rhoUStarR("rhoUStarR", rhoStarR*(U_negF+(Sstar-un_neg)*nHat));
        
        surfaceScalarField rhoHEStarL("rhoHEStarL", rhoStarL*HE_posF);
        surfaceScalarField rhoHEStarR("rhoHEStarR", rhoStarR*HE_negF);
        if (he.name()=="e")
        {
            surfaceScalarField rhoSLmun_pos( "rhoSLmun_pos", rho_posF * ( SL - un_pos ) );
            surfaceScalarField rhoSRmun_neg( "rhoSLmun_neg", rho_negF * ( SR - un_neg ) );

            const dimensionedScalar starRhoSmall("starRhoSmall", dimDensity*dimVelocity, SMALL);

            rhoSLmun_pos = (2.0*pos0(rhoSLmun_pos)-1.0)*max(mag(rhoSLmun_pos), starRhoSmall);
            rhoSRmun_neg = (2.0*pos0(rhoSRmun_neg)-1.0)*max(mag(rhoSRmun_neg), starRhoSmall);

            rhoHEStarL += rhoStarL*(Sstar-un_pos)*(Sstar+p_posF/rhoSLmun_pos);
            rhoHEStarR += rhoStarR*(Sstar-un_neg)*(Sstar+p_negF/rhoSRmun_neg);
        }
        
        surfaceScalarField Fmass_L("Fmass_L", rho_posF*un_pos);
        surfaceScalarField Fmass_R("Fmass_R", rho_negF*un_neg);
        surfaceScalarField Fmass_Lstar("Fmass_Lstar", rho_posF*un_pos + SL*(rhoStarL - rho_posF));
        surfaceScalarField Fmass_Rstar("Fmass_Rstar", rho_negF*un_neg + SR*(rhoStarR - rho_negF));

        surfaceVectorField Fmom_L("Fmom_L", rhoU_posF*un_pos + p_posF*nHat);
        surfaceVectorField Fmom_R("Fmom_R", rhoU_negF*un_neg + p_posF*nHat);
        surfaceVectorField Fmom_Lstar("Fmom_Lstar", Fmom_L + SL*(rhoUStarL - rhoU_posF));
        surfaceVectorField Fmom_Rstar("Fmom_Rstar", Fmom_R + SR*(rhoUStarR - rhoU_negF));
        
        surfaceScalarField Fnrg_L("Fnrg_L", rhoHE_posF*un_pos);
        surfaceScalarField Fnrg_R("Fnrg_R", rhoHE_negF*un_neg);
        if (he.name()=="e")
        {
            Fnrg_L += p_posF*un_pos;
            Fnrg_R += p_negF*un_neg;
        }
        surfaceScalarField Fnrg_Lstar("Fnrg_Lstar", Fnrg_L + SL*(rhoHEStarL - rhoHE_posF));
        surfaceScalarField Fnrg_Rstar("Fnrg_Rstar", Fnrg_R + SR*(rhoHEStarR - rhoHE_negF));
        
        PtrList<surfaceScalarField> Yi_Face(Y_.size());
        forAll(Y_, i)
        {
            Yi_Face.set(i, new surfaceScalarField("YiFace_" + Foam::name(i),
                    pos0(SL)*Yi_posF[i]
                + (1.0 - pos0(SL))*pos0(Sstar)*Yi_posF[i]
                + (1.0 - pos0(SL))*(1.0 - pos0(Sstar))*pos0(SR)*Yi_negF[i]
                + (1.0 - pos0(SR))*Yi_negF[i]
                )
            );
        }
        
        surfaceScalarField F_massHLLC("F_massHLLC", Fmass_L);
        surfaceVectorField F_momHLLC( "F_momHLLC",  Fmom_L);
        surfaceScalarField F_nrgHLLC( "F_nrgHLLC",  Fnrg_L);
        PtrList<surfaceScalarField> F_specHLLC(Y_.size());

        F_massHLLC = pos0(SL)*Fmass_L
            + (1.0 - pos0(SL))*pos0(Sstar)*Fmass_Lstar
            + (1.0 - pos0(SL))*(1.0 - pos0(Sstar))*pos0(SR)*Fmass_Rstar
            + (1.0 - pos0(SR))*Fmass_R;

        F_momHLLC = pos0(SL)*Fmom_L
            + (1.0 - pos0(SL))*pos0(Sstar)*Fmom_Lstar
            + (1.0 - pos0(SL))*(1.0 - pos0(Sstar))*pos0(SR)*Fmom_Rstar
            + (1.0 - pos0(SR))*Fmom_R;
        
        F_nrgHLLC = pos0(SL)*Fnrg_L
            + (1.0 - pos0(SL))*pos0(Sstar)*Fnrg_Lstar
            + (1.0 - pos0(SL))*(1.0 - pos0(Sstar))*pos0(SR)*Fnrg_Rstar
            + (1.0 - pos0(SR))*Fnrg_R;
        
        surfaceScalarField massFlux = F_massHLLC;
        forAll(Y_, i)
        {
            F_specHLLC.set(i, new surfaceScalarField("F_Yi_"+Foam::name(i),
                    Yi_Face[i] * massFlux
                )
            );
        }
        
        //================================================================
        //
        //  Calculate the flux terms for equations
        //
        //================================================================

        //
        //  Calculate the mass fluxes
        //
        phi_ = F_massHLLC*mesh.magSf();

        
        //
        //  Calculate the momentum fluxes
        //
        phiUp = F_momHLLC*mesh.magSf();

        //
        //  Calculate the energy fluxes
        //
        phiHEp = F_nrgHLLC*mesh.magSf();

        //
        //  Calculate the species fluxes
        //
        phiYi.setSize(Y_.size());
        forAll(Y_, i)
        {
            phiYi[i] = F_specHLLC[i] * mesh.magSf();
        }

    }
    else{

        //================================================================
        //
        //  Calculate the flux terms for equations
        //
        //================================================================

        //
        //  Calculate the mass fluxes
        //
        phi_ = aphiv_pos()*rho_posF + aphiv_neg()*rho_negF;

        //
        //  Calculate the momentum fluxes
        //
        phiUp = aphiv_pos()*rhoU_posF + aphiv_neg()*rhoU_negF
                + (a_pos()*p_posF + a_neg()*p_negF)*mesh.Sf();

        //
        //  Calculate the energy fluxes
        //
        phiHEp = aphiv_pos()*rhoHE_posF + aphiv_neg()*rhoHE_negF;
        
        if (he.name() == "e")
        {
            // Pressure jump term
            tmp<surfaceScalarField> phiP_jump = aSf()*(p_posF - p_negF);
            
            // Pressure upwind term
            tmp<surfaceScalarField> phiP_upwind = aphiv_pos()*p_posF + aphiv_neg()*p_negF;
            
            // Combine
            phiP = phiP_upwind + phiP_jump ;
            
            // Apply moving mesh correction
            if (mesh.moving())
            {
                phiP.ref() += mesh.phi()*(a_pos()*p_posF + a_neg()*p_negF);
            }

            phiHEp.ref() += phiP;
            
        }

        if (he.name() == "h")
        {
        
            // Pressure jump term
            tmp<surfaceScalarField> phiP = aSf()*(p_posF - p_negF);

            // Apply moving mesh correction
            if (mesh.moving())
            {
                phiP.ref() += mesh.phi()*(a_pos()*p_posF + a_neg()*p_negF);
            }

            phiHEp.ref() += phiP;
            
        }
        
        //
        //  Calculate the species fluxes
        //
        phiYi.setSize(Y_.size());
        forAll(Y_, i)
        {
            phiYi[i] = aphiv_pos()*rhoYi_posF[i] + aphiv_neg()*rhoYi_negF[i];
        }

    }
    
}


// ************************************************************************* //
