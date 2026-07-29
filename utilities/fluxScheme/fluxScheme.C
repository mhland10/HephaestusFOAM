/*---------------------------------------------------------------------------*\
         ________        |
        // \\   \\      |    HephaestusFOAM
       //   \\   \\     |      OpenFOAM addition for aerothermal chemistry
      //     \\   \\    |
      ||      \\   \\   |    From the University of Texas at San Antonio
      \\       \\___\\  |         the Laboratory for Turbulence, Sensing, 
       \\      ||___||  |              and Intelligent Systems
       /\\    //   //   |
      /  \\  //   //    |    Author:  Matthew Holland
     /   /\\//___//     |              matthew.holland@my.utsa.edu
    /   /  ^     ^      |
    \__/                |
                        |
-------------------------------------------------------------------------------
License
    This file is from HephaestusFOAM.

    HephaestusFOAM is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HephaestusFOAM is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HephaestusFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fluxScheme.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fluxScheme::fluxScheme
(
    const Time& runTime,
    const fvMesh& mesh,
    fluidMulticomponentThermo& thermo,
    volScalarField& p,
    volVectorField& U,
    volScalarField& rho,
    PtrList<volScalarField>& Y
)
:
    mesh(mesh),
    runTime_(runTime),
    thermo_(thermo),
    U_(U),
    rho_(rho),
    p_(p),
    T_(thermo_.T()),
    he_(thermo_.he()),
    Y_(Y),
    phi_
    (
        linearInterpolate(rho_*U_) & mesh.Sf()
    ),
    K_(0.5*magSqr(U_)),
    useEnthalpy_(false),
    fluxSchemeSelection_
    (
        mesh.schemes().dict().lookupOrDefault<word>("fluxScheme", "Kurganov")
    ),  
    physicalProperties_
    (
        IOobject
        (
            "physicalProperties",
            runTime_.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    phiRho_
    (
        IOobject
        (
            "phiRho",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimMass/dimTime,
            0.0
        )
    ),
    phiRhoUplusP_
    (
        IOobject
        (
            "phiRhoUplusP",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "zero",
            dimPressure*dimArea,
            vector::zero
        )
    ),
    phiRhoHE_
    (
        IOobject
        (
            "phiRhoHE",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimEnergy/dimTime,
            0.0
        )
    ),
    lambdaMax_
    (
        IOobject
        (
            "lambdaMax",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimVelocity,
            0.0
        )
    ),
    Uf_
    (
        IOobject
        (
            "lambdaMax",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "zero",
            dimVelocity,
            vector::zero
        )
    )
{
    Info<< "\n\n[fluxScheme] Constructor START\n" << endl;

    //
    //  Set up species allocation
    //
    phiRhoYi_.setSize(Y_.size());

    forAll(phiRhoYi_, i)
    {
        phiRhoYi_.set( i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phiRhoY" + Foam::name(i),
                    runTime_.name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar( "zero", dimMass/dimTime, 0.0 )
            )
        );
    }
    
    //
    //  Check the thermodynamics type
    //
    if (thermo_.he().name()=="h")
    {
        useEnthalpy_ = true;
    }
    else if (thermo_.he().name()=="e")
    {
        useEnthalpy_ = false;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown energy type: " << thermo.he().name()
            << exit(FatalError);
    }

    Info<< "[fluxScheme] Energy type: " << (useEnthalpy_ ? "Enthalpy" : "Internal Energy") << endl;

    //
    //  Initialize cantera
    //
    initializeCantera();

    Info<< "[fluxScheme] Cantera initialized..." << endl;

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fluxScheme::~fluxScheme()
{}

// * * * * * * * * * * * * *   Public Member Functions   * * * * * * * * * * //

void fluxScheme::corrector()
{
    clearTemporaryFields();
}

void fluxScheme::topoChange()
{
    pos.clear();
    neg.clear();
    clearTemporaryFields();
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void fluxScheme::clearTemporaryFields()
{
    rho_pos.clear();
    rho_neg.clear();

    rhoU_pos.clear();
    rhoU_neg.clear();

    U_pos.clear();
    U_neg.clear();

    p_pos.clear();
    p_neg.clear();
    
    T_pos.clear();
    T_neg.clear();

    he_pos.clear();
    he_neg.clear();

    a_pos.clear();
    a_neg.clear();

    c_pos.clear();
    c_neg.clear();

    aSf.clear();

    aphiv_pos.clear();
    aphiv_neg.clear();

    devTau.clear();
    
    rhoHE_pos.clear();
    rhoHE_neg.clear();
    
    rhoYi_pos.clear();
    rhoYi_neg.clear();

    phiUp.clear();
    phiHEp.clear();
    phiP.clear();
    phiYi.clear();

}

// * * * * * * * * * * * * * * Cantera Functions  * * * * * * * * * * * * * * //

using namespace Cantera;

void fluxScheme::initializeCantera()
{
    Foam::Info << "PsiCanteraThermo constructor from mesh & phaseName" << Foam::nl;
    
    // Access the dictionary via the mesh
    const dictionary& runDict = mesh.lookupObject<dictionary>("physicalProperties");
    const dictionary& canteraDict = runDict.subDict("cantera");
    const dictionary& setDict = runDict.subDict("thermoType");
    Foam::word energyType = setDict.lookupOrDefault<word>("energy", "sensibleEnthalpy");
    useEnthalpy_ = (energyType == "sensibleEnthalpy" || energyType == "enthalpy");
    Foam::Info << "Use Enthalpy? " << useEnthalpy_ << Foam::nl;
    
    // Read some properties
    word yamlFile = canteraDict.lookupOrDefault<word>("dataFile", "gri30.yaml");
    word phase = canteraDict.lookupOrDefault<word>("phaseName", "gri30");
    Foam::Info << "Cantera YAML file: " << yamlFile << ", phase: " << phase << Foam::nl;
    
    // Create the solution with default (ambient) state3
    canteraSolution_ = newSolution(yamlFile, phase);
}

// ************************************************************************* //

}