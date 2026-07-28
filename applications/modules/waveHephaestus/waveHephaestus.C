/*---------------------------------------------------------------------------*\
                        |
        //=\\===\\      |    HephaestusFOAM
       //   \\   \\     |      OpenFOAM addition for aerothermal chemistry
      //     \\   \\    |
      ||      \\   \\   |    From the University of Texas at San Antonio
      \\       \\===\\  |         the Laboratory for Turbulence, Sensing, 
       \\      ||   ||  |              and Intelligent Systems
       /\\    //   //   |
      /  \\  //   //    |    Author:  Matthew Holland
     /   /\\//   //     |              matthew.holland@my.utsa.edu
    /   /  ^-----^      |
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

#include "waveHephaestus.H"

#include "fvMeshStitcher.H"
#include "localEulerDdtScheme.H"
#include "hydrostaticInitialisation.H"
#include "fvcMeshPhi.H"
#include "fvcVolumeIntegrate.H"
#include "fvcReconstruct.H"
#include "fvcSnGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(waveFluid, 0);
    addToRunTimeSelectionTable(solver, waveFluid, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::waveFluid::correctCoNum(const surfaceScalarField& amaxSf)
{
    const scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

    CoNum_ =
        0.5*gMax(sumAmaxSf/mesh.V().primitiveField())*runTime.deltaTValue();

    const scalar meanCoNum =
        0.5
       *(gSum(sumAmaxSf)/gSum(mesh.V().primitiveField()))
       *runTime.deltaTValue();

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::solvers::waveFluid::clearTemporaryFields()
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::waveFluid::waveFluid(fvMesh& mesh)
:
    fluidSolver(mesh),
    
    //
    //  Thermodynamics
    //
    thermoPtr_(fluidMulticomponentThermo::New(mesh)),

    thermo_(*thermoPtr_),
    
    //
    //  CFD Data
    //
    U_
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    rho_
    (
        IOobject
        (
            "rho",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo_.renameRho()
    ),
    
    p_(thermo_.p()),
    
    Y_(thermo_.Y()),    
    
    //
    //  Derived Data
    //   
    phi_
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho_*U_) & mesh.Sf()
    ),
     
    dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.name(),
            mesh
        ),
        mesh,
        dimensionedScalar(p_.dimensions()/dimTime, 0)
    ),

    K("K", 0.5*magSqr(U_)),
    
    //
    //  Formulation data
    //
    fluxSchemeName_
    (
        mesh.schemes().dict().lookupOrDefault<word>("fluxScheme", "Kurganov")
    ),

    fluxScheme_
    (
        runTime,
        mesh,
        thermo_,
        p_,
        U_,
        rho_,
        Y_
    ),
    
    //
    //  Thermophysical data structure
    //    
    physicalProperties_
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    
    inviscid(physicalProperties_.lookupOrDefault<bool>("inviscid", false)),
    
    momentumTransport
    (
      compressible::momentumTransportModel::New
        (
            rho_,
            U_,
            phi_,
            thermo_
        )
    ),
    
    thermophysicalTransport
    (
      fluidMulticomponentThermophysicalTransportModel::New
        (
            momentumTransport(),
            thermo_
        )
    ),
    
    reaction(combustionModel::New(thermo_, momentumTransport())),
    
    //
    //  Stored field data
    //
    thermo(thermo_),
    p(p_),
    rho(rho_),
    U(U_),
    phi(phi_),
    Y(Y_)
{

    Info<< "\n\n[waveFluid] Constructor START\n" << endl;
    
    //
    //  Check the thermodynamics type
    //
    thermo.validate(type(), "h", "e");
    
    Info<< "[waveFluid] Thermodynamics model validated..." << endl;

    if (thermo.he().name()=="h")
    {
        useEnthalpy_ = true;
    }
    else if (thermo.he().name()=="e")
    {
        useEnthalpy_ = false;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown energy type: " << thermo.he().name()
            << exit(FatalError);
    }
    
    Info<< "[waveFluid] Energy type: " << (useEnthalpy_ ? "Enthalpy" : "Internal Energy") << endl;
    
    //
    //  Check momentum transport model  
    //
    if (momentumTransport.valid())
    {
        momentumTransport->validate();
        mesh.schemes().setFluxRequired(U.name());
    }
    
    Info<< "[waveFluid] momentum transport model validated..." << endl;

    //
    //  Initialize cantera
    //
    initializeCantera();

    Info<< "[waveFluid] Cantera initialized..." << endl;
    
    //
    //   Set up initial calculation
    //
    
    // Create initial flux set
    fluxPredictor();
    
    // Calculate the maximum wave speed
    if (transient())
    {
        const surfaceScalarField amaxSf
        (
            fluxScheme_.maxWavespeed()*mesh.magSf()
        );

        correctCoNum(amaxSf);
    }
    else if (LTS)
    {
        Info<< "Using LTS" << endl;

        trDeltaT = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    runTime.name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }
    
    Info<< "[waveFluid] fluxes initialized..." << endl;

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::waveFluid::~waveFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::waveFluid::preSolve()
{
    {
        const surfaceScalarField amaxSf
        (
            fluxScheme_.maxWavespeed()*mesh.magSf()
        );

        if (transient())
        {
            correctCoNum(amaxSf);
        }
        else if (LTS)
        {
            setRDeltaT(amaxSf);
        }
    }

    fvModels().preUpdateMesh();

    if (mesh.topoChanging() || mesh.stitcher().stitches())
    {
        pos.clear();
        neg.clear();

        clearTemporaryFields();

        fluxScheme_.topoChange();
    }

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
}

void Foam::solvers::waveFluid::postCorrector()
{
    if (!inviscid && pimple.correctTransport())
    {
        momentumTransport->correct();
        thermophysicalTransport->correct();
    }

    fluxScheme_.corrector();
}

void Foam::solvers::waveFluid::postSolve()
{}

// * * * * * * * * * * * * * * Cantera Functions * * * * * * * * * * * * * * //

using namespace Cantera;

void Foam::solvers::waveFluid::initializeCantera()
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
