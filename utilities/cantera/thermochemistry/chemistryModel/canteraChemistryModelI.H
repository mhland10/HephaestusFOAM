/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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

#include "UniformField.H"
#include "localEulerDdtScheme.H"
#include "cpuLoad.H"

// Cantera headers
#include "cantera/zerodim.h"
#include "cantera/thermo/IdealGasPhase.h" // defines class IdealGasPhase
#include "cantera/base/Solution.h"
#include "cantera/transport.h" // transport properties
#include <cantera/kinetics/Kinetics.h>


#include <iostream>

// * * * * * * * * * * * * * * * Inline Functions  * * * * * * * * * * * * * //

template<class ThermoType>
inline const Foam::multicomponentMixture<ThermoType>&
Foam::canteraChemistryModel<ThermoType>::mixture() const
{
    return mixture_;
}

template<class ThermoType>
inline const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::canteraChemistryModel<ThermoType>::RR() const
{
    return RR_;
}

// * * * * * * * * * * * * * * Pull Cantera Namespace  * * * * * * * * * * * //

using namespace Cantera;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::canteraChemistryModel<ThermoType>::canteraChemistryModel
(
    const fluidMulticomponentThermo& thermo
)
:
    chemistryModel<ThermoType>(thermo),
    cpuLoad_(this->lookupOrDefault("cpuLoad", false)),
    mixture_
    (
        dynamicCast<const multicomponentMixture<ThermoType>>(this->thermo())
    ),
    specieThermos_(mixture_.specieThermos()),
    Yvf_(this->thermo().Y()),
    nSpecie_(Yvf_.size()),
    RR_(nSpecie_),
    Y_(nSpecie_),
    c_(nSpecie_)
{
    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Yvf_[fieldi].name(),
                    this->mesh().time().name(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo.mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    Info<< "chemistryModel: Number of species = " << nSpecie_;
    
    const fvMesh& mesh = this->mesh();

    // Look up deltaTChem from the dictionary
    if (mesh.foundObject<dictionary>("chemistryProperties"))
    {
        const dictionary& chemDict = mesh.lookupObject<dictionary>("chemistryProperties");

        if (chemDict.found("deltaTChem"))
        {
            deltaTChem_ = scalarField(mesh.nCells(), chemDict.lookup("deltaTChem"));
        }
        else
        {
            FatalErrorInFunction
                << "deltaTChem must be specified in chemistryProperties"
                << exit(FatalError);
        }

        if (chemDict.found("deltaTChemMax"))
        {
            deltaTChemMax_ = chemDict.lookup("deltaTChemMax");
        }
        else
        {
            deltaTChemMax_ = mesh.time().deltaTValue(); // default to CFD timestep
        }
    }
    else
    {
        FatalErrorInFunction
            << "chemistryProperties dictionary not found in mesh"
            << exit(FatalError);
    }
    
    // ?? FIX: Use this->mesh() instead of bare mesh
    const dictionary& runDict = mesh.template lookupObject<dictionary>("chemistryProperties");
    const dictionary& canteraDict = runDict.subDict("cantera");

    // Read some properties
    const word yamlFile = canteraDict.lookupOrDefault<word>("dataFile", "gri30.yaml");
    const word phase    = canteraDict.lookupOrDefault<word>("phaseName", "gri30");
    Info << "Cantera YAML file: " << yamlFile << ", phase: " << phase << nl;
    
    // Create the solution with default (ambient) state3
    solution_ = newSolution(yamlFile, phase);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::canteraChemistryModel<ThermoType>::~canteraChemistryModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::canteraChemistryModel<ThermoType>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();

    const volScalarField& Tvf = this->thermo().T();
    const volScalarField& pvf = this->thermo().p();
    
    // Scratch space for mass fractions
    scalarField Yi(nSpecie_);
    
    // Set a thermo object for each cell to use
    auto thermo = solution_->thermo();
    
    // Allocate vector of net production rates
    std::vector<double> wdot(nSpecie_);

    forAll(rhovf, celli)
    {
        const scalar rho = rhovf[celli];
        const scalar T = Tvf[celli];
        const scalar p = pvf[celli];

        // Pull composition from OpenFOAM fields
        for (label i = 0; i < nSpecie_; i++)
        {
            Yi[i] = Yvf_[i][celli];
        }

        // Convert to Cantera composition string
        std::string canteraComposition = printCanteraMixture(Yi, mixture_);

        // Set the Cantera thermo state for this cell
        thermo->setState_TPY(T, p, canteraComposition);        

        // Pull back species production rates
        solution_->kinetics()->getNetProductionRates(wdot.data());
        for (label i = 0; i < nSpecie_; i++)
        {
            RR_[i][celli] = wdot[i] * specieThermos_[i].W(); // dY/dt
        }
        
    }
}

template<class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::canteraChemistryModel<ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    optionalCpuLoad& chemistryCpuLoad
    (
        optionalCpuLoad::New("canteraChemistryModel:cpuLoad", this->mesh(), cpuLoad_)
    );

    // CPU time logging
    cpuTime solveCpuTime;
    scalar totalSolveCpuTime = 0;

    if (!this->chemistry_)
    {
        return great;
    }

    const volScalarField& rho0vf =
        this->mesh().template lookupObject<volScalarField>
        (
            this->thermo().phasePropertyName("rho")
        ).oldTime();

    const volScalarField& T0vf = this->thermo().T().oldTime();
    const volScalarField& p0vf = this->thermo().p().oldTime();

    scalarField Y0(nSpecie_);

    // Minimum chemical timestep
    scalar deltaTMin = great;

    chemistryCpuLoad.resetCpuTime();

    // Set the Cantera thermo initialization
    auto& thermo = solution_->thermo();

    // Create a temporary reactor for this cell
    IdealGasReactor reactor(thermo);
    ReactorNet net;
    net.addReactor(reactor);

    // Optionally configure tolerances (CVODE/SUNDIALS)
    net.setTolerances(1e-8, 1e-15); // relative, then absolute
    
    forAll(rho0vf, celli)
    {
        const scalar rho0 = rho0vf[celli];

        scalar p = p0vf[celli];
        scalar T = T0vf[celli];
        
        for (label i=0; i<nSpecie_; i++)
        {
            Y_[i] = Y0[i] = Yvf_[i].oldTime()[celli];
        }

        // Convert to Cantera composition string
        std::string canteraComposition = printCanteraMixture(Y_, mixture_);

        // Set the Cantera thermo state for this cell
        thermo->setState_TPY(T, p, canteraComposition);
        
        // Initialise time progress
        scalar timeLeft = deltaT[celli];
        
        if (log_)
        {
            // Reset the solve time
            solveCpuTime.cpuTimeIncrement();
        }

        // Calculate the chemical source terms
        while (timeLeft > small)
        {
            scalar dt = min(timeLeft, deltaTChem_[celli]);
            
            // Integrate with Cantera for dt seconds
            net.advance(dt);             // advance the reactor network
            
            // Retrieve updated state
            T = reactor.temperature();
            for (label i=0; i<nSpecie_; i++)
                Y_[i] = reactor.massFraction(i);
            
            // Optionally update pressure if using a variable-P integration
            p = reactor.pressure();

            timeLeft -= dt;
        }

        if (log_)
        {
            totalSolveCpuTime += solveCpuTime.cpuTimeIncrement();
        }

        deltaTMin = min(deltaTChem_[celli], deltaTMin);
        deltaTChem_[celli] = min(deltaTChem_[celli], deltaTChemMax_);
        

        // Set the RR vector (used in the solver)
        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = rho0*(Y_[i] - Y0[i])/deltaT[celli];
        }

        if (cpuLoad_)
        {
            chemistryCpuLoad.cpuTimeIncrement(celli);
        }
    }

    if (log_)
    {
        cpuSolveFile_()
            << this->time().userTimeValue()
            << "    " << totalSolveCpuTime << endl;
    }

    return deltaTMin;
}


template<class ThermoType>
Foam::scalar Foam::canteraChemistryModel<ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class ThermoType>
Foam::scalar Foam::canteraChemistryModel<ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::canteraChemistryModel<ThermoType>::tc() const
{
    tmp<volScalarField> ttc
    (
        volScalarField::New
        (
            "tc",
            this->mesh(),
            dimensionedScalar(dimTime, small),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    scalarField& tc = ttc.ref();

    if (!this->chemistry_)
    {
        ttc.ref().correctBoundaryConditions();
        return ttc;
    }

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();

    const volScalarField& Tvf = this->thermo().T();
    const volScalarField& pvf = this->thermo().p();
    
    const scalar tiny = 1e-30; // threshold to avoid division by zero

    forAll(rhovf, celli)
    {
        const scalar rho = rhovf[celli];
        const scalar T = Tvf[celli];
        const scalar p = pvf[celli];
        
        for (label i=0; i<nSpecie_; i++)
        {
            Y_[i] = Yvf_[i][celli];
        }
        
        scalar tc_cell = 1/tiny;
        for (label i = 0; i < nSpecie_; i++)
        {
            if (Y_[i] > tiny && std::abs(RR_[i][celli]) > tiny)
            {
                scalar tc_i = Y_[i] / std::abs(RR_[i][celli]);
                tc_cell = std::min(tc_cell, tc_i);
            }
        }

        tc[celli] = tc_cell;
    }

    ttc.ref().correctBoundaryConditions();
    return ttc;
}


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::canteraChemistryModel<ThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (!this->chemistry_)
    {
        return tQdot;
    }

    //Cantera::reactionEvaluationScope scope(*this);
    
    auto& thermo = solution_->thermo();

    scalarField& Qdot = tQdot.ref();
    
    const auto& names = mixture_.specieNames();
    const double T_ref = 298.15; // Kelvin
    const double P_ref = 101325.0; // Pa

    forAll(Yvf_, i)
    {
        std::ostringstream streamHold;
        streamHold << names[i] << ":1.0";
        std::string canteraComposition = streamHold.str();
        thermo->setState_TPY(T_ref, P_ref, canteraComposition);
        const scalar hi = thermo->enthalpy_mass(); // J/kg
        
        forAll(Qdot, celli)
        {
            Qdot[celli] -= hi*RR_[i][celli];
        }
    }

    return tQdot;
}


// ************************************************************************* //
