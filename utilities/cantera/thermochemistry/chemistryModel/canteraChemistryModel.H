/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

Class
    Foam::canteraChemistryModel

Description
    This package creates a chemistry model that integrates Cantera into 
    OpenFOAM. This in turn puts much of the solving operation within Cantera's
    Sundials system.

SourceFiles
    canteraChemistryModelI.H
    canteraChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef canteraChemistryModel_H
#define canteraChemistryModel_H

#include "chemistryModel.H"

#include "volFields.H"
#include "multicomponentMixture.H"
#include "DynamicField.H"

#include "cantera/thermo/IdealGasPhase.h" // defines class IdealGasPhase
#include "cantera/base/Solution.h"
#include "cantera/transport.h" // transport properties
#include <iostream>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                  Class canteraChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class canteraChemistryModel
:
    public chemistryModel<ThermoType>
{

    // Private data
    
        //- Reference to the field of specie mass fractions
        const PtrList<volScalarField>& Yvf_;
    
        //- Switch to select performance logging
        Switch log_;

        //- Switch to enable per-cell CPU load caching for load-balancing
        Switch cpuLoad_;

        //- Type of the Jacobian to be calculated
        //const jacobianType jacobianType_;

        //- Reference to the multi component mixture
        const multicomponentMixture<ThermoType>& mixture_;
        
        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermos_;

        //- List of reaction rate per specie [kg/m^3/s]
        PtrList<volScalarField::Internal> RR_;

        //- Temporary mass fraction field
        mutable scalarField Y_;
        
        //- Temporary concentration field
        mutable scalarField c_;
        
        //- Number of species
        label nSpecie_;
        
        //- Log file for average time spent solving the chemistry
        autoPtr<OFstream> cpuSolveFile_;
        
        // Per-cell chemical timestep
        scalarField deltaTChem_;
    
        // Maximum allowed chemical timestep (optional)
        scalar deltaTChemMax_;

    // Private Member Functions

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        //  Variable number of species added
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);

protected:

    // Create pointer for Cantera solution
    Cantera::shared_ptr<Cantera::Solution> solution_;

public:

    //- Runtime type information
    TypeName("canteraChemistryModel");


    // Constructors

        //- Construct from thermo
        canteraChemistryModel(const fluidMulticomponentThermo& thermo);

        //- Disallow default bitwise copy construction
        canteraChemistryModel(const canteraChemistryModel&) = delete;


    //- Destructor
    virtual ~canteraChemistryModel();


    // Member Functions

        // Access

            //- Return reference to the mixture
            inline const multicomponentMixture<ThermoType>& mixture() const;
            
            /*
            //- The reactions
            inline const PtrList<Reaction<ThermoType>>& reactions() const;

            //- Thermodynamic data of the species
            inline const PtrList<ThermoType>& specieThermos() const;
            */
            //- The number of species
            inline virtual label nSpecie() const;


        // Overrides to basicChemistryModel functions

            //- Return reaction rates of the species [kg/m^3/s]
            virtual inline const PtrList<volScalarField::Internal>& RR() const;

            //- Calculates the reaction rates
            virtual void calculate() override;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return the heat release rate [kg/m/s^3]
            virtual tmp<volScalarField> Qdot() const;

        // ODE solution functions

            //- Solve the ODE system
            virtual void solve
            (
                scalar& p,
                scalar& T,
                scalarField& Y,
                const label li,
                scalar& deltaT,
                scalar& subDeltaT
            ) const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const canteraChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "canteraChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "canteraChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
