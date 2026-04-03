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

Class
    Foam::solvers::multicomponentFluid

Description
    Solver module for steady or transient turbulent flow of compressible
    multicomponent fluids with optional mesh motion and change.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, chemical reactions,
    combustion, Lagrangian particles, radiation, surface film etc. and
    constraining or limiting the solution.

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

SourceFiles
    gasHephaestus.C

See also
    Foam::solvers::fluidSolver
    Foam::solvers::isothermalFluid

\*---------------------------------------------------------------------------*/

#ifndef gasHephaestus_H
#define gasHephaestus_H

#include "isothermalFluid.H"
#include "fluidMulticomponentThermo.H"
#include "combustionModel.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"
#include "multivariateScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class multicomponentFluid Declaration
\*---------------------------------------------------------------------------*/

class multicomponentFluid
:
    public isothermalFluid
{

protected:

    // Thermophysical properties

        fluidMulticomponentThermo& thermo_;


    // Composition

        PtrList<volScalarField>& Y_;


    // Reactions

        autoPtr<combustionModel> reaction;
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;


    // Thermophysical transport

        autoPtr<fluidMulticomponentThermophysicalTransportModel>
            thermophysicalTransport;


private:

    // Private Member Functions

        //- Set rDeltaT for LTS
        virtual void setRDeltaT();


public:

    // Public Data

        //- Reference to the fluid thermophysical properties
        const fluidMulticomponentThermo& thermo;

        //- Reference to the composition
        const PtrList<volScalarField>& Y;
        
        //- The solver state
        struct SolverState
        {
            Field<vector> U;
            Field<scalar> p;
            Field<scalar> he;
            List<Field<scalar>> Yi;
        };
        /*
          This data is important to the higher order time stepping schemes like Runge-Kutta so that
        the solver can output and delta the state in the time stepping or iteration.
        
          We use Foam field objects to lighten the load of what is stored. OpenFOAM fv fields 
        contain a variety of additional information that is either useless or counterproductive to 
        store.
        
        */


    //- Runtime type information
    TypeName("GasHephaestus");


    // Constructors

        //- Construct from region mesh
        multicomponentFluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        multicomponentFluid(const multicomponentFluid&) = delete;


    //- Destructor
    virtual ~multicomponentFluid();


    // Member Functions

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Correct the momentum and thermophysical transport modelling
        virtual void postCorrector();
        
        
        /*
        
          Below are the functions that are required to do higher order time stepping, specifically
        Carpenter-style RK stepping where the PIMPLE loop is disconnected from the time stepping
        scheme. Remember that this stepping follows the formula:
        
        $$
        W^{(s)}=\alpha_sW^{(0)}+\beta_s\left( W^{(s-1)}+\Delta tR(W^{(s-1)} \right)
        $$
        
        */
        
        //- Add in a function to hold the current solver state
        SolverState outputState() const;
        
        //- Add in a function that differences two states
        SolverState deltaState( const SolverState& W0) const;
        
        //- Add in a function that blends two states
        SolverState blendState( const SolverState& W0, const SolverState& dW, scalar alpha, scalar beta ) const;
        
        //- Add in a function to close/restore the state for the time stepping method
        void restoreState(const SolverState& W);
        
        //- Add in a function that flags the solver as RK-compatible
        virtual bool supportsMultiStageTimeIntegration() const
        {
            return true;
        }
         


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multicomponentFluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
