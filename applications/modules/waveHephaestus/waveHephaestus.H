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
    Foam::solvers::waveFluid

Description
    Solver module for density-based solution of compressible flow that
    incorporates the multicomponent gas infrastructure, unlike the module it
    is derived from, `shockFluid`.

    For now, its flux schemes are based on central-upwind schemes of Kurganov 
    and Tadmor with support for mesh-motion and topology change.

    Reference:
    \verbatim
        Greenshields, C. J., Weller, H. G., Gasparini, L.,
        & Reese, J. M. (2010).
        Implementation of semi‐discrete, non‐staggered central schemes
        in a colocated, polyhedral, finite volume framework,
        for high‐speed viscous flows.
        International journal for numerical methods in fluids, 63(1), 1-21.
    \endverbatim

SourceFiles
    waveHephaestus.C

See also
    Foam::solvers::fluidSolver
    Foam::solvers::incompressibleFluid

\*---------------------------------------------------------------------------*/

#ifndef waveHephaestus_H
#define waveHephaestus_H

#include "fluidSolver.H"
#include "psiThermo.H"
#include "compressibleMomentumTransportModel.H"
#include "fluidThermoThermophysicalTransportModel.H"

// Include the following for multicomponent capability
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
                          Class waveFluid Declaration
\*---------------------------------------------------------------------------*/

class waveFluid
:
    public fluidSolver
{

protected:

    // Thermophysical properties

        //- Pointer to the fluid thermophysical properties
        autoPtr<fluidMulticomponentThermo> thermoPtr_;

        //- Reference to the fluid thermophysical properties
        fluidMulticomponentThermo& thermo_;

        //- Reference to the pressure field
        volScalarField& p_;

        //- The continuity density field
        volScalarField rho_;
        
        //- Reference to the species fraction fields
        PtrList<volScalarField>& Y_;
        
        //- The species concentrations
        PtrList<volScalarField> rhoYi_;
        
        //- Pointer to the combustion model  
        autoPtr<combustionModel> reaction;
        
        //- The multivariate interpolation scheme
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
        /*
          For those of us wondering what this is, it is the way we 
        couple the interpolation between different variables so we're 
        not using different methods for different variables that could
        result in bizarre distributions of nonsenese. For example, if 
        we had a linear interpolation of density, and a van Leer 
        reconstruction of velocity, it could greatly throw off other
        models trying to interact with the data.
        
        */
        
        //- The thermophysical transport model
        autoPtr<fluidMulticomponentThermophysicalTransportModel>
            thermophysicalTransport;


    // Kinematic properties

        //- Velocity field
        volVectorField U_;

        //- Mass-flux field
        surfaceScalarField phi_;

        //- Kinetic energy field
        //  Used in the energy equation
        volScalarField K;


    // Momentum transport

        bool inviscid;

        //- Pointer to the momentum transport model
        autoPtr<compressible::momentumTransportModel> momentumTransport;
        /*
          TODO: at some point we should generate some helper modules
        for this variable.
        
        */

    // Controls

        word fluxScheme;
        /*
          TODO: in the future when the flux scheme is broken into its
        own class, this will end up being a pointer to that class.
        
        */


    // Cached temporary fields

        tmp<surfaceScalarField> pos;
        tmp<surfaceScalarField> neg;

        tmp<surfaceScalarField> rho_pos;
        tmp<surfaceScalarField> rho_neg;

        tmp<surfaceVectorField> rhoU_pos;
        tmp<surfaceVectorField> rhoU_neg;

        tmp<surfaceVectorField> U_pos;
        tmp<surfaceVectorField> U_neg;

        tmp<surfaceScalarField> p_pos;
        tmp<surfaceScalarField> p_neg;

        tmp<surfaceScalarField> a_pos;
        tmp<surfaceScalarField> a_neg;

        tmp<surfaceScalarField> aSf;

        tmp<surfaceScalarField> aphiv_pos;
        tmp<surfaceScalarField> aphiv_neg;

        tmp<surfaceVectorField> devTau;
        
        tmp<surfaceScalarField> rhoHE_pos;
        tmp<surfaceScalarField> rhoHE_neg;
        
        tmp<surfaceScalarField> HE_pos;
        tmp<surfaceScalarField> HE_neg;
        
        List<tmp<surfaceScalarField>> rhoYi_pos;
        List<tmp<surfaceScalarField>> rhoYi_neg;
        
        List<tmp<surfaceScalarField>> Yi_pos;
        List<tmp<surfaceScalarField>> Yi_neg;

        //- Optional LTS reciprocal time-step field
        tmp<volScalarField> trDeltaT;


private:

    // Private Member Functions

        //- Set rDeltaT for LTS
        void setRDeltaT(const surfaceScalarField& amaxSf);

        //- Correct the cached Courant numbers
        void correctCoNum(const surfaceScalarField& amaxSf);

        //- Construct the continuity equation and correct the density
        void correctDensity();

        //- Interpolate field vf according to direction dir
        template<class Type>
        inline tmp<SurfaceField<Type>> interpolate
        (
            const VolField<Type>& vf,
            const surfaceScalarField& dir,
            const word& reconFieldName = word::null
        )
        {
            tmp<SurfaceField<Type>> tsf
            (
                fvc::interpolate
                (
                    vf,
                    dir,
                    "reconstruct("
                  + (reconFieldName != word::null ? reconFieldName : vf.name())
                  + ')'
                )
            );

            SurfaceField<Type>& sf = tsf.ref();

            sf.rename(vf.name() + '_' + dir.name());

            return tsf;
        }

        void fluxPredictor();
        /*
          TODO: Alter this to handle the additional flux terms
          
        */

        void clearTemporaryFields();


public:

    // Public Data

        //- Reference to the fluid thermophysical properties
        const fluidMulticomponentThermo& thermo;

        //- Reference to the pressure field
        const volScalarField& p;

        //- Reference to the continuity density field
        const volScalarField& rho;

        //- Reference to the velocity field
        const volVectorField& U;

        //- Reference to the mass-flux field
        const surfaceScalarField& phi;
        
        //- Reference to the composition
        const PtrList<volScalarField>& Y;
        
        //- Reference to total energy density (optional for public API)
        const volScalarField& he;


    //- Runtime type information
    TypeName("waveHephaestus");


    // Constructors

        //- Construct from region mesh
        waveFluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        waveFluid(const waveFluid&) = delete;


    //- Destructor
    virtual ~waveFluid();


    // Member Functions

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual void moveMesh();

        //- Corrections that follow mesh motion
        virtual void motionCorrector();

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum and thermophysical transport modelling
        virtual void postCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const waveFluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
