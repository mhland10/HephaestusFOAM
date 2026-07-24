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
    Foam::fluxScheme

Description
    This class is a utility for calculating the fluxes that go into the div()
    operator in the various conservation equations. This also includes the 
    construction and utilities needed for robust flux evaluation across 
    multiple schemes

SourceFiles
    fluxScheme.C

See also
    Foam::solvers::waveHephaestus

\*---------------------------------------------------------------------------*/

#ifndef fluxScheme_H
#define fluxScheme_H

#include "fluidSolver.H"
#include "psiThermo.H"
#include "compressibleMomentumTransportModel.H"
#include "fluidThermoThermophysicalTransportModel.H"

// Include the following for multicomponent capability
#include "fluidMulticomponentThermo.H"
#include "combustionModel.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"
#include "multivariateScheme.H"

// Cantera headers
#include "cantera/thermo/IdealGasPhase.h" // defines class IdealGasPhase
#include "cantera/base/Solution.h"
#include "cantera/transport.h" // transport properties
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class fluxScheme Declaration
\*---------------------------------------------------------------------------*/

class fluxScheme
{
    //=====================================================================
    //
    //  Object Data
    //
    //=====================================================================

    protected:

    // Dictionaries & controls - solver storage

        //- Reference to the mesh
        const fvMesh& mesh;

        //- Reference to the runtime
        const Time& runTime_;

    // Thermo models

        //- Reference to the fluid thermophysical properties
        fluidMulticomponentThermo& thermo_; 

    // Thermophysical properties
    
        //- Velocity field
        volVectorField& U_;
        
        //- The continuity density field
        volScalarField& rho_;

        //- Reference to the pressure field
        volScalarField& p_;

        //- Reference to the temperature field
        volScalarField& T_;

        //- Reference to the static energy field
        volScalarField& he_;
        
        //- Reference to the species fraction fields
        PtrList<volScalarField>& Y_;

     // Derived properties

        //- Mass-flux field
        surfaceScalarField phi_;

        //- Kinetic energy field
        //  Used in the energy equation
        volScalarField K_;

    // Cantera data 

        // Create pointer for Cantera solution
        Cantera::shared_ptr<Cantera::Solution> canteraSolution_;

        // Allow easy switching between enthalpy and internal energy
        bool useEnthalpy_;

    // Controls

        //- Which flux scheme the object will use
        word fluxSchemeSelection_;

    // Thermophysical models
    
        IOdictionary physicalProperties_;

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
        
        tmp<surfaceScalarField> T_pos;
        tmp<surfaceScalarField> T_neg;

        tmp<surfaceScalarField> he_pos;
        tmp<surfaceScalarField> he_neg;

        tmp<surfaceScalarField> a_pos;
        tmp<surfaceScalarField> a_neg;

        tmp<surfaceScalarField> c_pos;
        tmp<surfaceScalarField> c_neg;

        tmp<surfaceScalarField> aSf;

        tmp<surfaceScalarField> aphiv_pos;
        tmp<surfaceScalarField> aphiv_neg;

        tmp<surfaceVectorField> devTau;
        
        tmp<surfaceScalarField> rhoHE_pos;
        tmp<surfaceScalarField> rhoHE_neg;
        
        List<tmp<surfaceScalarField>> rhoYi_pos;
        List<tmp<surfaceScalarField>> rhoYi_neg;
        
        List<tmp<surfaceScalarField>> Yi_pos;
        List<tmp<surfaceScalarField>> Yi_neg;

        //- Optional LTS reciprocal time-step field
        tmp<volScalarField> trDeltaT;

        //- The flux fields
        tmp<surfaceVectorField> phiUp;
        tmp<surfaceScalarField> phiHEp;
        tmp<surfaceScalarField> phiP;
        List<tmp<surfaceScalarField>> phiYi;

    //- Flux terms to deliver to solver

        //- Mass flux over face area
        surfaceScalarField phiRho_;

        //- Momentum flux over face area
        surfaceVectorField phiRhoUplusP_;

        //- Energy flux over face area
        surfaceScalarField phiRhoHE_;

        //- Species fluxes over face area
        PtrList<surfaceScalarField> phiRhoYi_;

    //=====================================================================
    //
    //  Functions
    //
    //=====================================================================

    private:

    // Private Member Functions

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

        

        void clearTemporaryFields();

    public:

    // Public Functions

        //- Constructors
        fluxScheme
        (
            const Time& runTime,
            const fvMesh& mesh,
            fluidMulticomponentThermo& thermo,
            volScalarField& p,
            volVectorField& U,
            volScalarField& rho,
            PtrList<volScalarField>& Y
        );

        //- Disallow default bitwise copy construction
        fluxScheme(const fluxScheme&) = delete;

        //- Destructor
        virtual ~fluxScheme();

        //- Disallow default bitwise assignment
        void operator=(const fluxScheme&) = delete;

        //- The function that provides the predictor-stage calculation of the 
        //  fluxes
        void predictor();

        //- The function that provides the corrector-stage operations
        void corrector(); 

    // Accessing Functions

        //- Return the fluxes
        const surfaceScalarField& phiRho() const { return phiRho_; }
        const surfaceVectorField& phiRhoUplusP() const { return phiRhoUplusP_; }
        const surfaceScalarField& phiRhoHE() const { return phiRhoHE_; }
        const PtrList<surfaceScalarField>& phiRhoYi() const { return phiRhoYi_; }
        const surfaceScalarField& phiRhoY(const label i) const{ return phiRhoYi_[i]; }

    protected:

    // Cantera functions

        // Cantera initialization function
        void initializeCantera();        

        //- Function for face evaluation at the solver level
        virtual void faceEvaluate
        (
            const word& property,
            const surfaceScalarField& pFace,
            const surfaceScalarField& TFace,
            const PtrList<surfaceScalarField>& YFace,
            surfaceScalarField& field
        ) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}; // End class fluxScheme
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //