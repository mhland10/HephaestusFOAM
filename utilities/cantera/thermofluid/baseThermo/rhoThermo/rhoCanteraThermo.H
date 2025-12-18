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
Control
    Aspects of this software may be subject to CUI//EXPT restrictions. 
    Ensure you are using the public version, or confirm you are a U.S. 
    citizen or otherwise authorized under 22 CFR 120.62.

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
    Foam::rhoCanteraThermo

Description
    Base-class for fluid thermodynamic properties based on density using the
    Cantera implementation.

See also
    Foam::basicThermo

SourceFiles
    rhoCanteraThermo.C

\*---------------------------------------------------------------------------*/

#ifndef rhoCanteraThermo_H
#define rhoCanteraThermo_H

#include "RhoCanteraThermo.H"
#include "pureThermo.H"
#include "rhoThermo.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class rhoCanteraThermo Declaration
\*---------------------------------------------------------------------------*/

class rhoCanteraThermo
:
    virtual public rhoThermo,
    virtual public fluidThermo
{
public:

    // Public Classes

        //- Forward declare the composite class
        class composite;


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType =
            RhoCanteraThermo<BasicThermo<MixtureType, composite>>;

        //- The derived name
        static const word derivedThermoName;


    //- Runtime type information
    TypeName("rhoCanteraThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoCanteraThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<rhoCanteraThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~rhoCanteraThermo();


    // Member Functions

        // Derived Thermodynamic Properties

            //- Rename and return the thermodynamic density field [kg/m^3]
            //  This is used by solvers which create a separate continuity rho
            virtual tmp<volScalarField> renameRho();

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho);
};


/*---------------------------------------------------------------------------*\
                    Class rhoCanteraThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class rhoCanteraThermo::composite
:
    public basicThermo::implementation,
    public pureThermo,
    public rhoThermo::implementation,
    public fluidThermo::implementation,
    public rhoCanteraThermo
{
public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        template<class MixtureType>
        composite
        (
            const dictionary& dict,
            const MixtureType& mixture,
            const fvMesh& mesh,
            const word& phaseName
        )
        :
            basicThermo::implementation(dict, mesh, phaseName),
            rhoThermo::implementation(dict, mesh, phaseName),
            fluidThermo::implementation(dict, mesh, phaseName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //