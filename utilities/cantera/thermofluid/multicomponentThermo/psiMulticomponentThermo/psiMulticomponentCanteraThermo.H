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
    Foam::psiMulticomponentCanteraThermo

Description
    Base-class for multi-component fluid thermodynamic properties based on
    compressibility that uses Cantera for its engine.

See also
    Foam::basicThermo

SourceFiles
    psiMulticomponentCanteraThermo.C

\*---------------------------------------------------------------------------*/

#ifndef psiMulticomponentCanteraThermo_H
#define psiMulticomponentCanteraThermo_H

#include "psiCanteraThermo.H"
#include "fluidMulticomponentThermo.H"
#include "canteraFluidMulticomponentThermo.H"
#include "BasicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class psiMulticomponentCanteraThermo Declaration
\*---------------------------------------------------------------------------*/

class psiMulticomponentCanteraThermo
:
    virtual public psiCanteraThermo,
    virtual public fluidMulticomponentThermo
{
public:

    // Public Classes

        //- Forward declare the composite class
        class composite;
        
    // Public Methods
        virtual void calculate();// override;


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType =
            FluidMulticomponentThermo
            <
                MulticomponentThermo
                <
                    PsiCanteraThermo
                    <
                        BasicThermo<MixtureType, composite>
                    >
                >
            >;


    //- Runtime type information
    TypeName("psiMulticomponentCanteraThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiMulticomponentCanteraThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<psiMulticomponentCanteraThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~psiMulticomponentCanteraThermo();
};


/*---------------------------------------------------------------------------*\
          Class psiMulticomponentCanteraThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class psiMulticomponentCanteraThermo::composite
:
    public basicThermo::implementation,
    public fluidThermo::implementation,
    public psiCanteraThermo::implementation,
    public multicomponentThermo::implementation,
    public psiMulticomponentCanteraThermo
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
            fluidThermo::implementation(dict, mesh, phaseName),
            psiCanteraThermo::implementation(dict, mesh, phaseName),
            multicomponentThermo::implementation
            (
                dict,
                mixture.specieNames(),
                mesh,
                phaseName
            )
            
            
        {
          Info << "Created psiMulticomponentCanteraThermo" << nl;
        }
        
        template<class MixtureType>
        void calculate()
        {
            Info << "Using multicomponent calculate with species" << nl;
        }  
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
