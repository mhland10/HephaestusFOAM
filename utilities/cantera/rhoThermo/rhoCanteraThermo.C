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

\*---------------------------------------------------------------------------*/

#include "rhoCanteraThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rhoCanteraThermo, 0);
    defineRunTimeSelectionTable(rhoCanteraThermo, fvMesh);
}

const Foam::word Foam::rhoCanteraThermo::derivedThermoName("heRhoCanteraThermo");


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::rhoCanteraThermo> Foam::rhoCanteraThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<rhoCanteraThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rhoCanteraThermo::~rhoCanteraThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::rhoCanteraThermo::renameRho()
{
    rho().rename(phasePropertyName(Foam::typedName<rhoCanteraThermo>("rho")));
    return rho();
}


void Foam::rhoCanteraThermo::correctRho(const volScalarField& deltaRho)
{
    rho() += deltaRho;
}


// ************************************************************************* //
