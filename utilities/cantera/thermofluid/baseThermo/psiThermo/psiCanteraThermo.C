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

#include "psiCanteraThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(psiCanteraThermo, 0);
    defineRunTimeSelectionTable(psiCanteraThermo, fvMesh);
}

const Foam::word Foam::psiCanteraThermo::derivedThermoName("hePsiCanteraThermo");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::psiCanteraThermo::implementation::implementation
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& phaseName
)
{
    Foam::Info << "psiCanteraThermo constructor using dictionary, mesh, and phaseName" << Foam::nl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::psiCanteraThermo> Foam::psiCanteraThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<psiCanteraThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::psiCanteraThermo::~psiCanteraThermo()
{}


Foam::psiCanteraThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::psiCanteraThermo::renameRho()
{
    return rho();
}


void Foam::psiCanteraThermo::correctRho(const Foam::volScalarField& deltaRho)
{}


Foam::tmp<Foam::volScalarField> Foam::psiCanteraThermo::implementation::rho() const
{
    return p()*psi();
}


Foam::tmp<Foam::scalarField> Foam::psiCanteraThermo::implementation::rho
(
    const label patchi
) const
{
    return p().boundaryField()[patchi]*psi().boundaryField()[patchi];
}


// ************************************************************************* //
