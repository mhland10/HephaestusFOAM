/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "psiMulticomponentCanteraThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(psiMulticomponentCanteraThermo, 0);
    defineRunTimeSelectionTable(psiMulticomponentCanteraThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::psiMulticomponentCanteraThermo> Foam::psiMulticomponentCanteraThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<psiMulticomponentCanteraThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::psiMulticomponentCanteraThermo::~psiMulticomponentCanteraThermo()
{}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void Foam::psiMulticomponentCanteraThermo::calculate()
{
  Foam::Info << "Using calculate from composite definition" << Foam::nl;
  
  // Forward the call to the composite object
  //composite& comp = static_cast<composite&>(*this);
  
  /*
  //  Create an alias for the energy and pressure in the cells. This
  // this works by receiving the field (OpenFOAM) through the gofer 
  // function and points to it in memory.
  const volScalarField& hCells = this->he();
  const volScalarField& pCells = this->p();
  
  //  Create an alias for the temperature, flux dialation (psi), 
  // viscosity, and thermal conductivity. This is simply able to point
  // to the field that the composite object inherits.
  scalarField& TCells  = this->T_.primitiveFieldRef();
  scalarField& psiCells = this->psi_.primitiveFieldRef();
  scalarField& muCells = this->mu_.primitiveFieldRef();
  scalarField& kappaCells = this->kappa_.primitiveFieldRef();
  
  //  Create an alias for the constant pressure and volume specific 
  // heat. Like with energy, the field needs to be pulled from the 
  // gofer functions, but we use an auto pointer because apparently I 
  // have no clue how to properly argument this function.
  auto& CpCells  = this->Cp();
  auto& CvCells  = this->Cv();
  
  //  Pull a copy of the Yslicer in the object that was inherited.
  auto Yslicer = this->Yslicer();
  
  //  Now we loop through the cells to apply our data to the fluid in
  // the cells.
  forAll(TCells, celli)
  {
      auto composition = this->cellComposition(Yslicer, celli);
      
      
  }
  */


}

// ************************************************************************* //
