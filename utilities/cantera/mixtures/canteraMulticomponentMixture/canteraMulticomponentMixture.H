/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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
    Foam::canteraMulticomponentMixture

Description
    Multi component mixture to use Cantera

SourceFiles
    canteraMulticomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef canteraMulticomponentMixture_H
#define canteraMulticomponentMixture_H

#include "wordList.H"
#include "FieldListSlice.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class canteraMulticomponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class canteraMulticomponentMixture
{
public:

    // Public Typedefs
    
        

        //- The type of thermodynamics this mixture is instantiated for
        typedef ThermoType thermoType;

        //- Mixing type for thermodynamic properties
        typedef ThermoType thermoMixtureType;

        //- Mixing type for transport properties
        typedef ThermoType transportMixtureType;


/*
private:

    // Private Data

        //- Thermo model
        mutable autoPtr<ThermoType> mixture_;
*/

public:

    mutable thermoType mixture_;

    // Constructors

        //- Construct from a dictionary
        canteraMulticomponentMixture(const dictionary&);


    //- Destructor
    virtual ~canteraMulticomponentMixture();


    // Member Functions

        //- Return the instantiated type name
        //  Note: This defines the same name as pureMixture so that
        //  non-reacting and reacting solvers are compatible
        static word typeName()
        {
            return "multicomponentMixture<" + ThermoType::typeName + ">";
        }


        //- Return the specie names
        wordList specieNames() const
        {
            return wordList();
        }

        //- Return the mixture for thermodynamic properties
        inline const thermoMixtureType& thermoMixture
        (
            const scalarFieldListSlice&
        ) const;

        //- Return the mixture for transport properties
        inline const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&
        ) const;

        //- Return the mixture for transport properties
        inline const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&,
            const thermoMixtureType&
        ) const;

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        inline const ThermoType& specieThermo(const label speciei) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "canteraMulticomponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
