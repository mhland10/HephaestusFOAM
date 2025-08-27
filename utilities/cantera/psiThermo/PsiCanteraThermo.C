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

#include "PsiCanteraThermo.H"
#include "OSspecific.H"

// Cantera headers
#include "cantera/thermo/IdealGasPhase.h" // defines class IdealGasPhase
#include "cantera/base/Solution.h"
//#include "cantera/kinetics/GasKinetics.h"
#include "cantera/transport.h" // transport properties
#include <iostream>

// * * * * * * * Helpful Functions in an Anonymous Namespace * * * * * * * * //

// Put helpers in anonymous namespace so they don�t pollute global scope
namespace
{
    // Helper to detect specieNames()
    template <typename T>
    class hasSpecieNames
    {
    private:
        template <typename U> static auto test(int) -> decltype(std::declval<U>().specieNames(), std::true_type());
        template <typename U> static std::false_type test(...);
    
    public:
        static constexpr bool value = decltype(test<T>(0))::value;
    };
    
    template <typename ContainerType, typename MixtureType>
    typename std::enable_if<hasSpecieNames<MixtureType>::value, std::string>::type
    printCanteraMixture(const ContainerType& composition, const MixtureType& mixture)
    {
        //  Create a stream variable that allows for more efficient storag
        // through the iterations throught the composition
        std::ostringstream streamHold;
        
        //  Store the specie names from the mixture into a pointer for the 
        // names
        const auto& names = mixture.specieNames();
        
        //  Loop through the composition entries
        for (Foam::label i = 0; i < composition.size(); ++i)
        {
            if (composition[i]>=1e-12)
            {
              if (streamHold.tellp() > 0)  // only insert comma if not first entry
              {
                streamHold << ", " ;
              }
              streamHold << names[i] << ":" << composition[i] ;
            }
            
        }
        
        // Alter the string stream to an actual string to pass through
        std::string canteraComposition = streamHold.str();
        
        return canteraComposition;
    }
    
    /*
      This function prints the mixture that Cantera expects into a string format.
      
      **ONLY** if the mixture has no specie names**
      
    Args:
      composition (ref to ContainerType):  This is the output that comes from 
                                          the cellComposition method in the 
                                          object done with the Yslicer, also in
                                          the object.
                                          
      mixture (ref to a OpenFOAM Mixture):  Self explanatory, but this needs to
                                          be the mixture.
                                          
    Returns:
      canteraComposition (str):  The string that Cantera is expecting.
    */
    template <typename ContainerType, typename MixtureType>
    typename std::enable_if<!hasSpecieNames<MixtureType>::value, std::string>::type
    printCanteraMixture(const ContainerType& composition, const MixtureType& mixture)
    {
        // Return a Cantera value of 0
        std::string canteraComposition = "0";
        
        return canteraComposition;
    }

}

using namespace Cantera;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BaseThermo>
void Foam::PsiCanteraThermo<BaseThermo>::calculate()
{
    Info << "Using calculate method from PsiCanteraThermo" << nl;
    
    // Reference the energy and pressure in the cells as a scalarField 
    //  reference. The scalarField is an OpenFOAM type. We hold these constant
    //  to avoid accidentally changing them.
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;
    
    // Reference the other flow parameters in the cells as a scalarField
    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();
    
    // Pull a copy of the Yslicer method in this object that was inherited as
    //  part of the composite.
    auto Yslicer = this->Yslicer();

    forAll(TCells, celli)
    {
        auto composition = this->cellComposition(Yslicer, celli);
            
        // Print the composition
        //Info << "Cell " << celli << " composition: " << nl;
        const auto& mixtureObj = this->mixture();
        std::string canteraComposition = printCanteraMixture(composition, mixtureObj);
        //Info << "    " << canteraComposition << nl;
        
        // Set the state and mixture 
        auto thermo = solution_->thermo();
        thermo->setState_TPY(298.15, pCells[celli], canteraComposition);  
        double h_ref = thermo->enthalpy_mass();   
        thermo->setState_TPY(TCells[celli], pCells[celli], canteraComposition);  
        /*
          Note:  We do need to consider the reference enthalpy or enthalpy of 
                formation because OpenFOAM works with sensible enthalpy, 
                whereas Cantera works with absolute enthalpy. This also applies
                to internal energy, but we go through enthalpy for this.
        */
         
        if (useEnthalpy_)
        {  
            // Set the state and pull data into the cells
            thermo->setState_HP( hCells[celli]+h_ref, pCells[celli]);
            //Info << " For T=" << thermo->temperature() << "[K], P=" << pCells[celli] << nl;
            
            TCells[celli] = thermo->temperature();
            CpCells[celli] = thermo->cp_mass();
            CvCells[celli] = thermo->cv_mass();
            psiCells[celli] = thermo->density()/pCells[celli];
            muCells[celli] = solution_.get()->transport()->viscosity();
            kappaCells[celli] = solution_.get()->transport()->thermalConductivity();
        }
        else
        {    
            //Info << "Density = " << thermo->density() << " [SI?]" << nl;
            // Set the state via enthalpy and pull data into the cells
            double u = hCells[celli];
            double h = u + pCells[celli] / thermo->density() + h_ref;
            thermo->setState_HP( h, pCells[celli]);
            
            TCells[celli] = thermo->temperature();
            CpCells[celli] = thermo->cp_mass();
            CvCells[celli] = thermo->cv_mass();
            psiCells[celli] = thermo->density()/pCells[celli];
            muCells[celli] = solution_.get()->transport()->viscosity();
            kappaCells[celli] = solution_.get()->transport()->thermalConductivity();
            
            //Info << " For T=" << TCells[celli] << "[K], P=" << pCells[celli] << nl;

        }
        //Info << "Set cell to T=" << TCells[celli] << ", P=" <<  pCells[celli] << nl;
        
    }
    
    // Generate references to the boundary values associated with the 
    //  volScalarFields for the various flow parameters
    
    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& CpBf =
        this->Cp_.boundaryFieldRef();

    volScalarField::Boundary& CvBf =
        this->Cv_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& kappaBf =
        this->kappa_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];

        if (pT.fixesValue())
        {
            // Split based on if the state is defined by P, T
            
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);
                    
                // Print the composition
                //Info << "Cell " << celli << " composition: " << nl;
                const auto& mixtureObj = this->mixture();
                std::string canteraComposition = printCanteraMixture(composition, mixtureObj);
                //Info << "    " << canteraComposition << nl;
                
                // Set the reference state and the state of the face
                auto thermo = solution_->thermo();
                thermo->setState_TPY(298.15, pp[facei], canteraComposition);
                double h_ref = thermo->enthalpy_mass();
                thermo->setState_TPY(pT[facei], pp[facei], canteraComposition);        

                pCp[facei] = thermo->cp_mass();
                pCv[facei] = thermo->cv_mass();
                ppsi[facei] = thermo->density()/pp[facei];

                pmu[facei] = solution_.get()->transport()->viscosity();
                pkappa[facei] = solution_.get()->transport()->thermalConductivity();
                
                // Set the energy depending on what is being used
                if (useEnthalpy_)
                {
                  phe[facei] = thermo->enthalpy_mass()-h_ref;
                }
                else
                {
                  double u = thermo->enthalpy_mass()-h_ref-(pp[facei] / thermo->density());
                  phe[facei] = u;
                  //Info << "u_absolute = " << thermo->intEnergy_mass()-u_ref << ", phe = " << phe[facei] << ", u_ref = " << u_ref << nl;
                }
            }
        }
        else
        {
            // Now, if we have a state determined by energy, we set the state as usual
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);
                    
                // Print the composition
                //Info << "Cell " << celli << " composition: " << nl;
                const auto& mixtureObj = this->mixture();
                std::string canteraComposition = printCanteraMixture(composition, mixtureObj);
                //Info << "    " << canteraComposition << nl;
                
                // Set the reference state and the state of the face  
                auto thermo = solution_->thermo();
                thermo->setState_TPY(298.15, pp[facei], canteraComposition);
                double h_ref = thermo->enthalpy_mass();
                thermo->setState_TPY(pT[facei], pp[facei], canteraComposition);   
                    
                if (useEnthalpy_)
                {
                    thermo->setState_HP( phe[facei]+h_ref, pp[facei] );
                }
                else
                {
                    double u = thermo->intEnergy_mass();
                    double h = u + pp[facei] / thermo->density() + h_ref;
                    thermo->setState_HP( h, pp[facei] );
                }

                pT[facei] = thermo->temperature();

                pCp[facei] = thermo->cp_mass();
                pCv[facei] = thermo->cv_mass();
                ppsi[facei] = thermo->density()/pp[facei];

                pmu[facei] = solution_.get()->transport()->viscosity();
                pkappa[facei] = solution_.get()->transport()->thermalConductivity();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::PsiCanteraThermo<BaseThermo>::PsiCanteraThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    BaseThermo(mesh, phaseName)
{
    Foam::Info << "PsiCanteraThermo constructor from mesh & phaseName" << Foam::nl;
    
    // Access the dictionary via the mesh
    const dictionary& runDict = mesh.lookupObject<dictionary>("physicalProperties");
    const dictionary& canteraDict = runDict.subDict("cantera");
    const dictionary& setDict = runDict.subDict("thermoType");
    Foam::word energyType = setDict.lookupOrDefault<word>("energy", "sensibleEnthalpy");
    useEnthalpy_ = (energyType == "sensibleEnthalpy" || energyType == "enthalpy");
    Foam::Info << "Use Enthalpy? " << useEnthalpy_ << Foam::nl;
    
    // Read some properties
    word yamlFile = canteraDict.lookupOrDefault<word>("dataFile", "gri30.yaml");
    word phase = canteraDict.lookupOrDefault<word>("phaseName", "gri30");
    Foam::Info << "Cantera YAML file: " << yamlFile << ", phase: " << phase << Foam::nl;
    
    // Create the solution with default (ambient) state3
    solution_ = newSolution(yamlFile, phase);
    
    // Optionally set a default temperature & pressure
    auto thermo = solution_->thermo();
    //thermo->setState_TP(300.0, 101325.0);  // 300 K, 1 atm
    
    calculate();

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::PsiCanteraThermo<BaseThermo>::~PsiCanteraThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BaseThermo>
void Foam::PsiCanteraThermo<BaseThermo>::correct()
{
    if (BaseThermo::debug)
    {
        InfoInFunction << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();

    if (BaseThermo::debug)
    {
        Info<< "    Finished" << endl;
    }
}


// ************************************************************************* //
