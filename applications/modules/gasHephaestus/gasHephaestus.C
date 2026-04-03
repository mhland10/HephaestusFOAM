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

#include "gasHephaestus.H"
#include "localEulerDdtScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(multicomponentFluid, 0);
    addToRunTimeSelectionTable(solver, multicomponentFluid, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::multicomponentFluid::multicomponentFluid(fvMesh& mesh)
:
    isothermalFluid
    (
        mesh,
        autoPtr<fluidThermo>(fluidMulticomponentThermo::New(mesh).ptr())
    ),

    thermo_(refCast<fluidMulticomponentThermo>(isothermalFluid::thermo_)),

    Y_(thermo_.Y()),

    reaction(combustionModel::New(thermo_, momentumTransport())),

    thermophysicalTransport
    (
        fluidMulticomponentThermophysicalTransportModel::New
        (
            momentumTransport(),
            thermo_
        )
    ),

    thermo(thermo_),
    Y(Y_)
{
    thermo.validate(type(), "h", "e");

    forAll(Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(thermo.he());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::multicomponentFluid::~multicomponentFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::multicomponentFluid::prePredictor()
{
    isothermalFluid::prePredictor();

    if (pimple.predictTransport())
    {
        thermophysicalTransport->predict();
    }
}


void Foam::solvers::multicomponentFluid::postCorrector()
{
    isothermalFluid::postCorrector();

    if (pimple.correctTransport())
    {
        thermophysicalTransport->correct();
    }
}

/*
  This function pulls the data from the object (and ones it inherits from) and 
makes hard copies to store in the W solver state to allow more advanced time
stepping functionality, namely Carpenter-style Runge-Kutta schemes.

  Returns:
      W (SolverState):  This is a data structure that contains the data that 
                        defines the state for this specific CFD solver. 
                        Contains:
                        
                          U (vector Field):    Velocity vector
                          p (scalar Field):    Pressure field
                          he (scalar Field):   Energy field, depends on 
                                               selected thermo model
                          Yi[i] (scalar List): The list of mass fractions for
                                               the multicomponent composition 

*/
Foam::solvers::multicomponentFluid::SolverState
Foam::solvers::multicomponentFluid::outputState() const
{
    SolverState W;
    
    //
    // Copy vector/fields
    //
    
    // These two come from isothermalFluid
    W.U = Field<vector>(U_.internalField()); 
    W.p = Field<scalar>(p_.internalField());
    
    // This also comes from isothermalFluid via the thermo model
    W.he = Field<scalar>(thermo_.he().internalField());
    
    // Copy composition fields
    W.Yi.setSize(Y_.size());
    for (label i = 0; i < Y_.size(); ++i)
    {
        // Note that THIS solver contains the Y's
        W.Yi[i] = Field<scalar>(Y_[i].internalField());
    }
    
    
    
    return W;
    
}

/*
  This function takes one state (the input) and compares it to the current 
state of the solver. The difference is returned as dW. This allows a time 
stepping where the handling of the state is more flexible, and the time stepper
can merely get the state and the difference for its operations, like in a 
Carpenter-style RK scheme.

  Inputs:
    W0 (SolverState reference):  This is the state to compare against. See 
                                 header for details as to what this contains.
                               
  Returns:
    dW (SolverState):  This is the difference in the current state and the 
                       original state we are comparing against.
                       
*/
Foam::solvers::multicomponentFluid::SolverState 
Foam::solvers::multicomponentFluid::deltaState(const SolverState& W0) const
{
    // Get the current solver state
    SolverState W = this->outputState();

    SolverState dW;

    // Subtract each component to form the delta
    dW.U  = W.U - W0.U;
    dW.p  = W.p - W0.p;
    dW.he = W.he - W0.he;

    // Allocate composition list
    dW.Yi.setSize(W.Yi.size());

    // Subtract each species mass fraction
    for (label i = 0; i < W.Yi.size(); ++i)
    {
        dW.Yi[i] = W.Yi[i] - W0.Yi[i];
    }

    return dW;
}

/*
  This state takes the difference between the reference state and new state and
blends the difference with the reference state according to alpha (mulitplier 
for the reference state) and beta (multiplier for the state difference). This 
allows more flexible time stepping schemes, particularly the Carpenter-style
Runge-Kutta schemes, where we are interested in letting the PIMPLE handle the
hyperbolic-elliptic coupling to close the state, while allowing a high order
time stepping scheme. This follows the formula:

$$
W^{(s)}=\alpha_sW^{(0)}+\beta_s\left( W^{(s-1)}+\Delta tR(W^{(s-1)} \right)
$$

  Inputs:
    W0 (SolverState reference):  This is the reference state for the time step.
    
    dW (SolverState reference):  This is the difference in states for the time
                                 step.
    
    alpha (scalar):  This is the multiplier for the reference state.
    
    beta (scalar):   This is the multiplier for the difference in states.
    
  Returns:
    bW (SolverState):  This is the blended state that is the result of the time
                       stepping scheme.
                       
*/
Foam::solvers::multicomponentFluid::SolverState 
Foam::solvers::multicomponentFluid::blendState( 
  const SolverState& W0,
  const SolverState& dW, 
  scalar alpha, 
  scalar beta             )
  const
{
    // Define the solver state to blend into
    SolverState bW;
    
    // 
    //  Set the state
    //
    
    // Set velocity
    bW.U = alpha * W0.U + beta * dW.U;
    
    // Set pressure
    bW.p = alpha * W0.p + beta * dW.p;
    
    // Set energy
    bW.he = alpha * W0.he + beta * dW.he;
    
    // Set composition fields
    bW.Yi.setSize(W0.Yi.size());
    for (label i = 0; i < W0.Yi.size(); ++i)
    {
        bW.Yi[i] = alpha * W0.Yi[i] + beta * dW.Yi[i];
    }



    return bW;
        
}

/*
  This function re-sets the state based on the time stepping operations. This 
includes reconciling the thermophysical state and running the post-corrector
stage.

  Inputs:
    W (SolverState reference):  The solvers state to set the state from.
    
  Returns:
    **void**
    
*/
void Foam::solvers::multicomponentFluid::restoreState(const SolverState& W)
{
    
    //
    //  Set the OpenFOAM field values
    //
    
    U_.primitiveFieldRef() = W.U;
    p_.primitiveFieldRef() = W.p;
    thermo_.he().primitiveFieldRef() = W.he;
    
    for (label i = 0; i < Y_.size(); ++i)
    {
        Y_[i].primitiveFieldRef() = W.Yi[i];
    }
    
    // Correct the Bcs
    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    thermo_.he().correctBoundaryConditions();
    
    // Ensure thermo consistency
    thermo_.correct();

    // Re-run any necessary post corrections
    this->postCorrector();
    
}




// ************************************************************************* //
