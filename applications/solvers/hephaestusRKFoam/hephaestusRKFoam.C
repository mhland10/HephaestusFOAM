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

Application
    foamRun

Description
    Loads and executes an OpenFOAM solver module either specified by the
    optional \c solver entry in the \c controlDict or as a command-line
    argument.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

Usage
    \b hephaestusRKFoam [OPTION]

      - \par -solver <name>
        Solver name

      - \par -libs '(\"lib1.so\" ... \"libN.so\")'
        Specify the additional libraries loaded

    Example usage:
      - To run a \c rhoPimpleFoam case by specifying the solver on the
        command line:
        \verbatim
            foamRun -solver fluid
        \endverbatim

      - To update and run a \c rhoPimpleFoam case add the following entries to
        the controlDict:
        \verbatim
            application     foamRun;

            solver          fluid;
        \endverbatim
        then execute \c foamRun

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "solver.H"
#include "pimpleSingleRegionControl.H"
#include "setDeltaT.H"
#include "rkTableau.H"
#include "multistepWrapper.H"
#include <vector>
#include <type_traits> // for std::remove_reference
#include "string.H"
#include "stringOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "solver",
        "name",
        "Solver name"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    
    // pre-set the control dictionary pointer
    const dictionary& controlDict = runTime.controlDict();

    // Read the solverName from the optional solver entry in controlDict
    word solverName
    (
        runTime.controlDict().lookupOrDefault("solver", word::null)
    );

    // Optionally reset the solver name from the -solver command-line argument
    args.optionReadIfPresent("solver", solverName);

    // Check the solverName has been set
    if (solverName == word::null)
    {
        args.printUsage();

        FatalErrorIn(args.executable())
            << "solver not specified in the controlDict or on the command-line"
            << exit(FatalError);
    }
    else
    {
        // Load the solver library
        solver::load(solverName);
    }

    // Create the default single region mesh
    #include "createMesh.H"

    // Instantiate the selected solver
    autoPtr<solver> solverPtr(solver::New(solverName, mesh));
    solver& solver = solverPtr();  // Keep the name "solver"
    
    // Wrap using the type of the variable
    solvers::multistepWrapper<std::remove_reference<decltype(solver)>::type> wrappedSolver(solver);
    
    //
    //  Instantiate the Runge-Kutta scheme
    //
    word rkScheme;
    if (controlDict.found("rkScheme"))
    {
        // Take the specified RK scheme in the control dictionary under entry "rkScheme"
        controlDict.lookup("rkScheme") >> rkScheme;
    
        Info<< "Runge-Kutta scheme specified in controlDict: "
            << rkScheme << nl << endl;
    }
    else
    {
        // If there is no specification, default to the Carpenter Low-Storage 3-Stage scheme
        rkScheme = "CarpenterLS3";
    
        Info<< "No 'rkScheme' specified in controlDict." << nl
            << "Defaulting to CarpenterLS3 (3-stage low-storage RK)." 
            << nl << endl;
    }

    // Create the outer PIMPLE loop and control structure
    pimpleSingleRegionControl pimple(solver.pimple);

    // Set the initial time-step
    setDeltaT(runTime, solver);
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //                          Set up RK Coefficients                       //
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    word s = rkScheme;

    // Convert to lower case and remove spaces
    // Convert RK scheme to lowercase and remove spaces manually
    for (auto& c : rkScheme)
        c = std::tolower(c);
    
    rkScheme.erase(
        std::remove(rkScheme.begin(), rkScheme.end(), ' '),
        rkScheme.end()
    );
    
    rkTableau rk_;   // uses the null constructor

    if (s == "rk2" || s == "midpoint" || s == "rk2_midpoint")
    {
        Info << nl << "Using RK2 Midpoint scheme\n" << endl;
        rk_ = rkTableau::RK2_Midpoint();
    }
    else if (s == "heun" || s == "rk2_heun")
    {
        Info << nl << "Using RK2 Heun scheme\n" << endl;
        rk_ = rkTableau::RK2_Heun();
    }
    else if (s == "ralston" || s == "rk2_ralston")
    {
        Info << nl << "Using RK2 Ralston scheme\n" << endl;
        rk_ = rkTableau::RK2_Ralston();
    }
    else if (s == "rk4")
    {
        Info << nl << "Using RK4 Classic scheme\n" << endl;
        rk_ = rkTableau::RK4();
    }
    else if (s == "carpenterls3" || s == "lsrk33")
    {
        Info << nl << "Using RK3 Carpenter low-storage scheme\n" << endl;
        rk_ = rkTableau::CarpenterLSRK33();
    }
    else if (s == "carpenterls5" || s == "lsrk54")
    {
        Info << nl << "Using RK5 Carpenter low-storage scheme\n" << endl;
        rk_ = rkTableau::CarpenterLSRK54();
    }
    else
    {
        Info << nl << "Invalid scheme specification, defaulting back to RK3 Carpenter low-storage scheme\n" << endl;
        rk_ = rkTableau::CarpenterLSRK33();        
        
    }
    
    //
    //  Pull the coefficients
    //
    const auto& alpha = rk_.alpha();
    const auto& beta  = rk_.beta();
    const auto& c  = rk_.c();
    
    const label nStages = rk_.stages();
    
    Info << "Using RK scheme: " << rkScheme << nl
         << "Number of stages: " << nStages << nl << endl;
         
    // Generic fallback state structure
    struct SolverState
    {
      Field<vector> U;
      Field<scalar> p;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << nl << "Starting time loop\n" << endl;

    while (pimple.run(runTime))
    {  
        
        
        // --------------------------------------
        // Loop through RK stage data
        // --------------------------------------
        Info << "Low-storage RK coefficients for this time-step:\n";
    
        const auto& alpha = rk_.alpha();
        const auto& beta  = rk_.beta();
        const auto& c = rk_.c();
        const label nStages = rk_.stages();
    
        for (label stage=0; stage<nStages; ++stage)
        {
            Info << "  Stage " << stage+1 << ": "
                 << "alpha = " << alpha[stage] << ", "
                 << "beta  = " << beta[stage] << ", "
                 << "c  = " << c[stage] << nl;
        }
        
        // Update PIMPLE outer-loop parameters if changed
        pimple.read();

        solver.preSolve();

        // Adjust the time-step according to the solver maxDeltaT
        adjustDeltaT(runTime, solver);
        
        // --------------------------------------
        // Loop through RK stages
        // --------------------------------------
        scalar baseDeltaT = runTime.deltaTValue();
        for (label stage=0; stage<nStages; ++stage)
        {

            scalar stageDeltaT = c[stage] * baseDeltaT;
            
            Info << "**Stage " << stage+1 << ": "
                 << "dt = " << stageDeltaT << nl;
                 
            // Set the stage step according to the standard OpenFOAM 
            // methodology
            runTime.setDeltaT(stageDeltaT);
            runTime++;
            Info<< "Time = " << runTime.userTimeName() << nl << endl;
            
            // Store the solver state to separate the PIMPLE operations
            SolverState W0; 
            
        }
        
        runTime.setDeltaT(baseDeltaT);
        
        

        Info<< "Time = " << runTime.userTimeName() << nl << endl;
        
        
        
        // PIMPLE corrector loop
        while (pimple.loop())
        {
            solver.moveMesh();
            solver.motionCorrector();
            solver.fvModels().correct();
            solver.prePredictor();
            solver.momentumPredictor();
            solver.thermophysicalPredictor();
            solver.pressureCorrector();
            solver.postCorrector();
        }
        
        

        solver.postSolve();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
    
}


// ************************************************************************* //
