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
    rkTableau

Description
    rkTableau is an object that allows for the calculation of the Runge-Kutta
    tableau of order P.
    
SourceFiles
    rkTableau.C

\*---------------------------------------------------------------------------*/

#ifndef rkTableau_H
#define rkTableau_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                          Class rkTableau Declaration
\*---------------------------------------------------------------------------*/

class rkTableau
{

    //
    // Private data
    //
    
    // The stages of RK integration - of order P
    int stages_;
    
    // The stage coefficients - a[i][j]
    std::vector<std::vector<double>> a_;
    
    // Stage weights - b[j]
    std::vector<double> b_;
    
    // Nodes - c[i]
    std::vector<double> c_;

    // Low-storage coefficients
    std::vector<double> alpha_;
    std::vector<double> beta_;
    
    // Expliciticity
    bool explicit_;
    
    //
    // Constructors
    //
    
    // Construct from the number of stages
    rkTableau(const int stages);
    
  public:
    
    //
    // Constructors
    //
    
    // Construct from the A & c
    rkTableau(const std::vector<std::vector<double>>& a,
              const std::vector<double>& c);
              
    // Construct from the A: Explicit ONLY
    rkTableau(const std::vector<std::vector<double>>& a);
    
    // Construct from Null to let the factory functions buld
    rkTableau():stages_(0), alpha_(), beta_(){   }
    
    //
    // Functions
    //
    
    // Access to the stages
    int stages() const { return stages_; } // simple and safe
    
    // Access the a coefficients
    const std::vector<std::vector<double>>& a() const { return a_; }
    
    // Access the b weights
    const std::vector<double>& b() const { return b_; }
    
    // Access to the c nodes
    const std::vector<double>& c() const { return c_; }
    
    // Access the low storage coefficients
    const std::vector<double>& alpha() const { return alpha_; }
    const std::vector<double>& beta() const { return beta_; }
    
    // Calculate the standard explicit RK-2 schemes in factory methods
    static rkTableau RK2_Heun();
    static rkTableau RK2_Midpoint();
    static rkTableau RK2_Ralston();
    
    // Calculate the standard explicit RK-4 scheme in factory methods
    static rkTableau RK4();
    
    // Calculate the standard explicit Carpenter low-storage schemes in factory methods
    static rkTableau CarpenterLSRK33();
    static rkTableau CarpenterLSRK54();
    


};

#endif
