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

#include <stdexcept>
#include <cmath> // for std::abs
#include rkTableau.H

// * * * * * * * * * * * * *  Constructors * * * * * * * * * * * * * * * * * //

/*
  This constructor is meant for when we receive data on A and c. 
  
  Note that the constructor needs to check to make sure the inputs are valid 
  according to
  $$
    c_i=\sum_j a_{ij}
  $$ 
*/
rkTableau::rkTableau(const std::vector<std::vector<double>>& a,
                     const std::vector<double>& c)
    : a_(a), c_(c), stages_(a.size()), explicit_(true)
{
    // Validate explicit scheme: lower-triangular
    for (size_t i = 0; i < stages_; ++i)
        for (size_t j = i; j < stages_; ++j)
            if (std::abs(a_[i][j]) > 1e-12)
                throw std::runtime_error("A must be lower-triangular for explicit RK");

    // Optional: validate c against row sums
    for (size_t i = 0; i < stages_; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < i; ++j)
            sum += a_[i][j];
        if (std::abs(sum - c_[i]) > 1e-12)
            throw std::runtime_error("Row sum of A does not match c[i]");
    }

    // b_ can be set later or via factory methods
}

/*
  This is the constructor for when we receive A only. Note that this must be
  for explicit schemes only.
  
  
*/
rkTableau::rkTableau(const std::vector<std::vector<double>>& a)
    : a_(a), stages_(a.size()), explicit_(true)
{
    // Compute c automatically
    c_.resize(stages_);
    for (size_t i = 0; i < stages_; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < i; ++j)
            sum += a_[i][j];
        c_[i] = sum;
    }

    // b_ still needs to be set for complete RK
}

// * * * * * * * * * * * * * *  Factories  * * * * * * * * * * * * * * * * * //

/*
  The advantage to these factory methods is that we can call up the classical
  Runge-Kutta scheme without having to go through the computation.
*/

rkTableau rkTableau::RK2_Heun()
{
    // Coefficients for Heun's method
    std::vector<std::vector<double>> A = {{0.0, 0.0},
                                          {1.0, 0.0}};
    std::vector<double> c = {0.0, 1.0};

    // Create an rkTableau object using the constructor
    rkTableau tableau(A, c);

    // Set the weights
    tableau.b_ = {0.5, 0.5};

    // Set the low-storage coefficients
    tableau.alpha_ = {1.0, 1.0};
    tableau.beta_ = {0.0, 1.0};

    // Return the object
    return tableau;
}

rkTableau rkTableau::RK2_Midpoint()
{
    // Coefficients for the Midpoint method
    std::vector<std::vector<double>> A = {{0.0, 0.0},
                                          {0.5, 0.0}};
    std::vector<double> c = {0.0, 0.5};

    // Create an rkTableau object using the constructor
    rkTableau tableau(A, c);

    // Set the weights
    tableau.b_ = {0, 1.0};

    // Set the low-storage coefficients
    tableau.alpha_ = {1.0, 1.0};
    tableau.beta_ = {0.0, 0.5};

    // Return the object
    return tableau;
}

rkTableau rkTableau::RK2_Ralston()
{
    // Coefficients for Ralston's method
    std::vector<std::vector<double>> A = {{0.0, 0.0},
                                          {2.0/3.0, 0.0}};
    std::vector<double> c = {0.0, 2.0/3.0};

    // Create an rkTableau object using the constructor
    rkTableau tableau(A, c);

    // Set the weights
    tableau.b_ = {0.25, 0.75};

    // Set the low-storage coefficients
    tableau.alpha_ = {0.75, 0.75};
    tableau.beta_ = {0.0, 2.0/3.0};


    // Return the object
    return tableau;
}

rkTableau rkTableau::RK4()
{
    // Coefficients for the classic 4th-order Runge-Kutta method
    std::vector<std::vector<double>> A = {
        {0.0, 0.0, 0.0, 0.0},
        {0.5, 0.0, 0.0, 0.0},
        {0.0, 0.5, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0}
    };

    std::vector<double> c = {0.0, 0.5, 0.5, 1.0};

    // Create an rkTableau object using the constructor
    rkTableau tableau(A, c);

    // Set the weights
    tableau.b_ = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};

    // Set the low-storage coefficients
    tableau.alpha_ = {0.50, 0.50, 1.0, 0.0};
    tableau.beta_ = {0.0, 0.50, 0.50, 1.0};

    // Return the object
    return tableau;
}
