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

#include <stdexcept>
#include <cmath> // for std::abs
#include "rkTableau.H"

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

// Constructor from number of stages
rkTableau::rkTableau(const int stages)
    : stages_(stages), a_(stages), b_(stages), c_(stages),
      alpha_(stages), beta_(stages), explicit_(true)
{
    // Initialize the vectors
    for (int i = 0; i < stages_; ++i)
        a_[i].resize(stages_, 0.0);
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

rkTableau rkTableau::CarpenterLSRK33()
{
    rkTableau rk(3);

    rk.alpha_ = {
                0.0,
               -5.0/9.0,
               -153.0/128.0
            };

    rk.beta_ = {
                1.0/3.0,
                15.0/16.0,
                8.0/15.0
            };
    
    rk.c_ = {
                0.0,                // stage 1
                1.0/3.0,            // stage 2
                5.0/9.0             // stage 3
            };
    rk.c_[0] = 1.0 - rk.c_[1] - rk.c_[2];

    return rk;
}

rkTableau rkTableau::CarpenterLSRK54()
{
    rkTableau rk(5);

    rk.alpha_ = {
        0.0,
       -567301805773.0/1357537059087.0,
       -2404267990393.0/2016746695238.0,
       -3550918686646.0/2091501179385.0,
       -1275806237668.0/842570457699.0
    };

    rk.beta_ = {
        1432997174477.0/9575080441755.0,
        5161836677717.0/13612068292357.0,
        1720146321549.0/2090206949498.0,
        3134564353537.0/4481467310338.0,
        2277821191437.0/14882151754819.0
    };
    
    // Hard-coded stage times (c coefficients)
    rk.c_ = {
        0.0,                                     // stage 1
        1432997174477.0/9575080441755.0,        // stage 2
        2526269341429.0/6820363962896.0,        // stage 3
        2006345519317.0/3224310063776.0,        // stage 4
        2802321613138.0/2924317926251.0         // stage 5
    };
    rk.c_[0] = 1.0 - rk.c_[1] - rk.c_[2] - rk.c_[3] - rk.c_[4];

    return rk;
}
