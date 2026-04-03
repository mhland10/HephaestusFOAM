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
    multistepWrapper

Description
    To make a solver compatible with multistep time stepping, namely Carpenter-
    style RK time stepping, we need to wrap the solver so it is compatible with
    the multistep time stepping loop.
    
    The class templates the "SolverType" which allows us to put just about any
    OpenFOAM-compatible solver into the wrapper.

See also
    hephaestusRKFoam.C

\*---------------------------------------------------------------------------*/

#ifndef multistepWrapper_H
#define multistepWrapper_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

// Forward define a default SolverState structure for non-compatible solvers
struct DefaultSolverState
{
    Field<vector> U;
    Field<scalar> p;
};

/*---------------------------------------------------------------------------*\
  Trait to detect if a SolverType has a SolverState member
\*---------------------------------------------------------------------------*/

// C++14-compatible void_t
template<typename... Ts> struct make_void { using type = void; };
template<typename... Ts> using void_t = typename make_void<Ts...>::type;

// Trait to detect SolverState
template<typename, typename = void>
struct hasSolverState : std::false_type {};

template<typename T>
struct hasSolverState<T, void_t<typename T::SolverState>> : std::true_type {};

// Selector for SolverState
template<typename SolverType, bool hasState = hasSolverState<SolverType>::value>
struct SolverStateSelector;

// Specialization if SolverType has SolverState
template<typename SolverType>
struct SolverStateSelector<SolverType, true>
{
    using type = typename SolverType::SolverState;
};

// Specialization if SolverType does NOT have SolverState
template<typename SolverType>
struct SolverStateSelector<SolverType, false>
{
    using type = DefaultSolverState;
};

// Detect outputState()
template<typename T>
class hasOutputState
{
    typedef char yes[1];
    typedef char no[2];

    template<typename U, U> struct reallyHas;

    template<typename C>
    static yes& test(reallyHas<SolverState (C::*)() const, &C::outputState>*);

    template<typename>
    static no& test(...);

public:
    static constexpr bool value = sizeof(test<T>(0)) == sizeof(yes);
};

/*---------------------------------------------------------------------------*\
                      Class multistepWrapper Declaration
\*---------------------------------------------------------------------------*/

template<typename SolverType>
class multistepWrapper
{

public:

    //- Make the solver publicly accessible for external use
    SolverType& solver_;
    
    //- Here, we allow the object to switch between the default and inherited
    //  SolverState depending on if the incoming "solver" object was built for
    //  multitep time stepping
    using SolverState = typename SolverStateSelector<SolverType>::type;
    
    multistepWrapper(SolverType& s) : solver_(s) {}
    
    SolverState outputState() const
    {
        if constexpr (hasOutputState<SolverType>::value)
        {
            // Call the solver's own function
            return solver_.outputState();
        }
        else
        {
            // Use the wrapper's state_
            return state_;
        }
    }
    
    // outputState() wrapper
    SolverState outputState() const
    {
        return outputStateImpl(std::integral_constant<bool, hasOutputState<SolverType>::value>());
    }

private:
    // Case 1: solver has outputState()
    SolverState outputStateImpl(std::true_type) const
    {
        return solver_.outputState();
    }

    // Case 2: solver does NOT have outputState()
    SolverState outputStateImpl(std::false_type) const
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
    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //