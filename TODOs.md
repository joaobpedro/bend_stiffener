---
# BEND STIFFENER OPTIMIZER

  * [x] design the overall structure of the program
    * [x] the design has been made but needs improvement
    * [x] !!! need to separate the bend stiffener object from the bend stiffener get properties
  * [x] define inputs
  * [x] run a dummy example
    * [x] calculate the strain vs length given a tension and angle
    * [x] run maximum angles for a given tension up to a given strain value
  * [x] correct the EI calculation
    * [x] calculate Strain from the results I have
    * [x] implement non-linear material properties
  * [x] calculate the strain based on the results
  * [ ] add strain due to tension - to be discussed!!
  * [ ] !!! the convergence is broken right now, it seems that its not converging, and I do not understand how!!!
  * [x] increase the theta angle little by little to make sure we capture the non-linear behavior
  * [ ] implement optimization for a given set of load cases, minimize the strain distance.
  * [ ] implement input from the user on the bend stiffener dimensions
  * [ ] make possible for the user to should how to extract the tension/angles pairs, i.e., orcaflex post-processing
  or direct input from the user

  * [x] start with a E-mod - get strain - update E-mod - get_new_strain - (store the state in a state variable)
  * [x] apply the load in steps until equilibrium


  ## NOTES

  * [ ] the bending moment does not need to be calculated. Only the shear force.
        * [ ] knowing the shea force and the bending moment at the tip and can reverte back to the root bending moment
 
  ## PROBLEM
  * [ ] the current problem is that I am calculating the initial V0, M0 without non-linear EI
  * [ ] !!! PRIORITY is to include non-linear EI in the initial dtermination of the conditions
  * [ ] so i need an engine, which will import everything from the bend_stiffner , bs_physics , and the integrator
  so I can put together the puzzle
  * [ ] I need to pass inertia not as a values but as a function. Is that possible?
