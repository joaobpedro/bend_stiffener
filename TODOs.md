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
