Notes on code and organization
Author: GA
Last updated 7/13/2020

#######
GENERAL NOTES:
#######
-"tso" stands for "Teacher and Student Outcomes"

-The code was based off the model in "Teacher Model2019-taber.pdf," which appeared to be the most recent iteration of the model made in January 2019. Some modifications to the model
were made as I went along that I updated along the way. I will document these changes in short order. 

-The code is best viewed in a text editor with Unicode support, as the code features the frequent use of greek letters so that the code looks as close to the underlying text as possible.


#######
CODE NOTES:
#######

The following files are largely empty and fill be written later when appropriate:

------tso_master.jl: This will ultimately be the file that runs model estimation and all policy counterfactuals in one go.

------tso_simulate.jl: This will hold functions that perform data simulations after all value functions have been solved.

------tso_estimate.jl: This will hold functions that will assess the simulated data's fit of whatever auxiliary models we choose to include down the road. This will 
also define the objective function we will seek to minimize in estimation.

------tso_testing.jl: An internal file for me to play around in when testing various features/deficiencies of the code.

------tso_readin.jl: Will read in whatever external data, moments, or other useful stuff when we have it.

###

The following have had stuff written in them and are in something of a "first draft" state. They have not been debugged and so probably won't run, but they are heavily commented so that
interested readers can see how I'm going about this:

------tso_structs: Defines "structs" (Julia objects that can package a bunch of variables; similar to modules in FORTRAN) that will hold model parameters, state space vectors, and value functions

------tso_utilities: "Utilities" here means "useful things" instead of the typical economic meaning. Has some functions for discretizing distributions and forming interpolations if we ever use them.
Notably, also has functions that initialize model parameters for the first time and another that admits an arbitrary guess of model parameters and reformats them from a flat vector to something more manageable

------tso_background_functions.jl: Contains functions that govern wages, utility from occupations/major choices, teacher VA, costs of major/license acquisition, and teaching job offers.

------tso_model.jl: Contains functions that initialize model parameters, state spaces, and value functions before running backward recursion to solve for the value functions.

####
