# AutoMarkov
A C++ program that, given the description of a Markov model as a state machine, can numerically solve the model and also generate the differential equations in MATLAB.

# General Structure of the Code

The general structure of this code consists of a class named `AutoMarkov` which takes the initial state and a pointer to a function in its constructor. This function (`generateNextState`) takes the current state in the Markov model as input and generates the subsequent states along with their probabilities of occurrence. Using this function, the traversal from the start state begins, and the graph related to the Markov model is constructed.

The states in the Markov model are `struct`s that must contain a variable named `TAG`. When states are generated in the `generateNextState` function, this `TAG` should also be initialized; for instance, assign a value of 0 to states we consider as fail, and a value of 1 to other states.

Additionally, there is a function in this class named `getMatlabEquations` which takes the `TAG` value and the directory for saving the MATLAB file as inputs. This `TAG` is used, for example, when we want to calculate reliability, as the probability equations related to the healthy states need to be summed. If we have previously assigned a `TAG` value of 1 to healthy states, the function `getMatlabEquations` will sum the healthy states with the `TAG` value of 1 so that MATLAB can also compute the result of this sum of equations.

Moreover, another function in this class, `getGraph`, provides a textual description of the Markov model graph that can be rendered in tools such as https://csacademy.com/app/graph_editor/.

This class also facilitates the numerical calculation of state probabilities using the `updateProbabilities` function. This is done by calculating the next values based on the current probabilities with small time steps and, after advancing to the desired time, the sum of numerical states with a given `TAG` can be computed with the `getProbabilitySumWithTAG` function. In this way, reliability can be calculated for specific times with high speed and accuracy.
