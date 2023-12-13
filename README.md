# AutoMarkov
AutoMarkov is a C++ program designed to analyze and solve Markov models represented as state machines. This tool is capable of performing numerical solutions of probabilistic models and can generate the corresponding differential equations to be used with MATLAB for further analysis.

## Overview
Markov models are powerful tools for modeling the behavior and performance of systems subject to random state transitions, such as reliability and availability models in engineering. AutoMarkov automates the process of solving these models by numerically simulating the state transitions over time and providing differential equations that represent the model's dynamics.

## General Structure of the Code

### `AutoMarkov` Class

The core of the program is encapsulated in the `AutoMarkov` class. The class requires an initial state and a function that defines state transitions (`generateNextState`) upon instantiation.

### State Representation

States are defined as `struct`s containing at least one variable named `TAG`. This `TAG` categorizes states, typically distinguishing between failed states (`TAG = 0`) and operational states (`TAG = 1`). This categorization is used in subsequent calculations and the generation of MATLAB equations.

### Key Methods

- `getMatlabEquations(int TAG, const string& directory)`: Uses the specified `TAG` to filter states and generates MATLAB-friendly equations. For example, to calculate reliability, sum the probabilities of operational states (typically with `TAG = 1`) and generate corresponding MATLAB code.
- `getGraph(const string& directory)`: Outputs a textual description of the Markov model graph, which can be visualized using online graph editor tools.
- `updateProbabilities(double dt)`: Updates the probabilities of each state based on the defined transitions and a small timestep `dt`. This method is called iteratively to simulate the model over a specified duration.
- `getProbabilitySumWithTAG(int TAG)`: After simulating the model, this function calculates the sum of the probabilities of states with the specified `TAG`.

## Usage

1. **Define the State Structure**: Customize the `SystemState` struct to represent the states in your Markov model, including any necessary variables and a `TAG`.
2. **Implement State Transitions**: Write the `generateNextState` function to describe the possible transitions from any given state and their respective probabilities.
3. **Initialize the Model**: Create an `AutoMarkov` object with an initial state and the `generateNextState` function.
4. **Run the Simulation**: Use the `updateProbabilities` method to simulate the model over time.
5. **Output Results**: Utilize the `getProbabilitySumWithTAG` method to obtain reliability metrics, and the `getMatlabEquations` to generate MATLAB code for further analysis.

### Example

```cpp
// Define your initial state and transitions
SystemState initialState = {/* ... initialization ... */};
AutoMarkov<SystemState> autoMarkov(initialState, generateNextState);

// Run the model
autoMarkov.updateProbabilities(0.001); // Update probabilities with a timestep of 0.001

// Analyze the results
double reliability = autoMarkov.getProbabilitySumWithTAG(1);
cout << "Reliability after simulation: " << reliability << endl;

// Generate MATLAB equations
autoMarkov.getMatlabEquations(1, "MarkovModelEquations.m");
```

## Visualization

To visualize the state transition graph, run the `getGraph` method and import the resulting text file into a graph editor tool like the one available at https://csacademy.com/app/graph_editor/.
