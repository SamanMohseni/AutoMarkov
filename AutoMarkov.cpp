#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <map>

using namespace std;

// Generic AutoMarkov class capable of handling any state type.
template <typename State>
class AutoMarkov {
public:
	// Constructor takes an initial state and a function pointer for state transitions.
	AutoMarkov(State init_state, vector<pair<State, double>>(*next_states)(State)) {
		int index = 0;
		stateIndex[init_state] = index;
		stateByIndex[index] = init_state;
		index++;

		queue<State> states;
		states.push(init_state);

		while (states.size()) {
			vector<pair<State, double>> res = next_states(states.front());
			adjacencyList.resize(adjacencyList.size() + 1);

			for (int i = 0; i < res.size(); i++) {
				if (!stateIndex.count(res[i].first)) {
					states.push(res[i].first);
					stateIndex[res[i].first] = index;
					stateByIndex[index] = res[i].first;
					index++;
				}
				adjacencyList[stateIndex[states.front()]].push_back(make_pair(stateIndex[res[i].first], res[i].second));
			}

			states.pop();
		}

		probability.resize(index);
		fill(probability.begin(), probability.end(), 0);
		probability[0] = 1;
	}

	// Generates MATLAB equations for the Markov model, given a specific TAG.
	void getMatlabEquations(int TAG, const char* dir) {
		ofstream file_out;
		file_out.open(dir);

		file_out << "syms";
		for (int i = 0; i < adjacencyList.size(); i++) {
			file_out << " p" << i << "(t)";
		}
		file_out << ";" << endl << endl;

		vector<vector<pair<int, double>>> incomeList;
		incomeList.resize(adjacencyList.size());
		for (int i = 0; i < adjacencyList.size(); i++) {
			for (int j = 0; j < adjacencyList[i].size(); j++) {
				int neighbourIndex = adjacencyList[i][j].first;
				double lambda = adjacencyList[i][j].second;
				incomeList[neighbourIndex].push_back(make_pair(i, lambda));
			}
		}

		for (int i = 0; i < adjacencyList.size(); i++) {
			file_out << "eqn" << i << " = diff(p" << i << ") == ";

			for (int j = 0; j < incomeList[i].size(); j++) {
				int incomeIndex = incomeList[i][j].first;
				double lambda = incomeList[i][j].second;
				if (j) {
					file_out << " + ";
				}
				file_out << lambda << "*p" << incomeIndex;
			}

			double lambdaSum = 0;
			for (int j = 0; j < adjacencyList[i].size(); j++) {
				double lambda = adjacencyList[i][j].second;
				lambdaSum += lambda;
			}

			file_out << " - " << lambdaSum << "*p" << i << ";" << endl;
		}

		file_out << endl << "eqns = [";
		for (int i = 0; i < adjacencyList.size(); i++) {
			if (i) {
				file_out << "; ";
			}
			file_out << "eqn" << i;
		}
		file_out << "];" << endl << endl;

		file_out << "conds = [";
		file_out << "p0(0) == 1";
		for (int i = 1; i < adjacencyList.size(); i++) {
			file_out << ", p" << i << "(0) == 0";
		}
		file_out << "];" << endl << endl;

		file_out << "S = dsolve(eqns, conds);" << endl << endl;

		file_out << "R(t) = ";
		bool first = true;
		for (int i = 0; i < probability.size(); i++) {
			if (stateByIndex[i].TAG == TAG) {
				if (!first) {
					file_out << " + ";
				}
				else {
					first = false;
				}
				file_out << "S.p" << i;
			}
		}

		file_out << endl << endl;

		file_out.close();
	}

	// Updates the probabilities of each state based on the transition rates and time step.
	// Used for numerical computation of probabilities.
	void updateprobabilities(double dt) {
		vector<double> newProbability(probability.size(), 0);

		for (int i = 0; i < adjacencyList.size(); i++) {
			double lambdaSum = 0;
			for (int j = 0; j < adjacencyList[i].size(); j++) {
				int neighbourIndex = adjacencyList[i][j].first;
				double lambda = adjacencyList[i][j].second;
				newProbability[neighbourIndex] += lambda * dt * probability[i];
				lambdaSum += lambda;
			}
			newProbability[i] += (1 - lambdaSum * dt) * probability[i];
		}

		probability = newProbability;
	}
	
	// Outputs the graph of the Markov model in a textual format.
	void getGraph(const char* dir) {
		ofstream file_out;
		file_out.open(dir);
		for (int i = 0; i < adjacencyList.size(); i++) {
			for (int j = 0; j < adjacencyList[i].size(); j++) {
				int neighbourIndex = adjacencyList[i][j].first;
				double lambda = adjacencyList[i][j].second;
				file_out << i << "_" << stateByIndex[i].TAG << " ";
				file_out << neighbourIndex << "_" << stateByIndex[neighbourIndex].TAG << " ";
				file_out << lambda << "dt" << endl;
			}
		}
		file_out.close();
	}

	// Calculates the sum of probabilities for states with a given TAG.
	double getProbabilitySumWithTAG(int TAG) {
		double sum = 0;
		for (int i = 0; i < probability.size(); i++) {
			if (stateByIndex[i].TAG == TAG) {
				sum += probability[i];
			}
		}
		return sum;
	}

private:
	vector<vector<pair<int, double>>> adjacencyList; // {<index of neighbour, transition lambda>} [index of state]

	map<State, int> stateIndex;
	map<int, State> stateByIndex;

	vector<double> probability;
};

// The SystemState struct must represent a state in your Markov model.
// Modify the members of this struct to reflect the components and characteristics
// of the states in your model. For example, you might include system configuration,
// operational status, or any other relevant state descriptors.
struct SystemState {
	bool m[3]; // healthy active module count
	bool c[3]; // healthy comparator count
	bool sv; // switch or voter
	int s; // spare count
	int TAG; // Used to categorize states (e.g., operational vs. failed)

	// Serialize method must uniquely represent the state.
    	// If you modify the state variables, ensure this method generates
    	// a unique integer for each possible state configuration.
	int serialize() const {
		int serial = 0;
		for (int i = 0; i < 3; i++) {
			serial |= m[i] << i;
			serial |= c[i] << (i + 3);
		}
		serial |= sv << 6;
		serial |= s << 7;

		return serial;
	}

	bool operator==(const SystemState& other) const {
		return serialize() == other.serialize();
	}

	bool operator<(const SystemState& other)  const {
		return serialize() < other.serialize();
	}
};

// The generateNextState function defines the transition logic for your Markov model.
// Modify this function to reflect the possible transitions between states in your model,
// including the conditions that cause transitions and the rates at which they occur.
vector<pair<SystemState, double>> generateNextState(const SystemState state) {
	const double module_lambda = 0.1;
	const double others_lambda = 0.025;

	vector<pair<SystemState, double>> nextStates;
	SystemState next;

	if (state.TAG == 0) {
		return nextStates;
	}

	next = state;
	next.sv = 0;
	next.TAG = 0;
	nextStates.push_back(make_pair(next, 2 * others_lambda));

	for (int i = 0; i < 3; i++) {
		if (state.c[i]) {
			next = state;
			next.c[i] = 0;
			nextStates.push_back(make_pair(next, others_lambda));
		}
		if (state.m[i]) {
			next = state;
			if (state.s) {
				if (state.c[i]) {
					next.s--;
				}
				else {
					next.m[i] = 0;
				}
			}
			else {
				next.m[i] = 0;
			}
			next.TAG = (next.m[0] + next.m[1] + next.m[2]) > 1;
			nextStates.push_back(make_pair(next, module_lambda));
		}
	}

	return nextStates;
}

int main() {
	// Define the initial state of your model here by setting the appropriate state variables.
	SystemState initialState;
	for (int i = 0; i < 3; i++) {
		initialState.m[i] = 1;
		initialState.c[i] = 1;
	}
	initialState.sv = 1;
	initialState.s = 2;
	initialState.TAG = 1;

	// Create an instance of the AutoMarkov class with your initial state and transition function.
	AutoMarkov<SystemState> autoMarkov(initialState, generateNextState);

	// Modify the parameters of getMatlabEquations method call to match your desired output.
    	// TAG should correspond to a category of states you are interested in.
	autoMarkov.getMatlabEquations(1, "AutoMarkovSecondAssumption.m");

	autoMarkov.getGraph("BigGraph.txt");

	// Define the timestep for probability updates and the duration over which to update.
    	double dt = 0.001; // years
	// This method numerically solves the Markov model over the defined timestep 'dt'.
    	// Run a simulation loop to update the probabilities of each state at each timestep 'dt'.
   	// This loop advances the model's state probabilities over a simulated time of 5 years.
    	// Modify 'dt' and the loop range to match the desired timestep and simulation duration.
	for (double t = 0; t < 5; t += dt) {
		autoMarkov.updateprobabilities(dt);
	}
	
	// Calculate and output the sum of the probabilities of all states with a given TAG.
    	// This final probability sum provides a measure of the likelihood of being in any
    	// state categorized by the specified TAG (e.g., TAG = 1 for operational states) after
    	// the simulation has concluded. This can be used to assess the reliability or other
    	// performance metrics of the model over the simulated time period.
	cout << autoMarkov.getProbabilitySumWithTAG(1) << endl;

	return 0;
}
