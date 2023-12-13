#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <map>

using namespace std;

template <typename State>
class AutoMarkov {
public:
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

struct SystemState {
	bool m[3]; // healthy active module count
	bool c[3]; // healthy comparator count
	bool sv; // switch or voter
	int s; // spare count
	int TAG;

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

	SystemState initialState;
	for (int i = 0; i < 3; i++) {
		initialState.m[i] = 1;
		initialState.c[i] = 1;
	}
	initialState.sv = 1;
	initialState.s = 2;
	initialState.TAG = 1;

	AutoMarkov<SystemState> autoMarkov(initialState, generateNextState);

	autoMarkov.getMatlabEquations(1, "AutoMarkovSecondAssumption.m");

	autoMarkov.getGraph("BigGraph.txt");

	double dt = 0.001; // years
	for (double t = 0; t < 5; t += dt) {
		autoMarkov.updateprobabilities(dt);
	}

	cout << autoMarkov.getProbabilitySumWithTAG(1) << endl;

	return 0;
}