#ifndef STATE_H
#define STATE_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>


// a state machine is a collection of states
// each state has a list of transitions
class State
{
public:
	State(std::string name) :name(name) {};



	~State();

	std::unordered_map<std::string, State*> transitions;
	std::string name;


};

class StateMachine
{
public:
	StateMachine() {};
	~StateMachine() {};

	void add_state(State* state) {
		states.push_back(state);
	};



	void add_transition(State* from, State* to, std::string name) {
		from->transitions[name] = to;
	};
	void add_transition(std::string from, std::string to, std::string name) {
		for (int i = 0; i < states.size(); i++) {
			if (states[i]->name == from) {
				for (int j = 0; j < states.size(); j++) {
					if (states[j]->name == to) {
						states[i]->transitions[name] = states[j];
					}
				}
			}
		}
	};
	void add_transition(std::string from, State* to, std::string name) {
		for (int i = 0; i < states.size(); i++) {
			if (states[i]->name == from) {
				states[i]->transitions[name] = to;
			}
		}
	};
	void add_transition(State* from, std::string to, std::string name) {
		for (int i = 0; i < states.size(); i++) {
			if (states[i]->name == to) {
				from->transitions[name] = states[i];
			}
		}
	};

	void set_current_state(State* state) {
		currentState = state;
	}
	void set_current_state(std::string name) {
		for (int i = 0; i < states.size(); i++) {
			if (states[i]->name == name) {
				currentState = states[i];
				return;
			}
		}
		std::cout << "state not found" << std::endl;
	}

	std::string input; // the input to the state machine, this is the name of the transition, can only be one transition at a time
	std::vector<State*> states; // all the states in the state machine
	State* currentState; // the current state of the state machine
	void set_input(std::string input) {
		this->input = input;
	}


	State* get_current_state() {
		return currentState;
	}

	void update() {// according to the input, update new state
		if(currentState->transitions.count(input) != 0)// if the input is a valid transition
			currentState = currentState->transitions[input];
	}

	void print_state() {
		std::cout<< "current state:" << currentState->name << std::endl;
	}
};


class ExampleStateMachine : public StateMachine
{
public:
	ExampleStateMachine() {
		State* state1 = new State("state 1");
		State* state2 = new State("state 2");
		State* state3 = new State("state 3");



		add_state(state1);
		add_state(state2);
		add_state(state3);
		


		add_transition(state1, state2, "to_2");
		add_transition(state1, state3, "to_3");

		add_transition(state3, state1, "to_1");
		add_transition(state3, state2, "to_2");

		add_transition(state2, state1, "to_1");
		add_transition(state2, state3, "to_3");
		


		set_current_state(state1);
	}
	~ExampleStateMachine() {};

	void request_state2() {
		set_input("to_2");
	}
	void request_state1() {
		set_input("to_1");
	}
	void request_state3() {
		set_input("to_3");
	}





};




#endif