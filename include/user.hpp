#ifndef USER_HPP
#define USER_HPP

// #include <iostream>
// #include <fstream>
#include <ctime>
#include <stack>
#include <string>
// #include <sstream>
// #include <vector>
// #include <cstdlib> // for exit()

using namespace std;

////////////////========================================================////////////////
////////////////////////========  Gestion temps	========////////////////////////////////

stack<clock_t> tictoc_stack;

void tic() {
	tictoc_stack.push(clock());
}

void toc() {
    double time =((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC;
    cout << "Time elapsed: " << time << endl;
    tictoc_stack.pop();
}

void toc(vector<double>& times) {
	double time =((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << time << endl;
	times.push_back(time);
	tictoc_stack.pop();
}



////////////////========================================================////////////////
////////////////////////========  Conversions	========////////////////////////////////

template <typename nbr>
string NbrToStr(nbr N){
	ostringstream strs;
	strs << N;
	string str = strs.str();
	return str;
}

int StrToInt(string str){
	stringstream i(str);
	int  N;
	i >> N;
	return N;
}

Real StrToReal(string str){
	stringstream i(str);
	Real  N;
	i >> N;
	return N;
}

Cplx StrToCplx(string str){
	stringstream i(str);
	Cplx  N;
	i >> N;
	return N;
}

////////////////========================================================////////////////
////////////////////////========    String splitting	========////////////////////////

vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
vector<string> split(const string &s, char delim) {
	vector<string> elems;
	split(s, delim, elems);
	return elems;
}


#endif