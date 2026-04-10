#pragma once

#include <cstdio>
#include <emp-tool/emp-tool.h>
#include "emp-agmpc/emp-agmpc.h"
#include "mvzkp/test_mvzkp.h"
#ifdef __linux__
#include <unistd.h>
#endif

using namespace std;
using namespace emp;

// BristolFashion opens the path with fopen and does not check for failure; resolve a real file first.
static string resolve_ands_circuit() {
	const char* rel[] = {
		"../circuits/bristol_fashion/arith_circuit.txt",
		"circuits/bristol_fashion/arith_circuit.txt",
	};
	for (const char* p : rel) {
		FILE* f = fopen(p, "r");
		if (f) {
			fclose(f);
			return string(p);
		}
	}
#ifdef __linux__
	char buf[4096];
	ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
	if (n > 0) {
		buf[n] = '\0';
		string dir(buf);
		for (int up = 0; up < 8; ++up) {
			size_t slash = dir.rfind('/');
			if (slash == string::npos)
				break;
			dir.resize(slash);
			string cand = dir + "/circuits/bristol_fashion/arith_circuit.txt";
			FILE* f = fopen(cand.c_str(), "r");
			if (f) {
				fclose(f);
				return cand;
			}
		}
	}
#endif
	return "";
}

// Must match the launch script: run_3 for nP=3, run_4 for nP=4 (run_3 only starts 3 processes → deadlock if nP>3).
const static int nP = 4;
int party, port;
int main(int argc, char** argv) {
	parse_party_and_port(argv, &party, &port);
	if(party > nP)return 0;
	string circuit_path = resolve_ands_circuit();
	if (circuit_path.empty()) {
		cerr << "party " << party << ": cannot open circuits/bristol_fashion/arith_circuit.txt "
		        "(add that file or run from the project build/ directory).\n";
		return 1;
	}
	NetIOMP<nP> io(party, port);
#ifdef LOCALHOST
	NetIOMP<nP> io2(party, port+2*(nP+1)*(nP+1)+1);
#else
	NetIOMP<nP> io2(party, port+2*(nP+1));
#endif
	NetIOMP<nP> *ios[2] = {&io, &io2};
	ThreadPool pool(2*(nP-1)+2);	

	if (party == 1){
		cout << "Evaluate Arithmetic Circuit" << endl;
	}
	mvzkp_arith_bench_once<nP>(party, ios, &pool, circuit_path);

	return 0;
}