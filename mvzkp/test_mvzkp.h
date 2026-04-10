#pragma once

#include "mvzkp_bool_pcg.h"
#include "mvzkp_bool_iknp.h"
#include "mvzkp_arith.h"

using namespace std;
using namespace emp;

static const int MVZKP_BENCH_REPEATS = 1;

template<int nP>
int communication(NetIOMP<nP> * ios[2]) {
	return ios[0]->count() + ios[1]->count();
}

template<int nP>
void mvzkp_bool_pcg_bench_once(int party, NetIOMP<nP> * ios[2], ThreadPool * pool, string filename) {
	BristolFashion cf(filename.c_str());
	bool *in = new bool[cf.num_input]; bool *out = new bool[cf.num_output];
	memset(in, 0, cf.num_input*sizeof(bool));

	double sum_prep_ms = 0, sum_onl_ms = 0;
	double sum_prep_mb = 0, sum_onl_mb = 0;

	for (int iter = 0; iter < MVZKP_BENCH_REPEATS; ++iter) {
		int64_t c0 = communication<nP>(ios);

		auto start = clock_start();
		OneRound_SIF_Bool_PCG<nP>* sif = new OneRound_SIF_Bool_PCG<nP>(ios, pool, party, &cf);
		(void)time_from(start);
		ios[0]->sync();
		ios[1]->sync();

		start = clock_start();
		sif->Preprocess();
		double t_prep = time_from(start);
		int64_t c_prep = communication<nP>(ios);
		ios[0]->flush();
		ios[1]->flush();

		start = clock_start();
		sif->online(in, out);
		ios[0]->flush();
		ios[1]->flush();
		double t_onl = time_from(start);
		int64_t c_end = communication<nP>(ios);

		sum_prep_ms += t_prep / 1000.0;
		sum_onl_ms += t_onl / 1000.0;
		sum_prep_mb += (c_prep - c0) / 1000.0 / 1000.0;
		sum_onl_mb += (c_end - c_prep) / 1000.0 / 1000.0;

		delete sif;
	}

	delete[] in;
	delete[] out;

	double n = static_cast<double>(MVZKP_BENCH_REPEATS);
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / preprocessing time / "
		<< sum_prep_ms / n << " /ms\n" << flush;
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / preprocessing comm / "
		<< sum_prep_mb / n << " /MB\n" << flush;
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / online time / "
		<< sum_onl_ms / n << " /ms\n" << flush;
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / online comm / "
		<< sum_onl_mb / n << " /MB\n" << flush;
}

template<int nP>
void mvzkp_bool_iknp_bench_once(int party, NetIOMP<nP> * ios[2], ThreadPool * pool, string filename) {
	BristolFashion cf(filename.c_str());
	bool *in = new bool[cf.num_input]; bool *out = new bool[cf.num_output];
	memset(in, 0, cf.num_input*sizeof(bool));

	double sum_prep_ms = 0, sum_onl_ms = 0;
	double sum_prep_mb = 0, sum_onl_mb = 0;

	for (int iter = 0; iter < MVZKP_BENCH_REPEATS; ++iter) {
		int64_t c0 = communication<nP>(ios);

		auto start = clock_start();
		OneRound_SIF_Bool_IKNP<nP>* sif = new OneRound_SIF_Bool_IKNP<nP>(ios, pool, party, &cf);
		(void)time_from(start);
		ios[0]->sync();
		ios[1]->sync();

		start = clock_start();
		sif->Preprocess();
		double t_prep = time_from(start);
		int64_t c_prep = communication<nP>(ios);
		ios[0]->flush();
		ios[1]->flush();

		start = clock_start();
		sif->online(in, out);
		ios[0]->flush();
		ios[1]->flush();
		double t_onl = time_from(start);
		int64_t c_end = communication<nP>(ios);

		sum_prep_ms += t_prep / 1000.0;
		sum_onl_ms += t_onl / 1000.0;
		sum_prep_mb += (c_prep - c0) / 1000.0 / 1000.0;
		sum_onl_mb += (c_end - c_prep) / 1000.0 / 1000.0;

		delete sif;
	}

	delete[] in;
	delete[] out;

	double n = static_cast<double>(MVZKP_BENCH_REPEATS);
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / preprocessing time / "
		<< sum_prep_ms / n << " /ms\n" << flush;
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / preprocessing comm / "
		<< sum_prep_mb / n << " /MB\n" << flush;
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / online time / "
		<< sum_onl_ms / n << " /ms\n" << flush;
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / online comm / "
		<< sum_onl_mb / n << " /MB\n" << flush;
}

template<int nP>
void mvzkp_arith_bench_once(int party, NetIOMP<nP> * ios[2], ThreadPool * pool, string filename) {
	BristolFashion cf(filename.c_str());

	uint64_t *in = new uint64_t[cf.num_input]; uint64_t *out = new uint64_t[cf.num_output];
	memset(in, 0, cf.num_input*sizeof(uint64_t));
	//in[0] = 10;
    //in[1] = 20;
    //in[2] = 30;
    //in[3] = 40;

	double sum_prep_ms = 0, sum_onl_ms = 0;
	double sum_prep_mb = 0, sum_onl_mb = 0;

	for (int iter = 0; iter < MVZKP_BENCH_REPEATS; ++iter) {
		int64_t c0 = communication<nP>(ios);

		auto start = clock_start();
		OneRound_SIF_Arith<nP>* sif = new OneRound_SIF_Arith<nP>(ios, pool, party, &cf);
		(void)time_from(start);
		ios[0]->sync();
		ios[1]->sync();

		start = clock_start();
		sif->Preprocess();
		double t_prep = time_from(start);
		int64_t c_prep = communication<nP>(ios);
		ios[0]->flush();
		ios[1]->flush();

		start = clock_start();
		sif->online(in, out);
		ios[0]->flush();
		ios[1]->flush();
		double t_onl = time_from(start);
		int64_t c_end = communication<nP>(ios);

		sum_prep_ms += t_prep / 1000.0;
		sum_onl_ms += t_onl / 1000.0;
		sum_prep_mb += (c_prep - c0) / 1000.0 / 1000.0;
		sum_onl_mb += (c_end - c_prep) / 1000.0 / 1000.0;

		delete sif;
	}

	delete[] in;
	delete[] out;

	double n = static_cast<double>(MVZKP_BENCH_REPEATS);
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / preprocessing time / "
		<< sum_prep_ms / n << " /ms\n" << flush;
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / preprocessing comm / "
		<< sum_prep_mb / n << " /MB\n" << flush;
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / online time / "
		<< sum_onl_ms / n << " /ms\n" << flush;
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / online comm / "
		<< sum_onl_mb / n << " /MB\n" << flush;
}

template<int nP>
void bench_mpc_once(int party, NetIOMP<nP> * ios[2], ThreadPool * pool, string filename) {
	BristolFormat cf(filename.c_str());

	bool *in = new bool[cf.n1+cf.n2]; bool *out = new bool[cf.n3];
	memset(in, false, cf.n1+cf.n2);

	double sum_prep_ms = 0, sum_onl_ms = 0;
	double sum_prep_mb = 0, sum_onl_mb = 0;

	for (int iter = 0; iter < MVZKP_BENCH_REPEATS; ++iter) {
		int64_t c0 = communication<nP>(ios);

		auto start = clock_start();
		CMPC<nP>* mpc = new CMPC<nP>(ios, pool, party, &cf);
		ios[0]->flush();
		ios[1]->flush();
		(void)time_from(start);
		ios[0]->sync();
		ios[1]->sync();

		start = clock_start();
		mpc->function_independent();

		mpc->function_dependent();
		ios[0]->flush();
		ios[1]->flush();
		double t_prep = time_from(start);
		int64_t c_prep = communication<nP>(ios);
		ios[0]->flush();
		ios[1]->flush();

		start = clock_start();
		mpc->online(in, out);
		ios[0]->flush();
		ios[1]->flush();
		double t_onl = time_from(start);
		int64_t c_end = communication<nP>(ios);

		sum_prep_ms += t_prep / 1000.0;
		sum_onl_ms += t_onl / 1000.0;
		sum_prep_mb += (c_prep - c0) / 1000.0 / 1000.0;
		sum_onl_mb += (c_end - c_prep) / 1000.0 / 1000.0;

		delete mpc;
	}

	delete[] in;
	delete[] out;

	double n = static_cast<double>(MVZKP_BENCH_REPEATS);
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / preprocessing time / "
		<< sum_prep_ms / n << " /ms\n" << flush;
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / preprocessing comm / "
		<< sum_prep_mb / n << " /MB\n" << flush;
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / online time / "
		<< sum_onl_ms / n << " /ms\n" << flush;
	cout << "party " << party << "/ avg over " << MVZKP_BENCH_REPEATS << " runs / online comm / "
		<< sum_onl_mb / n << " /MB\n" << flush;
}
