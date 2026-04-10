#pragma once

#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include "emp-zk/emp-zk-arith/ostriple.h"
#include "emp-zk/emp-zk.h"
#include "emp-agmpc/netmp.h"
#include "emp-agmpc/helper.h"
#include "prog_vole.h"
#include <cstdint>
#include <vector>

using namespace std;

// via PCG
template<int nP>
class MVVOLEarith { public:
	ProgVOLE<NetIO> * vole[nP+1];
	NetIOMP<nP> *io;
	ThreadPool * pool;
	int party;
	PRG prg;
	__uint128_t Delta;
	Hash hash;
	bool is_setup = false;
	block fixed_seed;

	MVVOLEarith(NetIOMP<nP>* io, ThreadPool * pool, int party) {
		this->io = io;
		this->pool = pool;
		this->party = party;

		Delta = (__uint128_t)0;

		PRG prg_tmp;
		prg_tmp.random_block(&fixed_seed, 1);


		if (party == 1){
			for(int i = 2; i <= nP; ++i){
				// acts as VOLE sender
				vole[i] = new ProgVOLE<NetIO>(2, 1, &io->get(i, false), fixed_seed);
				vole[i]->setup_with_seed();
			}
		}
		else{
			// act as VOLE receiver
			prg.random_data(&Delta, sizeof(__uint128_t));
    		Delta = Delta & ((__uint128_t)0xFFFFFFFFFFFFFFFFLL);
    		Delta = mod(Delta, pr);
			
			vole[1] = new ProgVOLE<NetIO>(1, 1, &io->get(1, false), fixed_seed);
			vole[1]->setup_with_seed(Delta);
		}

	}
	~MVVOLEarith() {
		if (party == 1){
			for (int i=2; i <= nP; ++i){
				delete vole[i];
			}
		}
		else{
			delete vole[1];
		}
	}



	void computeVOLE_semi(__uint128_t * output[nP+1], int length){
		if (party == 1){
			// obtains authenticated values and corresponding MAC tags
			// output[i] = u[i] || MAC_tag[i]
			vector<future<void>> res;
			res.reserve((size_t)(nP - 1));
			for (int i = 2; i<=nP; ++i){
				int party2 = i;
				res.push_back(pool->enqueue([this, output, length, party2]() {
					vole[party2]->extend_with_seed(output[party2], length);

				}));
			}
			joinNclean(res);
		}
		else{
			vole[1]->extend_with_seed(output[1], length);
		
		}
		io->flush();
		io->sync();
	}

	// supposed that the length is already increased by 1! I.e, length := length + 1
	// The last element of output will be thrown away
	void computeVOLE_mal(__uint128_t * output[nP+1],  int length){

		// semi-honest security
		computeVOLE_semi(output, length);

		// malicious security: by adding the following consistency check
		block seed = sampleRandom(io, &prg, pool, party); 
		vector<uint64_t> s((size_t)length);
		PRG prg2(&seed);
		prg2.random_data(s.data(), (int)(s.size() * sizeof(uint64_t)));

		for(int i = 0; i < length; ++i){
			s[i] = mod(s[i]);
		}
		
		uint64_t u;
		uint64_t tmpM[nP+1];
		
		if(party == 1){
			// compute u
			u = 0;
			for(int i = 2; i <= nP; ++i){
				tmpM[i] = 0;
			}

			uint64_t tmpu, tmpm;
			for (int i = 0; i <length; ++i){
				const uint64_t si = s[i];
				const uint64_t v2 = HIGH64(output[2][i]);
				tmpu = mult_mod(si, v2);
				u = add_mod(u, tmpu);

				for(int j = 2; j <= nP; ++j){	
					tmpm = mult_mod(LOW64(output[j][i]), si);
					tmpM[j] = add_mod(tmpM[j], tmpm);
				}
			}	
			// send u and tmpM to verifiers
			vector<future<void>> res2;
			res2.reserve((size_t)(nP - 1));
			for (int i = 2; i<=nP; ++i){
				int party2 = i;
				res2.push_back(pool->enqueue([this, u, tmpM, party2]() {
					io->send_data(party2, &u, sizeof(uint64_t));
					io->send_data(party2, &tmpM[party2], sizeof(uint64_t));
					io->flush(party2);
				}));
			}
			joinNclean(res2);
			
		}
		else{
			// compute check_key for party 1
			uint64_t check_key = 0;
			int i = 0;
			for (; i + 1 < length; i += 2) {
				uint64_t a[2] = { mod(output[1][i], pr), mod(output[1][i + 1], pr) };
				uint64_t b[2] = { s[i], s[i + 1] };
				uint64_t r[2];
				mult_mod_bch2(r, a, b);
				check_key = add_mod(check_key, add_mod(r[0], r[1]));
			}
			for (; i < length; ++i) {
				uint64_t tmpk = mod(output[1][i], pr);
				tmpk = mult_mod(s[i], tmpk);
				check_key = add_mod(check_key, tmpk);
			}

			// receive u and tmpM from the dealer
			io->recv_data(1, &u, sizeof(uint64_t));
			io->recv_data(1, &tmpM[1], sizeof(uint64_t));

			
			
			uint64_t check_mac = mult_mod(u, Delta);
			check_mac = add_mod(check_mac, check_key);

			if (check_mac != tmpM[1])
				cout << "Consistency check for MV-sVOLE fails!" << endl;
		}
	}

	/* debug functions */
	void debug_VOLE(__uint128_t * output[nP+1],  int length){
		if (party == 1){

			// check the authenticated values
			for(int i = 0; i < length; ++i){
				uint64_t value_baseline =  HIGH64(output[2][i]);
				for(int j = 3; j <= nP; ++j){
					uint64_t value_to_compare = HIGH64(output[j][i]);
					if(value_baseline != value_to_compare){
						cout << "Inconsistent values!" << endl;
					}
				}
			}


			for(int i = 2; i <= nP; ++i){
				io->send_data(i, output[i], length * sizeof(__uint128_t));
			}
		}
		else{
			__uint128_t * recv_output = new __uint128_t[length];
			io->recv_data(1, recv_output, length * sizeof(__uint128_t));
			for(int i = 0; i < length; ++i){
				__uint128_t tmp = mod(Delta * HIGH64(recv_output[i]), pr);
				tmp = mod(tmp + output[1][i], pr);
				if(tmp != (recv_output[i] & 0xFFFFFFFFFFFFFFFFLL)){
					cout << "VOLE fails" << endl;
				}
			}
		}
	}
	
	
};
