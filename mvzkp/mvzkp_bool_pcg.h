#pragma once

#include "emp-agmpc/emp-agmpc.h"
#include "emp-tool/emp-tool.h"
#include "emp-tool/utils/f2k.h"
#include "mv-qsvole/mvqsvole_bool_pcg.h"
#include <cstdint>
#include <cstring>
#include <vector>

namespace mvzkp_bool_pcg_detail {
// Pack 64 bools (index 0 = LSB) into uint64_t; matches the original bit-by-bit loop.
static inline uint64_t pack_bool64_lsb(const bool *b) {
	uint64_t x = 0;
	for (int k = 0; k < 64; ++k)
		x |= (uint64_t)b[k] << k;
	return x;
}
} // namespace mvzkp_bool_pcg_detail

using namespace std;
using namespace emp;

//1-round SIF over binary field
template<int nP>
class OneRound_SIF_Bool_PCG{ public:
	
	MVSVOLEBOOL_PCG<nP>* fvolemp = nullptr;

	// single values
	block* mac[nP+1];
	block* key[nP+1];
	bool* value;

	// wires values
	block * wires_mac[nP+1];
	block * wires_key[nP+1];
	bool * wires;

	// vole wires values
	block* vole_wire_mac[nP + 1];
	bool* vole_wire_value;
	block* vole_wire_key;

	//vole quadratic values
	bool * vole_wire_y;
	block * vole_wire_y_key;
	block *vole_wire_y_mac[nP + 1];

	BristolFashion * cf;
	NetIOMP<nP> * io;
	int num_ands = 0, num_in;
	int party, total_pre, ssp;
	ThreadPool * pool;

	bool * mask_input;
	bool * mask_y;
	bool * d_alpha;
	bool * d_beta;
	bool * d_gamma;
	
	// for multiplication checks
	block * mt_A[nP+1];
	block * mt_B;
	block * y_A[nP+1];
	block * y_B;

	block Delta;
	block seed;
	PRP prp;
	PRG prg;

	// Reused in Preprocess() to avoid per-call heap traffic
	block *chi_buf = nullptr;
	block *ope_data_buf = nullptr;
	bool *ope_bool_buf = nullptr;

	// vole over extension fields (pair-wise)
	GaloisFieldPacking pack;
	block * value_ext[nP+1];
	block * mac_ext[nP+1];
	block * key_ext[nP+1];

	OneRound_SIF_Bool_PCG(NetIOMP<nP> * io[2], ThreadPool * pool, int party, BristolFashion * cf, bool * _delta = nullptr, int ssp=40) {
		this->party = party;
		this->io = io[0];
		// cf is the circuit
		this->cf = cf;
		this->pool = pool;
		
		for(int i = 0; i < cf->num_gate; ++i) {
			if (cf->gates[4*i+3] == AND_GATE)
				++num_ands;
		}
		num_in = cf->num_input;

		
		
		// 2*ssp for consistency check,  128*(np-1) for conversion
		total_pre = num_in + 2*num_ands + 2*ssp + (nP-1)*128+101;
		
		fvolemp = new MVSVOLEBOOL_PCG<nP>(io[1], pool, party);

		if (party != 1)
			Delta = fvolemp->Delta;

			for(int i  = 1; i <= nP; ++i) {
				key[i] = new block[total_pre];
				mac[i] = new block[total_pre];
				wires_key[i] = new block[cf->num_wire];
				wires_mac[i] = new block[cf->num_wire];
				mt_A[i] = new block[num_ands];
				y_A[i] = new block[num_ands*2];
				mac_ext[i] = new block[1];
				key_ext[i] = new block[1];
				value_ext[i] = new block[1];
				vole_wire_mac[i] = new block[cf->num_wire];
				vole_wire_y_mac[i] = new block[num_ands];
			}
			mt_B = new block[num_ands];
			y_B = new block[num_ands];
			value = new bool[total_pre];
			wires = new bool[cf->num_wire];
			mask_input = new bool[num_in];
			d_alpha = new bool[num_ands];
			d_beta = new bool[num_ands];
			d_gamma = new bool[num_ands];
			mask_y = new bool[num_ands];
			vole_wire_value = new bool[cf->num_wire];
			vole_wire_key = new block[cf->num_wire];
			vole_wire_y = new bool[num_ands];
			vole_wire_y_key = new block[num_ands];

			chi_buf = new block[num_ands];
			ope_data_buf = new block[128];
			ope_bool_buf = new bool[128];

	}
	~OneRound_SIF_Bool_PCG() {
		delete[] chi_buf;
		delete[] ope_data_buf;
		delete[] ope_bool_buf;
		delete[] mt_B;
		delete[] y_B;
		for(int i = 1; i <= nP; ++i) {
			delete[] key[i];
			delete[] mac[i];
			delete[] wires_key[i];
			delete[] wires_mac[i];
			delete[] mac_ext[i];
			delete[] key_ext[i];
			delete[] value_ext[i];
			delete[] vole_wire_mac[i];
			delete[] vole_wire_y_mac[i];
			delete[] mt_A[i];
			delete[] y_A[i];
		}
		delete[] value;
		delete[] wires;
		delete[] mask_input;
		delete[] d_alpha;
		delete[] d_beta;
		delete[] d_gamma;
		delete[] mask_y;
		delete[] vole_wire_value;
		delete[] vole_wire_key;
		delete[] vole_wire_y;
		delete[] vole_wire_y_key;
	}
    
	void Preprocess(){
		if (party == 1){
			fvolemp->computeVOLE_mal(mac, value, total_pre);
//generate extension vole
			for (int i = 2; i <= nP; ++i){
				int off_index = num_in + num_ands + (i-2)*128;
				memcpy(ope_data_buf, mac[i]+off_index, 128*sizeof(block));
				memcpy(ope_bool_buf, value+off_index, 128*sizeof(bool));

				uint64_t ch0 = mvzkp_bool_pcg_detail::pack_bool64_lsb(ope_bool_buf);
				uint64_t ch1 = mvzkp_bool_pcg_detail::pack_bool64_lsb(ope_bool_buf + 64);
				value_ext[i][0] = makeBlock(ch1, ch0);
				pack.packing(mac_ext[i], ope_data_buf);

			}

			for (int i = 0; i < num_in; ++i) {
				vole_wire_value[i] = value[i];
				for (int j = 2; j <= nP; ++j) {
					block *vwm = vole_wire_mac[j];
					block *m = mac[j];
					vwm[i] = m[i];
				}
			}
			int and_cnt = 0;
			for (int i = 0; i < cf->num_gate; ++i) {
				if (cf->gates[4*i+3] == AND_GATE) {
					const int g0 = cf->gates[4*i], g1 = cf->gates[4*i+1], g2 = cf->gates[4*i+2];
					vole_wire_value[g2] = value[num_in+and_cnt];
					bool u_alpha = vole_wire_value[g0];
					bool u_beta = vole_wire_value[g1];
					bool u_gamma = u_alpha & u_beta;
					vole_wire_y[and_cnt] = u_gamma;

                    //asign the mask for the quadratic value
					mask_y[and_cnt] = u_gamma ^ value[num_in+num_ands+and_cnt];
					for (int j = 2; j <= nP; ++j) {
						block *vwm = vole_wire_mac[j];
						block *mj = mac[j];
						block *yaj = y_A[j];
						block *vyym = vole_wire_y_mac[j];
						vwm[g2] = mj[num_in+and_cnt];
						vyym[and_cnt] = mj[num_in+num_ands+and_cnt];

						gfmul(vwm[g0], vwm[g1], &yaj[2*and_cnt]);
						yaj[2*and_cnt+1] = zero_block;
						if (u_alpha)
							yaj[2*and_cnt+1] = yaj[2*and_cnt+1] ^ vwm[g1];
						if (u_beta)
							yaj[2*and_cnt+1] = yaj[2*and_cnt+1] ^ vwm[g0];
						yaj[2*and_cnt+1] = yaj[2*and_cnt+1] ^ vyym[and_cnt];
					}
					and_cnt++;
				}
				else if (cf->gates[4*i+3] == XOR_GATE) {
					const int g0 = cf->gates[4*i], g1 = cf->gates[4*i+1], g2 = cf->gates[4*i+2];
					vole_wire_value[g2] = vole_wire_value[g0] ^ vole_wire_value[g1];
					for (int j = 2; j <= nP; ++j) {
						block *vwm = vole_wire_mac[j];
						vwm[g2] = vwm[g0] ^ vwm[g1];
					}
				}
				else {
					const int g0 = cf->gates[4*i], g2 = cf->gates[4*i+2];
					vole_wire_value[g2] = vole_wire_value[g0]^true;
					// compute the mac
					for (int j = 2; j <= nP; j++){
						block *vwm = vole_wire_mac[j];
						vwm[g2] = vwm[g0];
					}
			    }
		    }

			//debug_conversion();

			//check the quadratic relation
			seed = Hash::hash_for_block(mask_y, num_ands*sizeof(bool));
			PRG prg_chi;
			prg_chi.reseed(&seed);
			prg_chi.random_block(chi_buf, num_ands);

			for(int i = 2; i <= nP; ++i){
				block tmp_check;
				block com_y_A[2];
				com_y_A[0] = zero_block;
				com_y_A[1] = zero_block;
				block *yai = y_A[i];
				const block *chi = chi_buf;
				for(int j = 0; j < num_ands; ++j){
					gfmul(chi[j], yai[2*j], &tmp_check);
					com_y_A[0] = com_y_A[0] ^ tmp_check;
					gfmul(chi[j], yai[2*j+1], &tmp_check);
					com_y_A[1] = com_y_A[1] ^ tmp_check;
				}
				value_ext[i][0] = value_ext[i][0] ^ com_y_A[1];
				mac_ext[i][0] = mac_ext[i][0] ^ com_y_A[0];	
			}

			

			// send the proof (mask) to the verifiers
			vector<future<void>> res;
			res.reserve((size_t)(nP - 1));
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, party2]() {
					io->send_data(party2, mask_y, num_ands*sizeof(bool));
					io->send_data(party2, value_ext[party2], sizeof(block));
					io->send_data(party2, mac_ext[party2], sizeof(block));
					io->flush(party2);
				}));
			}
			joinNclean(res);

		}
		else{
			fvolemp->computeVOLE_mal(key, value, total_pre);
			
			int off_index = num_in + num_ands + (party-2)*128;
			memcpy(ope_data_buf, key[1]+off_index, 128*sizeof(block));
	
			pack.packing(key_ext[1], ope_data_buf);

			//debug_conversion();
			
			

			block recv_value_ext, recv_mac_ext;
			io->recv_data(1, mask_y, num_ands*sizeof(bool));
			io->recv_data(1, &recv_value_ext, sizeof(block));
			io->recv_data(1, &recv_mac_ext, sizeof(block));
			

			//input random wire
			block * const k1 = key[1];
			for (int i = 0; i < num_in; ++i) {
				vole_wire_key[i] = k1[i];
			}
			int and_cnt = 0;
			for (int i = 0; i < cf->num_gate; ++i) {
				if (cf->gates[4*i+3] == AND_GATE) {
					vole_wire_key[cf->gates[4*i+2]] = k1[num_in+and_cnt];
					vole_wire_y_key[and_cnt] = k1[num_in+num_ands+and_cnt];
					if(mask_y[and_cnt])
						vole_wire_y_key[and_cnt] = vole_wire_y_key[and_cnt] ^ Delta;

					//compute the values for quadratic multiplication relation
					block tmp_check;
					gfmul(vole_wire_key[cf->gates[4*i]], vole_wire_key[cf->gates[4*i+1]], &tmp_check);
					y_B[and_cnt] = tmp_check;
					gfmul(vole_wire_y_key[and_cnt], Delta, &tmp_check);
					y_B[and_cnt] = y_B[and_cnt] ^ tmp_check;
					
					and_cnt++;
				}
				else if (cf->gates[4*i+3] == XOR_GATE) {
					vole_wire_key[cf->gates[4*i+2]] = vole_wire_key[cf->gates[4*i]] ^ vole_wire_key[cf->gates[4*i+1]];
				}
				else {
					vole_wire_key[cf->gates[4*i+2]] = vole_wire_key[cf->gates[4*i]] ^ Delta;
				}
			}
			
			

			
//check the quadratic relation
			seed = Hash::hash_for_block(mask_y, num_ands*sizeof(bool));
			PRG prg_chi;
			prg_chi.reseed(&seed);
			prg_chi.random_block(chi_buf, num_ands);

			// compute the quadratic multiplication checks
			block tmp_check;
			block com_y_B = zero_block;
			block *yb = y_B;
			const block *chi = chi_buf;
			for(int i = 0; i < num_ands; ++i){
				gfmul(chi[i], yb[i], &tmp_check);
				com_y_B = com_y_B ^ tmp_check;
			}
			key_ext[1][0] = key_ext[1][0] ^ com_y_B;
			// check the multiplication gates
			gfmul(recv_value_ext, Delta, &tmp_check);
			tmp_check = tmp_check ^ recv_mac_ext;
			if(!cmpBlock(&tmp_check, &key_ext[1][0], 1))
				cout << "party " << party << " : quadratic multiplication checks fail" << endl;

		}

	}

	void online(bool * input, bool * output){
		block seed;

		if (party == 1){
			// Dealer generates the ''proof''

			// for each input wire
			
			for(int i = 0; i < num_in; ++i){
				mask_input[i] = input[i] ^ value[i];
				wires[i] = input[i];
				for (int j=2; j<=nP; ++j){
					block *wm = wires_mac[j];
					block *m = mac[j];
					wm[i] = m[i];
				}
			}
			
			// evaluate the circuit
			int counter_and = 0;
			for(int i = 0; i < cf->num_gate; ++i){
				const int *g = &cf->gates[4*i];
				if (g[3] == AND_GATE){

					bool w_alpha = wires[g[0]];
					bool w_beta = wires[g[1]];
					bool w_gamma = w_alpha & w_beta;

					// compute the output of AND gate
					wires[g[2]] = w_gamma; 

					// asign the mask for AND gate
					bool da = w_alpha ^ vole_wire_value[g[0]];
					bool db = w_beta ^ vole_wire_value[g[1]];
					bool dg = w_gamma ^ vole_wire_value[g[2]];
					d_alpha[counter_and] = da;
					d_beta[counter_and] = db;
					d_gamma[counter_and] = dg;
					
					// compute the values for multiplication checks
					for (int j = 2; j <= nP; j++){
						block *wm = wires_mac[j];
						block *vwm = vole_wire_mac[j];
						block *vyym = vole_wire_y_mac[j];
						block m_alpha = wm[g[0]];
						block m_beta = wm[g[1]];
						block m_gamma = vwm[g[2]];
						wm[g[2]] = m_gamma;

						block acc = zero_block;
						if (da)
							acc = acc ^ m_beta;
						if (db)
							acc = acc ^ m_alpha;
						acc = acc ^ vyym[counter_and];
						mt_A[j][counter_and] = acc ^ m_gamma;
					}

					counter_and++;
				}
				else if (g[3] == XOR_GATE){
					// compute the output of XOR gate
					wires[g[2]] = wires[g[0]] ^ wires[g[1]]; 

					// compute the mac
					for (int j = 2; j <= nP; j++){
						block *wm = wires_mac[j];
						wm[g[2]] = wm[g[0]] ^ wm[g[1]];
					}
				}
				else{
					wires[g[2]] = wires[g[0]] ^ true;

					// compute the mac
					for (int j = 2; j <= nP; j++){
						block *wm = wires_mac[j];
						wm[g[2]] = wm[g[0]];
					}
				}
			}

			// send the proof (mask) to the verifiers
			vector<future<void>> res;
			res.reserve((size_t)(nP - 1));
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, party2]() {
					block mt_A_or = zero_block;
					block *row = mt_A[party2];
					for (int k = 0; k < num_ands; ++k)
						mt_A_or = mt_A_or ^ row[k];
					block h = Hash::hash_for_block(&mt_A_or, sizeof(block));
					io->send_data(party2, mask_input, num_in*sizeof(bool));
					io->send_data(party2, d_alpha, num_ands*sizeof(bool));
					io->send_data(party2, d_beta, num_ands*sizeof(bool));
					io->send_data(party2, d_gamma, num_ands*sizeof(bool));
					//io->send_data(party2, wires, cf->num_wire*sizeof(bool));
					//io->send_data(party2, wires_mac[party2], cf->num_wire*sizeof(block));
					io->send_data(party2, &h, sizeof(block));
					io->flush(party2);
				}));
			}
			joinNclean(res);

		}
		else{
			// verifiers' protocol 

			// get the proof (mask) from the prover

			block mt_A_h;
			io->recv_data(1,mask_input,num_in*sizeof(bool));
			io->recv_data(1, d_alpha, num_ands*sizeof(bool));
			io->recv_data(1, d_beta, num_ands*sizeof(bool));
			io->recv_data(1, d_gamma, num_ands*sizeof(bool));
			//io->recv_data(1, wires, cf->num_wire*sizeof(bool));
			//io->recv_data(1, wires_mac[1], cf->num_wire*sizeof(block));
			io->recv_data(1, &mt_A_h, sizeof(block));
			

			// evaluate the circuit
			block * const wk = wires_key[1];
			block * const k1 = key[1];
			for (int i = 0; i < num_in; i++){
				block t = k1[i];
				if(mask_input[i])
					t = t ^ Delta;
				wk[i] = t;
			}

			int counter_and = 0;
			for(int i = 0; i < cf->num_gate; ++i){
				const int *g = &cf->gates[4*i];
				if (g[3] == AND_GATE){
					// compute the local mac key of the output of AND gate
					block wk_out = k1[num_in+counter_and];
					bool dg = d_gamma[counter_and];
					if(dg)
						wk_out = wk_out ^ Delta;
					wk[g[2]] = wk_out;

					bool da = d_alpha[counter_and];
					bool db = d_beta[counter_and];
					// compute the values for multiplication checks
					block tmp_check = zero_block;
					if(da && db)
						tmp_check = Delta;
					if(da)
					    tmp_check = tmp_check ^ wk[g[1]];
					if(db)
					    tmp_check =  tmp_check ^ wk[g[0]];
					tmp_check = tmp_check ^ vole_wire_y_key[counter_and];
					tmp_check = tmp_check ^ wk_out;
					mt_B[counter_and] = tmp_check;

					counter_and++;
				}
				else if (g[3] == XOR_GATE){

					// compute the local mac key
					wk[g[2]] = wk[g[0]] ^ wk[g[1]];
				}
				else{
					// compute the local mac key
					wk[g[2]] = wk[g[0]] ^ Delta;
				}
			}

			block mt_B_or = zero_block;
			block *mb = mt_B;
			for (int i = 0; i < num_ands; ++i)
				mt_B_or = mt_B_or ^ mb[i];
			block mt_B_h = Hash::hash_for_block(&mt_B_or, sizeof(block));
			if(!cmpBlock(&mt_B_h, &mt_A_h, 1))
				cout << "party " << party << " : multiplication checks fail" << endl;

		}

		
	}


	/*  debug function */
	void debug_conversion(){
		if (party == 1){
			for(int i = 2; i <= nP; ++i){
				io->send_data(i, value_ext[i], sizeof(block));
				io->send_data(i, mac_ext[i], sizeof(block));
			}
		}
		else{
			io->recv_data(1, value_ext[1], sizeof(block));
			io->recv_data(1, mac_ext[1], sizeof(block));


			block tmp;
			// tmp = value_ext[1][0] * Delta
			gfmul(value_ext[1][0], Delta, &tmp);
			tmp = tmp ^ mac_ext[1][0];
			if(!cmpBlock(&tmp, &key_ext[1][0], 1))
				cout << "conversion failed!" << endl;
		}
	}
	
	void debug_circuits_evaluation(){
		if(party == 1){
			for(int i = 2; i <= nP; ++i){
				io->send_data(i, wires, cf->num_wire*sizeof(bool));
				io->send_data(i, wires_mac[i], cf->num_wire*sizeof(block));
			}
		}
		else{
			bool * recv_wires = new bool[cf->num_wire];
			block * recv_wire_mac = new block[cf->num_wire];
			io->recv_data(1, recv_wires, cf->num_wire*sizeof(bool));
			io->recv_data(1, recv_wire_mac, cf->num_wire*sizeof(block));

			for(int i = 0; i < num_in; ++i){
				block tmp_check = recv_wire_mac[i];
				if(recv_wires[i])
					tmp_check = tmp_check ^ Delta;

				if(!cmpBlock(&tmp_check, &wires_key[1][i], 1))
					cout << "Input phase is wrong!" << endl;
			}

			for(int i = 0; i < cf->num_gate; ++i){
				if (cf->gates[4*i+3] == AND_GATE){
					block tmp_check = recv_wire_mac[cf->gates[4*i+2]];
					if(recv_wires[cf->gates[4*i+2]])
						tmp_check = tmp_check ^ Delta;
					
					if(!cmpBlock(&tmp_check, &wires_key[1][cf->gates[4*i+2]], 1)){
						cout << "AND gates are wrong!" << endl;
						// flag[0] = true;
						// first_wrong[0] = i;
					}
				}
				else if (cf->gates[4*i+3] == XOR_GATE){
					block tmp_check = recv_wire_mac[cf->gates[4*i+2]];
					if(recv_wires[cf->gates[4*i+2]])
						tmp_check = tmp_check ^ Delta;

					if(!cmpBlock(&tmp_check, &wires_key[1][cf->gates[4*i+2]], 1)){
						cout << "XOR gates are wrong!" << endl;
					}
						
				}
				else{
					block tmp_check = recv_wire_mac[cf->gates[4*i+2]];
					if(recv_wires[cf->gates[4*i+2]])
						tmp_check = tmp_check ^ Delta;

					if(!cmpBlock(&tmp_check, &wires_key[1][cf->gates[4*i+2]], 1)){
						cout << "NOT gates are wrong!" << endl;
					}
						
				}
			}
		}
	}


	void debug_process_vole_wire(){
		//cout<<"num_in: "<<num_in<<" "<<"num_gate: "<<cf->num_gate<<" "<<"num_wire: "<<cf->num_wire<<endl;
		if(party == 1){
			for(int i = 2; i <=nP; ++i){
				io->send_data(i, vole_wire_value, cf->num_wire*sizeof(bool));
				io->send_data(i, vole_wire_mac[i], cf->num_wire*sizeof(block));
				//io->send_data(i, vole_wire_y, num_ands*sizeof(uint64_t));
				//io->send_data(i, vole_wire_y_mac[i], num_ands*sizeof(uint64_t));
		
			}
		}
		else{
			bool * recv_vole_wire_value = new bool[cf->num_wire];
			block * recv_vole_wire_mac = new block[cf->num_wire];
			//uint64_t * recv_vole_wire_y = new uint64_t[num_ands];
			//uint64_t * recv_vole_wire_y_mac = new uint64_t[num_ands];
			io->recv_data(1, recv_vole_wire_value, cf->num_wire*sizeof(bool));
			io->recv_data(1, recv_vole_wire_mac, cf->num_wire*sizeof(block));
			//io->recv_data(1, recv_vole_wire_y, num_ands*sizeof(uint64_t));
			//io->recv_data(1, recv_vole_wire_y_mac, num_ands*sizeof(uint64_t));

		// Verifiers use recv_vole_wire_value (dealer broadcast); vole_wire_value[] is not filled here.
		for(int i = 0; i < num_in; ++i){
			block tmp = recv_vole_wire_mac[i];
			if(recv_vole_wire_value[i])
				tmp = tmp ^ Delta;
			if(!cmpBlock(&tmp, &vole_wire_key[i], 1))
				cout << "Party " << party << ": "<< i <<"-th input wire fails in the process" << endl;
		}

		for(int i = 0; i < cf->num_gate; ++i){
			if (cf->gates[4*i+3] == AND_GATE){
				block tmp = recv_vole_wire_mac[cf->gates[4*i+2]];
				if(recv_vole_wire_value[cf->gates[4*i+2]])
					tmp = tmp ^ Delta;
				if(!cmpBlock(&tmp, &vole_wire_key[cf->gates[4*i+2]], 1))
					cout << "Party " << party << ": "<< i <<"-th AND gate's output wire fails in the process" << endl;
			}
			else if (cf->gates[4*i+3] == XOR_GATE){
				block tmp = recv_vole_wire_mac[cf->gates[4*i+2]];
				if(recv_vole_wire_value[cf->gates[4*i+2]])
					tmp = tmp ^ Delta;
				if(!cmpBlock(&tmp, &vole_wire_key[cf->gates[4*i+2]], 1))
					cout << "Party " << party << ": "<< i <<"-th XOR gate's output wire fails in the process" << endl;
			}
			else if (cf->gates[4*i+3] == NOT_GATE){
				block tmp = recv_vole_wire_mac[cf->gates[4*i+2]];
				if(recv_vole_wire_value[cf->gates[4*i+2]])
					tmp = tmp ^ Delta;
				if(!cmpBlock(&tmp, &vole_wire_key[cf->gates[4*i+2]], 1))
					cout << "Party " << party << ": "<< i <<"-th NOT gate's output wire fails in the process" << endl;
			}
		}
		/*for(int i = 0; i < num_ands; ++i){
			uint64_t tmp = mult_mod(recv_vole_wire_y[i], Delta);
			tmp = add_mod(tmp, vole_wire_y_key[i]);
			if(tmp != recv_vole_wire_y_mac[i])
				cout << "Party " << party << ": "<< i <<"-th quadratic gate's output wire fails in the process" << endl;
		}*/
			delete[] recv_vole_wire_value;
			delete[] recv_vole_wire_mac;
		}
	}

	void debug_multiplication_check(){
		if(party == 1){
			for(int i = 2; i <= nP; ++i){
				io->send_data(i, mt_A[i], 2*num_ands*sizeof(block));
				io->send_data(i, value_ext[i], sizeof(block));
				io->send_data(i, mac_ext[i], sizeof(block));
			}
		}
		else{
			block * recv_mt_A = new block[2*num_ands];
			block recv_value_ext, recv_mac_ext;
			io->recv_data(1, recv_mt_A, 2*num_ands*sizeof(block));
			io->recv_data(1, &recv_value_ext, sizeof(block));
			io->recv_data(1, &recv_mac_ext, sizeof(block));

			// check extension VOLE
			block tmp_ext;
			gfmul(recv_value_ext, Delta, &tmp_ext);
			tmp_ext = tmp_ext ^ recv_mac_ext;
			if(!cmpBlock(&tmp_ext, &key_ext[1][0], 1))
				cout << "VOLE over extension field fails" << endl;

			for(int i = 0; i <num_ands; ++i){
				block tmp;
				gfmul(recv_mt_A[2*i+1], Delta, &tmp);
				tmp = tmp ^ recv_mt_A[2*i];
				if(!cmpBlock(&tmp, &mt_B[i], 1)){
					cout << "Each MT fails" << endl;
				}
			}

			PRG prg_chi;
			prg_chi.random_block(chi_buf, num_ands);
			block com_mt_A[2];
			com_mt_A[0] = zero_block;
			com_mt_A[1] = zero_block;
			block com_mt_b = zero_block;

			for(int i = 0; i < num_ands; ++i){
				block tmp;
				gfmul(recv_mt_A[2*i], chi_buf[i], &tmp);
				com_mt_A[0] = com_mt_A[0] ^ tmp;
				gfmul(recv_mt_A[2*i+1], chi_buf[i], &tmp);
				com_mt_A[1] = com_mt_A[1] ^ tmp;
				gfmul(mt_B[i], chi_buf[i], &tmp);
				com_mt_b = com_mt_b ^ tmp;
			}

			block tmp_check;
			gfmul(com_mt_A[1], Delta, &tmp_check);
			tmp_check = tmp_check ^ com_mt_A[0];
			if(!cmpBlock(&tmp_check, &com_mt_b, 1))
				cout << "Random combination of MT fails" << endl;

			delete[] recv_mt_A;
		}
	}
};