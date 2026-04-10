#pragma once

#include "emp-agmpc/emp-agmpc.h"
#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-vole/utility.h"
#include "emp-zk/emp-zk-arith/ostriple.h"
#include "mv-qsvole/mvqsvole_arith.h"
#include <cstdint>
#include <cstring>

using namespace std;
using namespace emp;

//1-round SIF over large prime field
template<int nP>
class OneRound_SIF_Arith{ public:
	
	MVVOLEarith<nP>* fvolemp = nullptr;
	
	
	// single values
	__uint128_t* mac[nP+1];
	__uint128_t* key[nP+1];


	// wires values
	uint64_t * wires;
	uint64_t * wires_mac[nP+1];
	uint64_t * wires_key;

	// vole wires values
	uint64_t *vole_wire_mac[nP + 1];
	uint64_t *vole_wire_value;
	uint64_t *vole_wire_key;

	//vole quadratic values
	uint64_t *vole_wire_y;
	uint64_t *vole_wire_y_key;
	uint64_t *vole_wire_y_mac[nP + 1];

	BristolFashion * cf;
	NetIOMP<nP> * io;
	int num_ands = 0, num_in;
	int party, total_pre, ssp;
	ThreadPool * pool;

	uint64_t * mask_input;
	uint64_t * d_gamma;
    uint64_t * mask_y;
	uint64_t * d_alpha;
	uint64_t * d_beta;
	
	// for multiplication checks
	uint64_t * mt_A[nP+1];
	uint64_t * mt_B;
	uint64_t * y_A[nP+1];
	uint64_t * y_B;

	__uint128_t Delta;
	PRP prp;
	PRG prg;

	// vole over extension fields (pair-wise)
	uint64_t vole_mac_ext[nP+1];
	uint64_t vole_value_ext[nP+1];
	uint64_t vole_key_ext;

    OneRound_SIF_Arith(NetIOMP<nP> * io[2], ThreadPool * pool, int party, BristolFashion * cf, bool * _delta = nullptr) {
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

		// 1 for consistency check, 1*nP for conversion
		total_pre = num_in + 2*num_ands + 1 + (nP-1);
		
		fvolemp = new MVVOLEarith<nP>(io[1], pool, party);


		Delta = fvolemp->Delta;


		for(int i  = 1; i <= nP; ++i) {
			key[i] = new __uint128_t[total_pre];
			mac[i] = new __uint128_t[total_pre];
			wires_mac[i] = new uint64_t[cf->num_wire];
			mt_A[i] = new uint64_t[num_ands];
			y_A[i] = new uint64_t[num_ands*2];
			vole_wire_mac[i] = new uint64_t[cf->num_wire];
			vole_wire_y_mac[i] = new uint64_t[num_ands];
		}
		mt_B = new uint64_t[num_ands];
		y_B = new uint64_t[num_ands];
		wires = new uint64_t[cf->num_wire];
		wires_key = new uint64_t[cf->num_wire];
		mask_input = new uint64_t[num_in];
		d_alpha = new uint64_t[num_ands];
		d_beta = new uint64_t[num_ands];
		d_gamma = new uint64_t[num_ands];
		mask_y = new uint64_t[num_ands];
		vole_wire_value = new uint64_t[cf->num_wire];
		vole_wire_key = new uint64_t[cf->num_wire];
		vole_wire_y = new uint64_t[num_ands];
		vole_wire_y_key = new uint64_t[num_ands];
	}

    ~OneRound_SIF_Arith() {
		for(int i = 1; i <= nP; ++i) {
			delete[] key[i];
			delete[] mac[i];	
			delete[] wires_mac[i];
			delete[] vole_wire_mac[i];
			delete[] vole_wire_y_mac[i];
			delete[] mt_A[i];
			delete[] y_A[i];
		}
		delete[] mt_B;
		delete[] y_B;
		delete[] mask_input;
		delete[] d_alpha;
		delete[] d_beta;
		delete[] d_gamma;
		delete[] mask_y;
		delete[] wires;
		delete[] wires_key;
		delete[] vole_wire_key;
		delete[] vole_wire_value;
		delete[] vole_wire_y;
		delete[] vole_wire_y_key;
	}
	


	void Preprocess(){


		block seed;
		uint64_t * chi = new uint64_t[num_ands];
		if (party == 1){
			fvolemp->computeVOLE_mal(mac, total_pre);

/************************************************************************* */


       for (int i = 0; i < num_in; ++i) {
	       vole_wire_value[i] = HIGH64(mac[2][i]);
		   for(int j = 2; j <= nP; ++j){
			   vole_wire_mac[j][i] = LOW64(mac[j][i]);
		   }
		}
        //cout<<"party: "<<party<<" vole_wire_value_wire_0_initial: "<<vole_wire_value[0]<<endl;

      int and_cnt = 0;
      // 按门遍历输出线
      for (int i = 0; i < cf->num_gate; ++i) {
	        int type = cf->gates[4*i+3];
	        int out_wire = cf->gates[4*i+2];
	     if (type == AND_GATE) {
		// VOLE 为随机VOLE，即mac[2][num_in+and_cnt]
		    int in0 = cf->gates[4*i+0];
		    int in1 = cf->gates[4*i+1];
		    vole_wire_value[out_wire] = HIGH64(mac[2][num_in + and_cnt]);
			uint64_t u_alpha = vole_wire_value[in0];
			uint64_t u_beta = vole_wire_value[in1];
			uint64_t u_gamma = mult_mod(u_alpha, u_beta); 
			//generate quadratic relation
			vole_wire_y[and_cnt] = u_gamma;

			//asign the mask for the quadratic value
			mask_y[and_cnt] = PR - u_gamma;
			mask_y[and_cnt] = add_mod(mask_y[and_cnt], HIGH64(mac[2][num_in + num_ands + and_cnt]));
			
			for(int j = 2; j <= nP; ++j){
				vole_wire_mac[j][out_wire] = LOW64(mac[j][num_in + and_cnt]);
				vole_wire_y_mac[j][and_cnt] = LOW64(mac[j][num_in + num_ands + and_cnt]);

				y_A[j][2*and_cnt] = mult_mod(vole_wire_mac[j][in0], vole_wire_mac[j][in1]);
				y_A[j][2*and_cnt+1] = PR - vole_wire_y_mac[j][and_cnt];
				y_A[j][2*and_cnt+1] = add_mod( mult_mod(u_alpha, vole_wire_mac[j][in1]) , y_A[j][2*and_cnt+1]);
				y_A[j][2*and_cnt+1] = add_mod( mult_mod(u_beta, vole_wire_mac[j][in0]) , y_A[j][2*and_cnt+1]);
			}
		    and_cnt++;
	    } else if (type == XOR_GATE) {
		   int in0 = cf->gates[4*i+0];
		   int in1 = cf->gates[4*i+1];
		   vole_wire_value[out_wire] = add_mod(vole_wire_value[in0], vole_wire_value[in1]);
		   for(int j = 2; j <= nP; ++j){
			   vole_wire_mac[j][out_wire] = add_mod(vole_wire_mac[j][in0], vole_wire_mac[j][in1]);
		   }
	    } else {
		// 对于直接连线等其它类型，直接转抄输入
		   int in0 = cf->gates[4*i+0];
		   vole_wire_value[out_wire] = vole_wire_value[in0];
		   for(int j = 2; j <= nP; ++j){
			   vole_wire_mac[j][out_wire] = vole_wire_mac[j][in0];
		   }
	}
}

/************************************************************************** */
	        //check quadratic relation
			//using seed to generate different chi's
			seed = Hash::hash_for_block(mask_y, num_ands * sizeof(uint64_t));
			PRG prg_chi;
			prg_chi.reseed(&seed);
			prg_chi.random_data(chi, num_ands * sizeof(uint64_t));
			for(int i = 0; i < num_ands; ++i){
				chi[i] = mod(chi[i]);
			}
			
			for(int i = 2; i <= nP; ++i){
				int off_index = num_in + 2*num_ands + i-2;
				vole_value_ext[i] = HIGH64(mac[i][off_index]);
				vole_mac_ext[i] = LOW64(mac[i][off_index]);
			}
			

			for(int j = 0; j < num_ands; ++j){
				uint64_t cj = chi[j];
				for(int i = 2; i <= nP; ++i){
					uint64_t tmp = mult_mod(y_A[i][2*j], cj);
					vole_mac_ext[i] = add_mod(vole_mac_ext[i], tmp);
					tmp = mult_mod(y_A[i][2*j+1], cj);
					vole_value_ext[i] = add_mod(vole_value_ext[i], tmp);
				}
			}

			// send the proof (mask) to the verifiers
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, party2]() {
					io->send_data(party2, mask_y, num_ands*sizeof(uint64_t));
					io->send_data(party2, &vole_value_ext[party2], sizeof(uint64_t));
					io->send_data(party2, &vole_mac_ext[party2], sizeof(uint64_t));
					io->flush(party2);
				}));
			}
			joinNclean(res);

			//fvolemp->debug_VOLE(mac, total_pre);
			//debug_process_vole_wire();
		}
		else{
			fvolemp->computeVOLE_mal(key, total_pre);
			uint64_t recv_vole_value_ext, recv_vole_mac_ext;
			io->recv_data(1, mask_y, num_ands * sizeof(uint64_t));
			io->recv_data(1, &recv_vole_value_ext, sizeof(uint64_t));
			io->recv_data(1, &recv_vole_mac_ext, sizeof(uint64_t));
			
/******************************************************** */
			// 输入线
			for (int i = 0; i < num_in; ++i) {
				vole_wire_key[i] = LOW64(key[1][i]);
			}

			int and_cnt = 0;
			// 按门遍历输出线
			for (int i = 0; i < cf->num_gate; ++i) {
				int type = cf->gates[4*i+3];
				int out_wire = cf->gates[4*i+2];
				if (type == AND_GATE) {
					int in0 = cf->gates[4*i+0];
					int in1 = cf->gates[4*i+1];
					vole_wire_key[out_wire] = LOW64(key[1][num_in + and_cnt]);
					uint64_t tmp = mult_mod(mask_y[and_cnt], Delta);
				    vole_wire_y_key[and_cnt] = add_mod(LOW64(key[1][num_in + num_ands + and_cnt]), tmp);
					uint64_t tmp_key = mult_mod(vole_wire_key[in0], vole_wire_key[in1]);

					//generate quadratic relation
					uint64_t tmp_key2 = mult_mod(vole_wire_y_key[and_cnt], Delta);
					y_B[and_cnt] = add_mod(tmp_key, tmp_key2);
					//cout<<"party: "<<party<<"and_key: "<<vole_wire_key[out_wire]<<endl;
					and_cnt++;
				} else if (type == XOR_GATE) {
					int in0 = cf->gates[4*i+0];
					int in1 = cf->gates[4*i+1];
					vole_wire_key[out_wire] = add_mod(vole_wire_key[in0], vole_wire_key[in1]);
				} else {
					int in0 = cf->gates[4*i+0];
					vole_wire_key[out_wire] = vole_wire_key[in0];
				}
			}
		
		    //check quadratic relation
			/*********************************************** */
			int off_index = num_in + 2*num_ands + party-2;
			vole_key_ext = LOW64(key[1][off_index]);

			// using seed to generate different chi's
			seed = Hash::hash_for_block(mask_y, num_ands*sizeof(uint64_t));
			PRG prg_chi;
			prg_chi.reseed(&seed);
			prg_chi.random_data(chi, num_ands*sizeof(uint64_t));
			for(int i = 0; i < num_ands; ++i){
				chi[i] = mod(chi[i]);
			}

			// compute the multiplication checks
			for (int i = 0; i < num_ands; ++i){
				uint64_t tmp = mult_mod(y_B[i], chi[i]);
				vole_key_ext = add_mod(tmp, vole_key_ext);  
			}

			// check the multiplication gates
			uint64_t tmp_check = mult_mod(recv_vole_value_ext, Delta);
			tmp_check = add_mod(vole_key_ext, tmp_check);
			if(tmp_check != recv_vole_mac_ext)
				cout << "party " << party << " : quadratic multiplication checks fail" << endl;
			
			//fvolemp->debug_VOLE(key, total_pre);
			//debug_process_vole_wire();
		}
		delete[] chi;
	}

	
	void online(uint64_t * input, uint64_t * output){
		
		block seed; 
		if (party == 1){
			// Prover generates the ''proof''

			// for each input wire
			for(int i = 0; i < num_in; ++i){
				wires[i] = input[i];
				uint64_t minus_input = PR - input[i];
				mask_input[i] = add_mod(vole_wire_value[i], minus_input);
				for (int j=2; j<=nP; ++j){
					wires_mac[j][i] = vole_wire_mac[j][i];
				}
			}
			
			
			// evaluate the circuit
			int counter_and = 0;
			for(int i = 0; i < cf->num_gate; ++i){
				if (cf->gates[4*i+3] == AND_GATE){
					// compute the output of AND gate
					int wa = cf->gates[4*i];
					int wb = cf->gates[4*i+1];
					int wout = cf->gates[4*i+2];
					uint64_t w_alpha = wires[wa];
					uint64_t w_beta = wires[wb];
					uint64_t w_gamma = mult_mod(w_alpha, w_beta); 
					//cout<<"gate: "<<i<<" "<<w_alpha<<" "<<w_beta<<" "<<w_gamma<<" "<<endl;
					wires[wout] = w_gamma;
					
					// asign the mask for AND gate
					d_alpha[counter_and] = PR - w_alpha;
					d_alpha[counter_and] = add_mod(d_alpha[counter_and], vole_wire_value[wa]);
					d_beta[counter_and] = PR - w_beta;
					d_beta[counter_and] = add_mod(d_beta[counter_and], vole_wire_value[wb]);
					d_gamma[counter_and] = PR - w_gamma;
					d_gamma[counter_and] = add_mod(d_gamma[counter_and], vole_wire_value[wout]);

					// compute the values for multiplication checks
					for (int j = 2; j <= nP; j++){
						uint64_t m_alpha = wires_mac[j][wa];
						uint64_t m_beta = wires_mac[j][wb];
						uint64_t m_gamma = vole_wire_mac[j][wout];

						wires_mac[j][wout] = m_gamma; 

						uint64_t term1 = mult_mod(d_beta[counter_and], m_alpha);
						uint64_t term2 = mult_mod(d_alpha[counter_and], m_beta);
						mt_A[j][counter_and] = add_mod(term1, term2);
						mt_A[j][counter_and] = add_mod(mt_A[j][counter_and], m_gamma);
						mt_A[j][counter_and] = add_mod(mt_A[j][counter_and], (PR - vole_wire_y_mac[j][counter_and]));
						//cout<<"party: "<<party<<" mt_A["<<j<<"]["<<counter_and<<"] = "<<mt_A[j][counter_and]<<endl;
					}
					counter_and++;
				}
				else if (cf->gates[4*i+3] == XOR_GATE){
					// compute the output of XOR gate
					wires[cf->gates[4*i+2]] = add_mod(wires[cf->gates[4*i]], wires[cf->gates[4*i+1]]);
					//cout<<"gate: "<<i<<" "<<wires[cf->gates[4*i]]<<" "<<wires[cf->gates[4*i+1]]<<" "<<wires[cf->gates[4*i+2]]<<" "<<endl;
					for (int j = 2; j <= nP; ++j){
						wires_mac[j][cf->gates[4*i+2]] = add_mod(wires_mac[j][cf->gates[4*i]], wires_mac[j][cf->gates[4*i+1]]);
						  
					}
				}
			}
		
			// send the proof (mask) to the verifiers (per-party mt_A hash folded by OR)
			vector<future<void>> res;
			for (int party2 = 2; party2 <= nP; ++party2) {
				uint64_t mt_A_or = 0;
				for (int k = 0; k < num_ands; ++k)
					mt_A_or |= mt_A[party2][k];
				block mt_A_h = Hash::hash_for_block(&mt_A_or, sizeof(uint64_t));
				res.push_back(pool->enqueue([this, party2, mt_A_h]() {
					io->send_data(party2, mask_input, num_in*sizeof(uint64_t));
					io->send_data(party2, d_alpha, num_ands*sizeof(uint64_t));
					io->send_data(party2, d_beta, num_ands*sizeof(uint64_t));
					io->send_data(party2, d_gamma, num_ands*sizeof(uint64_t));
					io->send_data(party2, wires, cf->num_wire*sizeof(uint64_t));
					io->send_data(party2, wires_mac[party2], cf->num_wire*sizeof(uint64_t));
					io->send_data(party2, &mt_A_h, sizeof(block));
					io->flush(party2);
				}));
			}
			joinNclean(res);
			//debug_circuit_evaluation();

		}
		else{
			
			// verifiers' protocol 
			uint64_t * rev_wires = new uint64_t[cf->num_wire];
			uint64_t * rev_wires_mac = new uint64_t[cf->num_wire];// get the proof (mask) from the dealer
			block mt_A_h;
			io->recv_data(1, mask_input, num_in*sizeof(uint64_t));
			io->recv_data(1, d_alpha, num_ands*sizeof(uint64_t));
			io->recv_data(1, d_beta, num_ands*sizeof(uint64_t));
			io->recv_data(1, d_gamma, num_ands*sizeof(uint64_t));
			io->recv_data(1, rev_wires, cf->num_wire*sizeof(uint64_t));
			io->recv_data(1, rev_wires_mac, cf->num_wire*sizeof(uint64_t));
			io->recv_data(1, &mt_A_h, sizeof(block));

			// evaluate the circuit
			for (int i = 0; i < num_in; ++i){
				uint64_t tmp = mult_mod(mask_input[i], Delta);
				wires_key[i] =  add_mod(vole_wire_key[i], tmp);
				if(wires_key[i] != add_mod(rev_wires_mac[i],(PR - mult_mod(Delta, rev_wires[i]))))
					cout << "Party " << party << ": "<< i <<"-th input wire fails online" << endl;
			}

			int counter_and = 0;
			for(int i = 0; i < cf->num_gate; ++i){
				if (cf->gates[4*i+3] == AND_GATE){
					// compute the multiplicaiton gate
					int wa = cf->gates[4*i];
					int wb = cf->gates[4*i+1];
					int wout = cf->gates[4*i+2];

					uint64_t tmp = mult_mod(d_gamma[counter_and], Delta);
					wires_key[wout] = add_mod(vole_wire_key[wout], tmp);
					tmp = mult_mod(rev_wires[wout], Delta);
					tmp = add_mod(tmp, wires_key[wout]);
					if(tmp != rev_wires_mac[wout])
						cout << "Party " << party << ": "<< i <<"-th multiplication gate's output wire fails online" << endl;

					// compute the values for multiplication checks
					uint64_t tmp_key = mult_mod(vole_wire_key[wa], d_beta[counter_and]);
					uint64_t tmp_key1 = mult_mod(d_alpha[counter_and], vole_wire_key[wb]);
					tmp = add_mod(tmp_key, tmp_key1);
					tmp = add_mod(tmp, (PR - vole_wire_y_key[counter_and]));
					uint64_t tmp_key2 = mult_mod(mult_mod(d_alpha[counter_and], d_beta[counter_and]), Delta);
					tmp = add_mod(tmp, tmp_key2);
					mt_B[counter_and] = add_mod(tmp, wires_key[wout]);
					//cout<<"party: "<<party<<" mt_B["<<counter_and<<"] = "<<mt_B[counter_and]<<endl;
					counter_and++;
				}
				else if (cf->gates[4*i+3] == XOR_GATE){
					// compute the add gate
					wires_key[cf->gates[4*i+2]] = add_mod(wires_key[cf->gates[4*i]], wires_key[cf->gates[4*i+1]]);
				}
			}
			

			uint64_t mt_B_or = 0;
			for (int i = 0; i < num_ands; ++i) {
					//compute OR operation
				mt_B_or |= mt_B[i];
			}
			block mt_B_h = Hash::hash_for_block(&mt_B_or, sizeof(uint64_t));
			if (!cmpBlock(&mt_B_h, &mt_A_h, 1))
				cout << "party " << party << " : multiplication checks fail_new" << endl;
				// 可以将mt_A_hash进行后续使用，比如打印或存储
				// cout << "mt_A[" << j << "] OR hash: " << *(uint64_t*)&mt_A_hash << endl;	
			//debug_circuit_evaluation();
			delete[] rev_wires;
			delete[] rev_wires_mac;
		}
	}

	
/* debug functions */
   void debug_extension_vole(){
		if(party == 1){
			for(int i = 2; i <= nP; ++i){
				io->send_data(i, &vole_value_ext[i], sizeof(uint64_t));
				io->send_data(i, &vole_mac_ext[i], sizeof(uint64_t));
			}
		}
		else{
			uint64_t recv_vole_value, recv_vole_mac;
			io->recv_data(1, &recv_vole_value, sizeof(uint64_t));
			io->recv_data(1, &recv_vole_mac, sizeof(uint64_t));

			uint64_t tmp = mult_mod(recv_vole_value, Delta);
			tmp = add_mod(tmp, vole_key_ext);
			if(tmp != recv_vole_mac)
				cout << "Extension VOLE fails" << endl;
		}
	}

	/*************************************************** */
	void debug_process_vole_wire(){
		//cout<<"num_in: "<<num_in<<" "<<"num_gate: "<<cf->num_gate<<" "<<"num_wire: "<<cf->num_wire<<endl;
		if(party == 1){
			for(int i = 2; i <=nP; ++i){
				io->send_data(i, vole_wire_value, cf->num_wire*sizeof(uint64_t));
				io->send_data(i, vole_wire_mac[i], cf->num_wire*sizeof(uint64_t));
				io->send_data(i, vole_wire_y, num_ands*sizeof(uint64_t));
				io->send_data(i, vole_wire_y_mac[i], num_ands*sizeof(uint64_t));
		
			}
		}
		else{
			uint64_t * recv_vole_wire_value = new uint64_t[cf->num_wire];
			uint64_t * recv_vole_wire_mac = new uint64_t[cf->num_wire];
			uint64_t * recv_vole_wire_y = new uint64_t[num_ands];
			uint64_t * recv_vole_wire_y_mac = new uint64_t[num_ands];
			io->recv_data(1, recv_vole_wire_value, cf->num_wire*sizeof(uint64_t));
			io->recv_data(1, recv_vole_wire_mac, cf->num_wire*sizeof(uint64_t));
			io->recv_data(1, recv_vole_wire_y, num_ands*sizeof(uint64_t));
			io->recv_data(1, recv_vole_wire_y_mac, num_ands*sizeof(uint64_t));

		// check the input phase
		for(int i = 0; i < num_in; ++i){
			uint64_t tmp = mult_mod(recv_vole_wire_value[i], Delta);
			tmp = add_mod(tmp, vole_wire_key[i]);
			if(tmp != recv_vole_wire_mac[i])
				cout << "Party " << party << ": "<< i <<"-th input wire fails in the process" << endl;
		}

		for(int i = 0; i < cf->num_gate; ++i){
			if (cf->gates[4*i+3] == AND_GATE){
				//cout<<"party: "<<party<<" and_key "<<vole_wire_key[cf->gates[4*i+2]]<<endl;
				uint64_t tmp = mult_mod(recv_vole_wire_value[cf->gates[4*i+2]], Delta);
				tmp = add_mod(tmp, vole_wire_key[cf->gates[4*i+2]]);
				if(tmp != recv_vole_wire_mac[cf->gates[4*i+2]])
					cout << "Party " << party << ": "<< i <<"-th multiplication gate's output wire fails in the process" << endl;
			}
			else if (cf->gates[4*i+3] == XOR_GATE){
				uint64_t tmp = mult_mod(recv_vole_wire_value[cf->gates[4*i+2]], Delta);
				tmp = add_mod(tmp, vole_wire_key[cf->gates[4*i+2]]);
				if(tmp != recv_vole_wire_mac[cf->gates[4*i+2]])
					cout << "Party " << party << ": "<< i <<"-th addition gate's output wire fails in the process" << endl;
			}
		}
		for(int i = 0; i < num_ands; ++i){
			uint64_t tmp = mult_mod(recv_vole_wire_y[i], Delta);
			tmp = add_mod(tmp, vole_wire_y_key[i]);
			if(tmp != recv_vole_wire_y_mac[i])
				cout << "Party " << party << ": "<< i <<"-th quadratic gate's output wire fails in the process" << endl;
		}
		}
	}

	/************************************************************************* */
	void debug_circuit_evaluation(){
		if(party == 1){
			for(int i = 2; i <= nP; ++i){
				io->send_data(i, wires, cf->num_wire*sizeof(uint64_t));
				io->send_data(i, wires_mac[i], cf->num_wire*sizeof(uint64_t));
			}
		}
		else{
			uint64_t * recv_wire = new uint64_t[cf->num_wire];
			uint64_t * recv_wire_mac = new uint64_t[cf->num_wire];
			io->recv_data(1, recv_wire, cf->num_wire*sizeof(uint64_t));
			io->recv_data(1, recv_wire_mac, cf->num_wire*sizeof(uint64_t));

			// check the input phase
			for(int i = 0; i < num_in; ++i){
				uint64_t tmp = mult_mod(recv_wire[i], Delta);
				tmp = add_mod(tmp, wires_key[i]);
				if(tmp != recv_wire_mac[i]){
					cout << "Party " << party << ": "<< i <<"-th input wire fails" << endl;
				}
			}

			// check the multiplication gates
			for(int i = 0; i < cf->num_gate; ++i){
				if (cf->gates[4*i+3] == AND_GATE){
					uint64_t tmp = mult_mod(recv_wire[cf->gates[4*i+2]], Delta);
					tmp = add_mod(tmp, wires_key[cf->gates[4*i+2]]);
					if(tmp != recv_wire_mac[cf->gates[4*i+2]])
						cout << "Party " << party << ": "<< i <<"-th multiplication gate's output wire fails" << endl;
				}
				else if (cf->gates[4*i+3] == XOR_GATE){
					uint64_t tmp = mult_mod(recv_wire[cf->gates[4*i+2]], Delta);
					tmp = add_mod(tmp, wires_key[cf->gates[4*i+2]]);
					if(tmp != recv_wire_mac[cf->gates[4*i+2]])
						cout << "Party " << party << ": "<< i <<"-th addition gate's output wire fails" << endl;
				}
			}


		}
	}
};