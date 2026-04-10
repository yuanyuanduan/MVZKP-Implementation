#pragma once

#include "emp-agmpc/emp-agmpc.h"
#include "emp-tool/emp-tool.h"
#include "emp-tool/utils/f2k.h"
#include "mv-qsvole/mvqsvole_bool_iknp.h"
#include <cstdint>
#include <cstring>
#include <vector>

using namespace std;
using namespace emp;

//1-round SIF over binary field
template<int nP>
class OneRound_SIF_Bool_IKNP{ public:
	
	MVSVOLEBOOL_IKNP<nP>* fvolemp = nullptr;

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

	// vole over extension fields (pair-wise)
	GaloisFieldPacking pack;
	block * value_ext[nP+1];
	block * mac_ext[nP+1];
	block * key_ext[nP+1];

	OneRound_SIF_Bool_IKNP(NetIOMP<nP> * io[2], ThreadPool * pool, int party, BristolFashion * cf, bool * _delta = nullptr, int ssp=40) {
		this->party = party;
		this->io = io[0];
		// cf is the circuit
		this->cf = cf;
		this->pool = pool;
		
		{
			const int *g = cf->gates.data();
			for(int i = 0; i < cf->num_gate; ++i, g += 4)
				if (g[3] == AND_GATE)
					++num_ands;
		}
		num_in = cf->num_input;

		//cout << "num_in: " << num_in << " num_ands: " << num_ands <<"wire_num: " << cf->num_wire << endl;
		// 2*ssp for consistency check,  128*(np-1) for conversion
		total_pre = num_in + 2*num_ands + 2*ssp + (nP-1)*128+101;
		
		fvolemp = new MVSVOLEBOOL_IKNP<nP>(io[1], pool, party);


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
	}
	~OneRound_SIF_Bool_IKNP() {
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
		delete[] mt_B;
		delete[] y_B;
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
		if (party == 1) {
			fvolemp->computeVOLE_mal(mac, value, total_pre);
			bool_iknp_preprocess_prover_extension_vole<nP>(
				BoolIknpPreprocessProverExtCtx<nP>(num_in, num_ands, mac, value, pack, value_ext, mac_ext));
			bool_iknp_preprocess_prover_mac_on_circuit<nP>(BoolIknpPreprocessProverMacCtx<nP>(cf, num_in, num_ands,
				value, mac, vole_wire_value, vole_wire_mac, vole_wire_y, mask_y, vole_wire_y_mac, y_A));
			bool_iknp_preprocess_prover_quadratic_fold_and_send<nP>(
				BoolIknpPreprocessProverQuadCtx<nP>(io, pool, num_ands, mask_y, y_A, value_ext, mac_ext, &seed));
		} else {
			fvolemp->computeVOLE_mal(key, value, total_pre);
			Delta = fvolemp->Delta;
			bool_iknp_preprocess_verifier_extension_recv_and_check<nP>(BoolIknpPreprocessVerifierCtx<nP>(io, party,
				num_in, num_ands, Delta, cf, key[1], pack, key_ext, mask_y, vole_wire_key, vole_wire_y_key, y_B,
				&seed));
		}
	}

	void online(bool * input, bool * output){
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
			const int *gate_row = cf->gates.data();
			int counter_and = 0;
			for(int i = 0; i < cf->num_gate; ++i, gate_row += 4){
				const int wa = gate_row[0], wb = gate_row[1], wo = gate_row[2];
				const int gtype = gate_row[3];
				if (gtype == AND_GATE){

					bool w_alpha = wires[wa];
					bool w_beta = wires[wb];
					bool w_gamma = w_alpha & w_beta;

					wires[wo] = w_gamma;

					d_alpha[counter_and] = w_alpha ^ vole_wire_value[wa];
					d_beta[counter_and] = w_beta ^ vole_wire_value[wb];
					d_gamma[counter_and] = w_gamma ^ vole_wire_value[wo];

					const bool da = d_alpha[counter_and];
					const bool db = d_beta[counter_and];
					for (int j = 2; j <= nP; j++){
						block *wm = wires_mac[j];
						block *vwm = vole_wire_mac[j];
						block *vyym = vole_wire_y_mac[j];
						block m_alpha = wm[wa];
						block m_beta = wm[wb];
						block m_gamma = vwm[wo];
						wm[wo] = m_gamma;

						block t = vyym[counter_and] ^ m_gamma;
						if (da) t = t ^ m_beta;
						if (db) t = t ^ m_alpha;
						mt_A[j][counter_and] = t;
					}

					counter_and++;
				}
				else if (gtype == XOR_GATE){
					wires[wo] = wires[wa] ^ wires[wb];
					for (int j = 2; j <= nP; j++){
						block *wm = wires_mac[j];
						wm[wo] = wm[wa] ^ wm[wb];
					}
				}
				else{
					wires[wo] = wires[wa] ^ true;
					for (int j = 2; j <= nP; j++){
						block *wm = wires_mac[j];
						wm[wo] = wm[wa];
					}
				}
			}

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
			

			block * const wk = wires_key[1];
			block * const k1 = key[1];
			for (int i = 0; i < num_in; i++){
				block t = k1[i];
				if(mask_input[i])
					t = t ^ Delta;
				wk[i] = t;
			}

			const int *gate_row = cf->gates.data();
			int counter_and = 0;
			for(int i = 0; i < cf->num_gate; ++i, gate_row += 4){
				const int wa = gate_row[0], wb = gate_row[1], wo = gate_row[2];
				const int gtype = gate_row[3];
				if (gtype == AND_GATE){
					block outk = k1[num_in+counter_and];
					if(d_gamma[counter_and])
						outk = outk ^ Delta;
					wk[wo] = outk;

					const bool da = d_alpha[counter_and];
					const bool db = d_beta[counter_and];
					block tmp_check = vole_wire_y_key[counter_and] ^ outk;
					if (da && db)
						tmp_check = tmp_check ^ Delta;
					if(da)
						tmp_check = tmp_check ^ wk[wb];
					if(db)
						tmp_check = tmp_check ^ wk[wa];
					mt_B[counter_and] = tmp_check;

					counter_and++;
				}
				else if (gtype == XOR_GATE){
					wk[wo] = wk[wa] ^ wk[wb];
				}
				else{
					wk[wo] = wk[wa] ^ Delta;
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

			block * chi = new block[num_ands];
			PRG prg_chi;
			prg_chi.random_block(chi, num_ands);
			block com_mt_A[2];
			com_mt_A[0] = zero_block;
			com_mt_A[1] = zero_block;
			block com_mt_b = zero_block;

			for(int i = 0; i < num_ands; ++i){
				block tmp;
				gfmul(recv_mt_A[2*i], chi[i], &tmp);
				com_mt_A[0] = com_mt_A[0] ^ tmp;
				gfmul(recv_mt_A[2*i+1], chi[i], &tmp);
				com_mt_A[1] = com_mt_A[1] ^ tmp;
				gfmul(mt_B[i], chi[i], &tmp);
				com_mt_b = com_mt_b ^ tmp;
			}

			block tmp_check;
			gfmul(com_mt_A[1], Delta, &tmp_check);
			tmp_check = tmp_check ^ com_mt_A[0];
			if(!cmpBlock(&tmp_check, &com_mt_b, 1))
				cout << "Random combination of MT fails" << endl;
			
		}
	}
};