#pragma once

#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include "emp-agmpc/netmp.h"
#include "emp-agmpc/helper.h"
#include "emp-tool/utils/block.h"
#include "emp-tool/circuits/circuit_file.h"
#include "emp-tool/utils/f2k.h"
#include <vector>
#include <cstring>

using namespace std;
using namespace emp;

// mv-qsvole over binary field via IKNP OTE
template<int nP>
class MVSVOLEBOOL_IKNP { public:
	IKNP<NetIO> *abit[nP+1];
	NetIOMP<nP> *io;
	ThreadPool * pool;
	int party;
	PRG prg;
	block Delta;
	Hash hash;
	int ssp;

	MVSVOLEBOOL_IKNP(NetIOMP<nP>* io, ThreadPool * pool, int party, int ssp = 40) {
		this->io = io;
		this->pool = pool;
		this->party = party;
    	this->ssp = ssp;

		if (party == 1){
			for(int i = 2; i <= nP; ++i){
				// acts as COT receiver
				abit[i] = new IKNP<NetIO>(io->get(i, false));
			}
		}
		else{
			// act as COT sender
			abit[1] = new IKNP<NetIO>(io->get(1, false));
		}

       

	}
	~MVSVOLEBOOL_IKNP() {
		if (party == 1){
			for (int i=2; i <= nP; ++i){
				delete abit[i];
			}
		}
		else{
			delete abit[1];
		}
	}
	
	void computeVOLE_semi(block * output[nP+1], bool* data, int length) {

		if (party == 1){
            // acts as COT receiver

            prg.random_bool(data, length);

			vector<future<void>> res;
			for (int i = 2; i<=nP; ++i){
				int party2 = i;
				res.push_back(pool->enqueue([this, output, data, length, party2]() {
					abit[party2]->recv_cot(output[party2], data, length);
					io->flush(party2);
				}));
			}
			joinNclean(res);
		}
		else{
            // acts as COT sender
			abit[1]->send_cot(output[1], length);
            Delta = abit[1]->Delta;
			io->flush(1);

		}
		
	}

    void computeVOLE_mal(block * output[nP+1], bool* data, int length){

        computeVOLE_semi(output, data, length);

        block coin = sampleRandom(io, &prg, pool, party); 
		bool * s = new bool[length*ssp];
		PRG prg2(&coin);
		prg2.random_bool(s, length*ssp);

        bool * u = new bool[ssp];
        
        if (party == 1){

            // compute u
            for(int i = 0; i < ssp; ++i){
                u[i] = 0;
                for (int j = 0; j < length; ++j){
                    int index = i*length + j;
                    if(s[index]){
                        u[i] = u[i] ^ data[j];
                    }
                }
            }

            block * tmpM[nP+1];
            for(int i = 2; i <=nP;++i){
                tmpM[i] = new block[ssp];
            }

            // compute tmpM
            for (int i = 2; i <= nP; ++i){
                for (int j = 0; j < ssp; ++j){
                    tmpM[i][j] = zero_block;
                    for(int k = 0; k < length; ++k){
                        int index = j*length + k;
                        if(s[index]){
                            tmpM[i][j] = tmpM[i][j] ^ output[i][k];
                        }

                    }
                }
            }

            // send to the verifiers
            vector<future<void>> res;
			for (int i = 2; i<=nP; ++i){
				int party2 = i;
				res.push_back(pool->enqueue([this, u, tmpM, party2]() {
					io->send_data(party2, u, ssp*sizeof(bool));
					io->send_data(party2, tmpM[party2], ssp*sizeof(block));
				}));
			}
			joinNclean(res);

            for(int i = 2; i <=nP;++i){
                delete[] tmpM[i];
            }
        }
        else{
            block * tmpK = new block[ssp];
            block * tmpM = new block[ssp];

            // compute tmpK for party 1
			for(int i = 0; i <ssp; ++i){
				tmpK[i] = zero_block;
				for(int j = 0; j < length; ++j){
					int index = i*length + j;
					if (s[index]){
						tmpK[i] = tmpK[i] ^ output[1][j];
					}	
				}
			}

            io->recv_data(1, u, ssp*sizeof(bool));
			io->recv_data(1, tmpM, ssp*sizeof(block));

            // check the MAC
            bool flag = true;
            block tmp;
            for(int i = 0; i < ssp; ++i){
                tmp = tmpM[i];
                if(u[i])
                    tmp = tmp ^ Delta;
                if(!cmpBlock(&tmp, &tmpK[i], 1))
                    flag = false;
            }
            if (!flag)
                cout << "Consistency checks for mv-qsvole fail" << endl;
        }

    }

    /* debug function*/
    void debug_check(block * output[nP+1], bool* data, int length){
        if(party == 1){
            for(int i = 2; i <= nP; ++i){
                io->send_data(i, output[i], length*sizeof(block));
                io->send_data(i, data, length*sizeof(bool));
            }
        }
        else{
            block * recv_mac = new block[length];
            io->recv_data(1, recv_mac, length*sizeof(block));
            io->recv_data(1, data, length*sizeof(bool));

            for(int i = 0; i < length; ++i){
                block tmp = recv_mac[i];
                if(data[i])
                    tmp = tmp ^ Delta;
                if(!cmpBlock(&tmp, &output[1][i], 1))
                    cout << "sVOLE fails!" << endl;
            }
        }
    }
	
};

// --- IKNP bool SIF: preprocessing after computeVOLE_mal (bundled ctx → arity 1) ---

template<int nP>
struct BoolIknpPreprocessProverExtCtx {
	int num_in = 0;
	int num_ands = 0;
	block** mac = nullptr;
	bool* value = nullptr;
	GaloisFieldPacking* pack = nullptr;
	block** value_ext = nullptr;
	block** mac_ext = nullptr;

	BoolIknpPreprocessProverExtCtx() = default;
	BoolIknpPreprocessProverExtCtx(int num_in_, int num_ands_, block** mac_, bool* value_,
			GaloisFieldPacking& pack_, block** value_ext_, block** mac_ext_)
		: num_in(num_in_), num_ands(num_ands_), mac(mac_), value(value_), pack(&pack_),
		  value_ext(value_ext_), mac_ext(mac_ext_) {}
};

template<int nP>
struct BoolIknpPreprocessProverMacCtx {
	BristolFashion* cf = nullptr;
	int num_in = 0;
	int num_ands = 0;
	bool* value = nullptr;
	block** mac = nullptr;
	bool* vole_wire_value = nullptr;
	block** vole_wire_mac = nullptr;
	bool* vole_wire_y = nullptr;
	bool* mask_y = nullptr;
	block** vole_wire_y_mac = nullptr;
	block** y_A = nullptr;

	BoolIknpPreprocessProverMacCtx() = default;
	BoolIknpPreprocessProverMacCtx(BristolFashion* cf_, int num_in_, int num_ands_, bool* value_, block** mac_,
			bool* vole_wire_value_, block** vole_wire_mac_, bool* vole_wire_y_, bool* mask_y_,
			block** vole_wire_y_mac_, block** y_A_)
		: cf(cf_), num_in(num_in_), num_ands(num_ands_), value(value_), mac(mac_),
		  vole_wire_value(vole_wire_value_), vole_wire_mac(vole_wire_mac_), vole_wire_y(vole_wire_y_),
		  mask_y(mask_y_), vole_wire_y_mac(vole_wire_y_mac_), y_A(y_A_) {}
};

template<int nP>
struct BoolIknpPreprocessProverQuadCtx {
	NetIOMP<nP>* io = nullptr;
	ThreadPool* pool = nullptr;
	int num_ands = 0;
	bool* mask_y = nullptr;
	block** y_A = nullptr;
	block** value_ext = nullptr;
	block** mac_ext = nullptr;
	block* seed = nullptr;

	BoolIknpPreprocessProverQuadCtx() = default;
	BoolIknpPreprocessProverQuadCtx(NetIOMP<nP>* io_, ThreadPool* pool_, int num_ands_, bool* mask_y_,
			block** y_A_, block** value_ext_, block** mac_ext_, block* seed_)
		: io(io_), pool(pool_), num_ands(num_ands_), mask_y(mask_y_), y_A(y_A_), value_ext(value_ext_),
		  mac_ext(mac_ext_), seed(seed_) {}
};

template<int nP>
struct BoolIknpPreprocessVerifierCtx {
	NetIOMP<nP>* io = nullptr;
	int party = 0;
	int num_in = 0;
	int num_ands = 0;
	block Delta = zero_block;
	BristolFashion* cf = nullptr;
	block* key1 = nullptr;
	GaloisFieldPacking* pack = nullptr;
	block** key_ext = nullptr;
	bool* mask_y = nullptr;
	block* vole_wire_key = nullptr;
	block* vole_wire_y_key = nullptr;
	block* y_B = nullptr;
	block* seed = nullptr;

	BoolIknpPreprocessVerifierCtx() = default;
	BoolIknpPreprocessVerifierCtx(NetIOMP<nP>* io_, int party_, int num_in_, int num_ands_, block Delta_,
			BristolFashion* cf_, block* key1_, GaloisFieldPacking& pack_, block** key_ext_,
			bool* mask_y_, block* vole_wire_key_, block* vole_wire_y_key_, block* y_B_, block* seed_)
		: io(io_), party(party_), num_in(num_in_), num_ands(num_ands_), Delta(Delta_), cf(cf_), key1(key1_),
		  pack(&pack_), key_ext(key_ext_), mask_y(mask_y_), vole_wire_key(vole_wire_key_),
		  vole_wire_y_key(vole_wire_y_key_), y_B(y_B_), seed(seed_) {}
};

template<int nP>
void bool_iknp_preprocess_prover_extension_vole(const BoolIknpPreprocessProverExtCtx<nP>& c) {
	block* ope_data = new block[128];
	bool* ope_bool = new bool[128];
	for (int i = 2; i <= nP; ++i) {
		int off_index = c.num_in + c.num_ands + (i - 2) * 128;
		memcpy(ope_data, c.mac[i] + off_index, 128 * sizeof(block));
		memcpy(ope_bool, c.value + off_index, 128 * sizeof(bool));
		uint64_t ch[2];
		for (int j = 0; j < 2; ++j) {
			if (ope_bool[64 * j + 63])
				ch[j] = 1;
			else
				ch[j] = 0;
			for (int k = 62; k >= 0; k--) {
				ch[j] <<= 1;
				if (ope_bool[64 * j + k])
					ch[j]++;
			}
		}
		c.value_ext[i][0] = makeBlock(ch[1], ch[0]);
		c.pack->packing(c.mac_ext[i], ope_data);
	}
	delete[] ope_data;
	delete[] ope_bool;
}

template<int nP>
void bool_iknp_preprocess_prover_mac_on_circuit(const BoolIknpPreprocessProverMacCtx<nP>& c) {
	for (int i = 0; i < c.num_in; ++i) {
		c.vole_wire_value[i] = c.value[i];
		for (int j = 2; j <= nP; ++j)
			c.vole_wire_mac[j][i] = c.mac[j][i];
	}
	int and_cnt = 0;
	for (int i = 0; i < c.cf->num_gate; ++i) {
		if (c.cf->gates[4 * i + 3] == AND_GATE) {
			c.vole_wire_value[c.cf->gates[4 * i + 2]] = c.value[c.num_in + and_cnt];
			bool u_alpha = c.vole_wire_value[c.cf->gates[4 * i]];
			bool u_beta = c.vole_wire_value[c.cf->gates[4 * i + 1]];
			bool u_gamma = u_alpha & u_beta;
			c.vole_wire_y[and_cnt] = u_gamma;
			c.mask_y[and_cnt] = u_gamma ^ c.value[c.num_in + c.num_ands + and_cnt];
			for (int j = 2; j <= nP; ++j) {
				c.vole_wire_mac[j][c.cf->gates[4 * i + 2]] = c.mac[j][c.num_in + and_cnt];
				c.vole_wire_y_mac[j][and_cnt] = c.mac[j][c.num_in + c.num_ands + and_cnt];
				gfmul(c.vole_wire_mac[j][c.cf->gates[4 * i]], c.vole_wire_mac[j][c.cf->gates[4 * i + 1]],
					&c.y_A[j][2 * and_cnt]);
				c.y_A[j][2 * and_cnt + 1] = zero_block;
				if (u_alpha)
					c.y_A[j][2 * and_cnt + 1] =
						c.y_A[j][2 * and_cnt + 1] ^ c.vole_wire_mac[j][c.cf->gates[4 * i + 1]];
				if (u_beta)
					c.y_A[j][2 * and_cnt + 1] =
						c.y_A[j][2 * and_cnt + 1] ^ c.vole_wire_mac[j][c.cf->gates[4 * i]];
				c.y_A[j][2 * and_cnt + 1] = c.y_A[j][2 * and_cnt + 1] ^ c.vole_wire_y_mac[j][and_cnt];
			}
			and_cnt++;
		} else if (c.cf->gates[4 * i + 3] == XOR_GATE) {
			c.vole_wire_value[c.cf->gates[4 * i + 2]] =
				c.vole_wire_value[c.cf->gates[4 * i]] ^ c.vole_wire_value[c.cf->gates[4 * i + 1]];
			for (int j = 2; j <= nP; ++j) {
				c.vole_wire_mac[j][c.cf->gates[4 * i + 2]] =
					c.vole_wire_mac[j][c.cf->gates[4 * i]] ^ c.vole_wire_mac[j][c.cf->gates[4 * i + 1]];
			}
		} else {
			c.vole_wire_value[c.cf->gates[4 * i + 2]] = c.vole_wire_value[c.cf->gates[4 * i]] ^ true;
			for (int j = 2; j <= nP; ++j)
				c.vole_wire_mac[j][c.cf->gates[4 * i + 2]] = c.vole_wire_mac[j][c.cf->gates[4 * i]];
		}
	}
}

template<int nP>
void bool_iknp_preprocess_prover_quadratic_fold_and_send(const BoolIknpPreprocessProverQuadCtx<nP>& c) {
	*c.seed = Hash::hash_for_block(c.mask_y, c.num_ands * sizeof(bool));
	vector<block> chi(c.num_ands);
	PRG prg_chi;
	prg_chi.reseed(c.seed);
	prg_chi.random_block(chi.data(), c.num_ands);
	for (int i = 2; i <= nP; ++i) {
		block tmp_check;
		block com_y_A[2];
		com_y_A[0] = zero_block;
		com_y_A[1] = zero_block;
		for (int j = 0; j < c.num_ands; ++j) {
			gfmul(chi[j], c.y_A[i][2 * j], &tmp_check);
			com_y_A[0] = com_y_A[0] ^ tmp_check;
			gfmul(chi[j], c.y_A[i][2 * j + 1], &tmp_check);
			com_y_A[1] = com_y_A[1] ^ tmp_check;
		}
		c.value_ext[i][0] = c.value_ext[i][0] ^ com_y_A[1];
		c.mac_ext[i][0] = c.mac_ext[i][0] ^ com_y_A[0];
	}
	vector<future<void>> res;
	for (int i = 2; i <= nP; ++i) {
		int party2 = i;
		NetIOMP<nP>* io_ptr = c.io;
		bool* my = c.mask_y;
		int na = c.num_ands;
		block** ve = c.value_ext;
		block** me = c.mac_ext;
		res.push_back(c.pool->enqueue([io_ptr, my, na, ve, me, party2]() {
			io_ptr->send_data(party2, my, na * sizeof(bool));
			io_ptr->send_data(party2, ve[party2], sizeof(block));
			io_ptr->send_data(party2, me[party2], sizeof(block));
			io_ptr->flush(party2);
		}));
	}
	joinNclean(res);
}

template<int nP>
void bool_iknp_preprocess_verifier_extension_recv_and_check(const BoolIknpPreprocessVerifierCtx<nP>& c) {
	block* ope_data = new block[128];
	int off_index = c.num_in + c.num_ands + (c.party - 2) * 128;
	memcpy(ope_data, c.key1 + off_index, 128 * sizeof(block));
	c.pack->packing(c.key_ext[1], ope_data);
	delete[] ope_data;

	block recv_value_ext, recv_mac_ext;
	c.io->recv_data(1, c.mask_y, c.num_ands * sizeof(bool));
	c.io->recv_data(1, &recv_value_ext, sizeof(block));
	c.io->recv_data(1, &recv_mac_ext, sizeof(block));

	for (int i = 0; i < c.num_in; ++i)
		c.vole_wire_key[i] = c.key1[i];
	int and_cnt = 0;
	for (int i = 0; i < c.cf->num_gate; ++i) {
		if (c.cf->gates[4 * i + 3] == AND_GATE) {
			c.vole_wire_key[c.cf->gates[4 * i + 2]] = c.key1[c.num_in + and_cnt];
			c.vole_wire_y_key[and_cnt] = c.key1[c.num_in + c.num_ands + and_cnt];
			if (c.mask_y[and_cnt])
				c.vole_wire_y_key[and_cnt] = c.vole_wire_y_key[and_cnt] ^ c.Delta;
			block tmp_check;
			gfmul(c.vole_wire_key[c.cf->gates[4 * i]], c.vole_wire_key[c.cf->gates[4 * i + 1]], &tmp_check);
			c.y_B[and_cnt] = tmp_check;
			gfmul(c.vole_wire_y_key[and_cnt], c.Delta, &tmp_check);
			c.y_B[and_cnt] = c.y_B[and_cnt] ^ tmp_check;
			and_cnt++;
		} else if (c.cf->gates[4 * i + 3] == XOR_GATE) {
			c.vole_wire_key[c.cf->gates[4 * i + 2]] =
				c.vole_wire_key[c.cf->gates[4 * i]] ^ c.vole_wire_key[c.cf->gates[4 * i + 1]];
		} else {
			c.vole_wire_key[c.cf->gates[4 * i + 2]] = c.vole_wire_key[c.cf->gates[4 * i]] ^ c.Delta;
		}
	}

	*c.seed = Hash::hash_for_block(c.mask_y, c.num_ands * sizeof(bool));
	vector<block> chi(c.num_ands);
	PRG prg_chi;
	prg_chi.reseed(c.seed);
	prg_chi.random_block(chi.data(), c.num_ands);
	block tmp_check;
	block com_y_B = zero_block;
	for (int i = 0; i < c.num_ands; ++i) {
		gfmul(chi[i], c.y_B[i], &tmp_check);
		com_y_B = com_y_B ^ tmp_check;
	}
	c.key_ext[1][0] = c.key_ext[1][0] ^ com_y_B;
	gfmul(recv_value_ext, c.Delta, &tmp_check);
	tmp_check = tmp_check ^ recv_mac_ext;
	if (!cmpBlock(&tmp_check, &c.key_ext[1][0], 1))
		cout << "party " << c.party << " : quadratic multiplication checks fail" << endl;
}
