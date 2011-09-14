/**
 * PNC - Change of basis
 *
 */

#include <itpp/itstat.h>
#include <itpp/itcomm.h>

#include "sisomac.h"
#include "ant1relay.h"

using namespace itpp;
using namespace std;

Array<cvec> to_mimo_input(Array<cvec> &user_input) {

	int sym_len = user_input(0).size();
	Array<cvec> mimo_in(sym_len);
	for(int i=0; i<sym_len; i++) {
		cvec in(user_input.size());
		for(int k=0; k<user_input.size(); k++) {
			in(k) = user_input(k).get(i);
		}
		mimo_in(i) = in;
	}

	return mimo_in;
}

int qam_sym_err(bvec bv_src, ivec iv_dem) {
	int err = 0;
	for(int i=0, k=0; i<iv_dem.size(); i++, k+=2) {
		int src = (bv_src(k).value()<<1) + bv_src(k+1).value();
		if(src!=iv_dem(i))
			err++;
	}
	return err;
}

int bpsk_sym_err(bvec bv_src, ivec iv_dem) {
	int err = 0;
	for(int i=0; i<iv_dem.size(); i++) {
		if(bv_src(i)!=iv_dem(i))
			err++;
	}
	return err;
}


int main(int argc, char *argv[])
{
	RNG_randomize();
	Real_Timer tt;
	
	bool is_bpsk = true;

	char buf[1024];

	// log file
	ofstream log;
	log.open("log_pnc_cord.txt", ios::app|ios::out);
	log<<"==================================================================================================="<<endl;

	//////////////////////////////////////////////////
	// Simulation parameters
	//////////////////////////////////////////////////
	int num_user = 2;
	//int num_rx_ant = 1;
	
	
	/////////////////////////////////////////////////
	// Users' modulators
	/////////////////////////////////////////////////
	QAM qam(4);		// 4-QAM
	BPSK_c bpsk;	// BPSK
	Modulator<complex<double> > *mod = NULL;
	if(is_bpsk) {
		mod = &bpsk;
	} else {
		mod = &qam;
	}
	cvec syms = mod->get_symbols();
	cout<<"constellation: "<<syms<<endl;

	int block_num = 10000;
	int msg_len = 1000;
	int bits_per_symbol = mod->bits_per_symbol();
	int sym_len = msg_len/bits_per_symbol;
	
	double Es = 1;	
	// vec EsN0dB  = linspace(30,30,1);
	vec EsN0dB  = linspace(10,30,21);
	vec EsN0    = inv_dB(EsN0dB);
	vec N0      = Es * pow(EsN0, -1.0);
	vec sqrt_N0 = sqrt(N0);
	vec sigma2  = N0/2;
	vec sigma   = sqrt(sigma);


	bvec bv_msg_u1, bv_msg_u2;
	bvec bv_msg_xor;

	cvec cv_tx_sym_u1, cv_tx_sym_u2;
	cvec cv_rx_sym;

	/////////////////////////////////////////////////
	// Channel initialization
	/////////////////////////////////////////////////
	SisoMac sisomac(num_user);

	/////////////////////////////////////////////////
	// PNC Relay
	/////////////////////////////////////////////////
	Ant1Relay relay;


	/////////////////////////////////////////////////
	// Simulation
	/////////////////////////////////////////////////
	cout<<"Start simulation"<<endl;
	log<<"#EsN0dB       #err         SER"<<endl;
	for(int i=0; i<EsN0dB.size(); i++) {

		int tot_sym = 0, err = 0;
		for(int bk=1; bk<=block_num; bk++) {

			//======================================
			// generate message bits
			//======================================
			bv_msg_u1 = randb(msg_len);
			bv_msg_u2 = randb(msg_len);
			bv_msg_xor = bv_msg_u1 + bv_msg_u2;

			//======================================
			// modulation
			//======================================
			cvec cv_txsig_u1 = mod->modulate_bits(bv_msg_u1);
			cvec cv_txsig_u2 = mod->modulate_bits(bv_msg_u2);

			//======================================
			// prepare proper format for mimo input
			//======================================
			Array<cvec> user_input(2);
			user_input(0) = cv_txsig_u1;
			user_input(1) = cv_txsig_u2;
			Array<cvec> mimo_input = to_mimo_input(user_input);

			//======================================
			// channel
			//======================================
			sisomac.gen_h();
			sisomac.set_N0(N0(i));
			Array<complex<double> > miso_output = sisomac.channel(mimo_input);

			// relay
			cvec h = sisomac.get_h();
			relay.set_channel(h);


			//======================================
			// MIMO-PNC Demapping
			//======================================
			ivec dem_sym;
			// dem_sym = relay.bpsk_chcord_demapping(miso_output);
			// dem_sym = relay.bpsk_nc_ml_demapping(miso_output, syms);
			dem_sym = relay.bpsk_mrc_pnc_demapping(miso_output);

			// dem_sym = relay.bpsk_vecproj_pnc_demapping(miso_output);

			// cout<<"bv_msg_xor="<<bv_msg_xor<<endl;
			// cout<<"dem_sym="<<dem_sym<<endl;

			if(is_bpsk)
				err += bpsk_sym_err(bv_msg_xor, dem_sym);
			else
				err += qam_sym_err(bv_msg_xor, dem_sym);

			tot_sym += sym_len;
			printf("EsN0dB=%4.1f, #bk=%6i, #err=%6i, SER=%1.3e\r",
					EsN0dB(i), bk, err, ((double)err)/tot_sym);

		}
		printf("\n");

		// log
		sprintf(buf, "  %4.1f      %6i    %1.3e", EsN0dB(i), err, ((double)err)/tot_sym);
		log<<buf<<endl;
	}

	log<<"===================================================================================================\n\n"<<endl;
	log.close();


	return 0;
}
