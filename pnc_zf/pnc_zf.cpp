/**
 *	Simulation of PNC-ZF using IT++
 *	Two antennas in the relay
 *
 */

#include <itpp/itcomm.h>

#include "../../comm/comm.h"
#include "../../comm/debug_log.h"
#include "mimomac.h"
#include "zftwrelay.h"

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

int sym_err(bvec bv_src, ivec iv_dem) {
	int err = 0;
	for(int i=0, k=0; i<iv_dem.size(); i++, k+=2) {
		int src = (bv_src(k).value()<<1) + bv_src(k+1).value();
		if(src!=iv_dem(i))
			err++;
	}
	return err;
}


int main(int argc, char *argv[]) 
{
	RNG_randomize();
	
	Real_Timer tt;

	///////////////////////////
	// Open debuging log file
	///////////////////////////
	OpenDbg1("dbglog1.txt");
	OpenDbg2("dbglog2.txt");

	//////////////////////////////////////////////////
	// Simulation parameters
	//////////////////////////////////////////////////
	double Es = 1;
	double EsN0dB = 15;
	double EsN0 = pow(10, EsN0dB/10);
	double N0 = 1/EsN0;
	double sqrt_N0 = sqrt(N0);
	double sigma2 = N0/2;
	double sigma  = sqrt(sigma2);

	int num_user = 2;
	int num_rx_ant = 2;


	// block
	int block_num = 100;
	int msg_len = 10000;
	int sym_len = msg_len/2;  // QPSK

	bvec bv_msg_u1, bv_msg_u2;
	bvec bv_msg_xor;
	
	cvec cv_tx_sym_u1, cv_tx_sym_u2;
	cvec cv_rx_sym;
	
	/////////////////////////////////////////////////
	// Linear combination coefficient
	/////////////////////////////////////////////////
	vec a = "1 1";
	cout<<a(0)<<","<<a(1)<<endl;
	
	/////////////////////////////////////////////////
	// Users' modulators
	/////////////////////////////////////////////////
	QAM qam(4);                     //The 4-QAM modulator class
	cvec syms = qam.get_symbols();
	cout<<"4-QAM constellation: "<<syms<<endl;

	/////////////////////////////////////////////////
	// Channel initialization
	/////////////////////////////////////////////////
	MimoMac mimomac(num_user, num_rx_ant, N0);
	mimomac.genH();
	cmat H = mimomac.get_H();
	cout<<"H="<<H<<endl;

	/////////////////////////////////////////////////
	// PNC Relay
	// Demapping region (x1 + x2)
	/////////////////////////////////////////////////
	ZfTwRelay relay;
	relay.set_H(H);
	relay.init_dem_region(a, syms, syms);
	relay.show_sp_constellation();
	relay.show_dem_regions();


	/////////////////////////////////////////////////
	// Simulation
	/////////////////////////////////////////////////
	int bk;
	int tot_sym = 0;
	int err = 0;
	double ser;

	for(bk=1; bk<=block_num; bk++) {
		
		//======================================
		// generate message bits
		//======================================
		bv_msg_u1 = randb(msg_len);
		bv_msg_u2 = randb(msg_len);
		bv_msg_xor = bv_msg_u1 + bv_msg_u2;
		
		//======================================
		// modulation
		//======================================
		cvec cv_txsig_u1 = qam.modulate_bits(bv_msg_u1);
		cvec cv_txsig_u2 = qam.modulate_bits(bv_msg_u2);
		
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
		Array<cvec> mimo_output = mimomac.channel(mimo_input);
		
		
		//======================================
		// PNC Demapping
		//======================================
		ivec dem_sym = relay.pnc_demapping(mimo_output);
		
		// cout<<"bv_msg_xor="<<bv_msg_xor<<endl;
		// cout<<"dem_sym="<<dem_sym<<endl;
		
		err += sym_err(bv_msg_xor, dem_sym);
		
		tot_sym += sym_len;
		printf("#bk=%6i, #err=%6i, ser=%1.3e\r", bk, err, ((double)err)/tot_sym);

		// cout<<bv_msg_u1<<endl;
		// cout<<cv_txsig_u1<<endl;
		
		// cout<<bv_msg_u2<<endl;
		// cout<<cv_txsig_u2<<endl;
		
		// cout<<user_input<<endl;
		// cout<<mimo_input<<endl;

		// cout<<mimo_output<<endl;
		
	}
	printf("\nEnd\n");



	///////////////////////////
	// Close debuging log file
	CloseDbg1();
	CloseDbg2();

	return 0;
}
