/**
 *	Simulation of PNC-ZF using IT++
 *	Two antennas in the relay
 *
 */

#include <itpp/itcomm.h>


#include "../../comm/comm.h"
// #include "../../comm/random.h"
// #include "../../comm/modulator.h"
// #include "../../comm/channel.h"
#include "../../comm/debug_log.h"
#include "CRelay.h"


using namespace itpp;
using namespace std;


int main(int argc, char *argv[]) 
{
	//srand(1);
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


	// block
	int block_num = 1;
	int msg_len = 10;
	int sym_len = msg_len/2;  // QPSK

	bvec bv_msg_u1, bv_msg_u2;
	bvec bv_msg_xor;
	
	cvec cv_tx_sym_u1, cv_tx_sym_u2;
	cvec cv_rx_sym;


	/////////////////////////////////////////////////
	// Channel initialization
	/////////////////////////////////////////////////
	//AWGN_Channel awgn_channel;     //The AWGN channel class
	complex<double> H[2][2] = {
		{1.6, 0.3},
		{0.5, 1.2}
	};


	/////////////////////////////////////////////////
	// Users' modulators
	/////////////////////////////////////////////////
	QPSK qpsk;                     //The QPSK modulator class


	/////////////////////////////////////////////////
	// PNC Relay
	// Demapping region (x1 + x2)
	/////////////////////////////////////////////////





	/////////////////////////////////////////////////
	// Simulation
	/////////////////////////////////////////////////
	int bk;
	int tot_sym = 0;
	int err = 0;
	double ser;

	for(bk=1; bk<=block_num; bk++) {
		
		// generate message bits
		bv_msg_u1 = randb(msg_len);
		bv_msg_u2 = randb(msg_len);
		bv_msg_xor = bv_msg_u1 + bv_msg_u2;
		
		// modulation
		
		
		cout<<bv_msg_u1<<endl;
		cout<<bv_msg_u2<<endl;
		cout<<bv_msg_xor<<endl;
		
	}
	printf("\nEnd\n");



	///////////////////////////
	// Close debuging log file
	CloseDbg1();
	CloseDbg2();

	return 0;
}
