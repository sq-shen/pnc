/**
 *	Simulation of PNC-ZF
 *
 */

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <cmath>
#include <sstream>
#include <limits>

#include <sys/time.h>

#include "comm.h"
#include "random.h"
#include "modulator.h"
#include "channel.h"
#include "debug_log.h"
#include "CRelay.h"



using namespace std;

int main(int argc, char *argv[]) 
{
	//////////////////////////////////////////////////
	// Simulation parameters
	//////////////////////////////////////////////////
	double Es = 1;
	double EsN0dB = 20;
	double EsN0 = pow(10, EsN0dB/10);
	double N0 = 1/EsN0;
	double sigma2 = N0/2;


	// block
	int block_num = 1;
	int msg_len = 100;
	int sym_len = msg_len/2;  // QPSK

	// sequences
	char *bk_msg_u1 = new char[msg_len];
	char *bk_msg_u2 = new char[msg_len];

	complex<double> *bk_tx_sym_u1 = new complex<double>[sym_len];
	complex<double> *bk_tx_sym_u2 = new complex<double>[sym_len];

	//complex<double>
	complex<double> *bk_rx_sym = new complex<double>[sym_len];


	/////////////////////////////////////////////////
	// Channel initialization
	/////////////////////////////////////////////////

	AwgnChannel awgn_channel(N0); // AWGN


	/////////////////////////////////////////////////
	// Users' modulators
	/////////////////////////////////////////////////
	CQpsk qpsk(GrayNonNor);

	/////////////////////////////////////////////////
	// PNC Relay
	/////////////////////////////////////////////////
	CRelay relay;
	relay.set_rx_sym_len(sym_len);

	numeric_limits<double>::max();





	/////////////////////////////////////////////////
	// Simulation
	/////////////////////////////////////////////////
	int bk;
	for(bk=1; bk<=block_num; bk++) {
		
		// generate message bits
		for(int i=0; i<msg_len; i++) {
			bk_msg_u1[i] = RandBin();
			bk_msg_u2[i] = RandBin();
		}

		// modulation


	}





	// free memory
	delete [] bk_msg_u1;
	delete [] bk_msg_u2;
	delete [] bk_tx_sym_u1;
	delete [] bk_tx_sym_u2;

	return 0;
}
