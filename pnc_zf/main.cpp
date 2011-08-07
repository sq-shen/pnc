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

#include "../../comm/comm.h"
#include "../../comm/random.h"
#include "../../comm/modulator.h"
#include "../../comm/channel.h"
#include "../../comm/debug_log.h"
#include "CRelay.h"



using namespace std;

int main(int argc, char *argv[]) 
{
	srand(1);

	///////////////////////////
	// Open debuging log file
	///////////////////////////
	OpenDbg1("dbglog1.txt");
	OpenDbg2("dbglog2.txt");

	//////////////////////////////////////////////////
	// Simulation parameters
	//////////////////////////////////////////////////
	double Es = 1;
	double EsN0dB = 20;
	double EsN0 = pow(10, EsN0dB/10);
	double N0 = 1/EsN0;
	double sigma2 = N0/2;
	double sigma  = sqrt(sigma2);


	// block
	int block_num = 1;
	int msg_len = 10;
	int sym_len = msg_len/2;  // QPSK

	// sequences
	char *bk_msg_u1  = new char[msg_len];
	char *bk_msg_u2  = new char[msg_len];
	char *bk_msg_xor = new char[msg_len];

	complex<double> *bk_tx_sym_u1 = new complex<double>[sym_len];
	complex<double> *bk_tx_sym_u2 = new complex<double>[sym_len];

	//complex<double>
	complex<double> *bk_rx_sym = new complex<double>[sym_len];


	/////////////////////////////////////////////////
	// Channel initialization
	/////////////////////////////////////////////////
	complex<double> h1 = 1.0;
	complex<double> h2 = 1.0;


	/////////////////////////////////////////////////
	// Users' modulators
	/////////////////////////////////////////////////
	CQpsk qpsk(GrayNonNor);

	/////////////////////////////////////////////////
	// PNC Relay
	/////////////////////////////////////////////////
	dem_region_1d x_dem_rgn[3];
	dem_region_1d y_dem_rgn[3];
	int dec_label[3][3];

	x_dem_rgn[0].low_db = numeric_limits<double>::min();
	x_dem_rgn[0].up_bd  = -sqrt(2.0)/2;
	x_dem_rgn[0].idx    = 0;

	x_dem_rgn[1].low_db = -sqrt(2.0)/2;
	x_dem_rgn[1].up_bd  = sqrt(2.0)/2;
	x_dem_rgn[1].idx    = 1;

	x_dem_rgn[2].low_db = sqrt(2.0)/2;
	x_dem_rgn[2].up_bd  = numeric_limits<double>::max();
	x_dem_rgn[2].idx    = 2;

	y_dem_rgn[0].low_db = numeric_limits<double>::min();
	y_dem_rgn[0].up_bd  = -sqrt(2.0)/2;
	y_dem_rgn[0].idx    = 0;

	y_dem_rgn[1].low_db = -sqrt(2.0)/2;
	y_dem_rgn[1].up_bd  = sqrt(2.0)/2;
	y_dem_rgn[1].idx    = 1;

	y_dem_rgn[2].low_db = sqrt(2.0)/2;
	y_dem_rgn[2].up_bd  = numeric_limits<double>::max();
	y_dem_rgn[2].idx    = 2;

	dec_label[0][0] = 0;	// 00
	dec_label[0][2] = 0;	// 00
	dec_label[2][0] = 0;	// 00
	dec_label[2][2] = 0;	// 00

	dec_label[0][1] = 1;	// 01
	dec_label[2][1] = 1;	// 01

	dec_label[1][0] = 2;	// 10
	dec_label[1][2] = 2;	// 10

	dec_label[1][1] = 3;	// 11


	CRelay relay;
	relay.set_rx_sym_len(sym_len);
	relay.init_dem_region(3, x_dem_rgn, 3, y_dem_rgn);




	/////////////////////////////////////////////////
	// Simulation
	/////////////////////////////////////////////////
	int bk;

	for(bk=1; bk<=block_num; bk++) {
		
		// generate message bits
		for(int i=0; i<msg_len; i++) {
			bk_msg_u1[i] = RandBin();
			bk_msg_u2[i] = RandBin();
			bk_msg_xor[i] = (bk_msg_u1[i] ^ bk_msg_u2[i]) & 0x01;
		}

		// modulation
		qpsk.modulate_bits(msg_len, bk_msg_u1, bk_tx_sym_u1, 1.0);
		qpsk.modulate_bits(msg_len, bk_msg_u2, bk_tx_sym_u2, 1.0);

		// transmit through the channel
		for(int i=0; i<sym_len; i++) {
			bk_rx_sym[i] = h1 * bk_tx_sym_u1[i] + h2 * bk_tx_sym_u2[i] + sigma * randn_c();
			//bk_rx_sym[i] = h1 * bk_tx_sym_u1[i] + h2 * bk_tx_sym_u2[i];
		}


		WriteArrayDbg1("bk_msg_u1", msg_len, bk_msg_u1, int);
		WriteArrayDbg1("bk_msg_u2", msg_len, bk_msg_u2, int);
		WriteArrayDbg2("bk_msg_xor", msg_len, bk_msg_xor, int);
		WriteArrayDbg2("bk_rx_sym", sym_len, bk_rx_sym, complex<double>);

	}





	// free memory
	delete [] bk_msg_u1;
	delete [] bk_msg_u2;
	delete [] bk_tx_sym_u1;
	delete [] bk_tx_sym_u2;

	///////////////////////////
	// Close debuging log file
	CloseDbg1();
	CloseDbg2();

	return 0;
}
