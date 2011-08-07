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

#include <sys/time.h>

#include "comm.h"
#include "channel.h"
#include "debug_log.h"


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

	int msg_len = 100;
	int sym_len = msg_len/2;


	/////////////////////////////////////////////////
	// Channel initialization
	/////////////////////////////////////////////////
	complex<double> *ch_gain = new complex<double>[sym_len];

	// AWGN
	AwgnChannel awgn_channel(sigma2);
	for(int i=0; i<num_user; i++)
		ch_gain[i] = 1.0;


	cout<<EsN0<<endl;

	return 0;
}
