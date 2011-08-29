/**
 *	Simulation of PNC-MMSE using IT++
 *	Two antennas in the relay
 *
 */
#include <sstream>
 
#include <itpp/itstat.h>
#include <itpp/itcomm.h>

#include "mimomac.h"
#include "zftwrelay.h"

using namespace itpp;
using namespace std;



int main(int argc, char *argv[]) 
{
	// Fixed channel matrix from file
	cmat read_H;
	bool is_fixed_H = false;
	if(argc>=2) {
		it_file ff;
		ff.open(argv[1]);
		ff>>Name("H")>>read_H;
		ff.close();	
		cout<<"Read H: "<<read_H<<endl;
		is_fixed_H = true;
	} else {
		cout<<"draw_const [H file]"<<endl;
		return 0;
	}

	Real_Timer tt;

	char buf[1024];

	//////////////////////////////////////////////////
	// Simulation parameters
	//////////////////////////////////////////////////
	int num_user = 2;
	int num_rx_ant = 2;

	int block_num = 1;
	int msg_len = 10;
	int sym_len = msg_len/2;  // QPSK

	double EsN0dB  = 15;
	double Es 	   = 1;
	double EsN0    = inv_dB(EsN0dB);
	double N0      = Es * pow(EsN0, -1.0);
	double sqrt_N0 = sqrt(N0);
	double sigma2  = N0/2;
	double sigma   = sqrt(sigma2);
	
	cvec cv_tx_sym_u1, cv_tx_sym_u2;
	cvec cv_rx_sym;
	

	/////////////////////////////////////////////////
	// Users' modulators
	/////////////////////////////////////////////////
	QAM qam(4);                     //The 4-QAM modulator class
	cvec qam_syms = qam.get_symbols();
	//cout<<"4-QAM constellation: "<<syms<<endl;

	/////////////////////////////////////////////////
	// Channel initialization
	/////////////////////////////////////////////////
	MimoMac mimomac(num_user, num_rx_ant);
	mimomac.set_H(read_H);


	/////////////////////////////////////////////////
	// PNC Relay
	// Demapping region (x1 + x2)
	/////////////////////////////////////////////////
	ZfTwRelay relay;
	
	cvec w_conj, a;
	cmat H, Ry;
	
	H = mimomac.get_H();
	relay.set_H(H);
	Ry = relay.calc_Ry(N0);
	a = relay.pnc_mmse_a(N0);
	
	w_conj = relay.pnc_mmse_detector(a, Ry);
	//w_conj = ones_c(2);
	
	relay.init_sp_const(w_conj, qam_syms, qam_syms);
	relay.show_sp_constellation();
	
	
	///////////////////////////////////////////////////
	// Save to file
	///////////////////////////////////////////////////
	Array<cvec> pts(4);
	for(int i=0; i<pts.size(); i++)
		pts(i).set_size(4);
	ivec idxs = zeros_i(4);
	
	Array<sp_pt> sp_const = relay.get_sp_constellation();
	for(int i=0; i<sp_const.size(); i++) {
		int label = sp_const(i).label;
		pts(label)(idxs(label)++) = sp_const(i).pt;
	}
	
	// log file
	ofstream file;
	file.open("sp_constellation.txt", ios::out);
	for(int i=0; i<4; i++) {
		stringstream stream;
		for(int j=0; j<pts.size(); j++) {
			complex<double> pt = round_to_zero(pts(j)(i));
			stream<<pt.real()<<"\t"<<pt.imag()<<"\t";
		}
		file<<stream.str()<<endl;
	}
	file.flush();
	file.close();
	
	return 0;
}
