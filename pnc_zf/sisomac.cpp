#include "sisomac.h"

using namespace itpp;
using namespace std;

SisoMac::SisoMac(int nuser) {

	num_user = nuser;
	
	// h: 1 x nuser
	h = randn_c(num_user);
}

SisoMac::SisoMac(int nuser, double var_2d) {

	num_user = nuser;
	N0 = var_2d;

	// h: 1 x nuser
	h = randn_c(num_user);
}

SisoMac::~SisoMac() {
}


// channel
Array<complex<double> > SisoMac::channel(const Array<cvec> &tx_signals) {
	
	Array<complex<double> > rx_signals(tx_signals.size());
	double sigma = sqrt(N0);
	for(int i=0; i<tx_signals.size(); i++) {
		cvec in = tx_signals(i);
		complex<double> out = h * in + sigma*randn_c();	// y = hx + z
		// complex<double> out = h * in;	// zero noise: y = hx
		rx_signals(i) = out;
	}

	return rx_signals;
}
	
// Generate new channel
void SisoMac::gen_h() {
	h = randn_c(num_user);
}

