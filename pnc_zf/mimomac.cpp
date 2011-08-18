#include "mimomac.h"

using namespace itpp;
using namespace std;

MimoMac::MimoMac(int nuser, int nrx) {

	num_user = nuser;
	num_rx_ant = nrx;
	
	// H: nrx x nuser
	H = randn_c(num_rx_ant, num_user);
}

MimoMac::MimoMac(int nuser, int nrx, double var_2d) {

	num_user = nuser;
	num_rx_ant = nrx;
	N0 = var_2d;

	// H: nrx x nuser
	H = randn_c(num_rx_ant, num_user);
}

MimoMac::~MimoMac() { 
}


// channel
Array<cvec> MimoMac::channel(const Array<cvec> &tx_signals) {
	
	Array<cvec> rx_signals(tx_signals.size());
	double sigma = sqrt(N0);
	for(int i=0; i<tx_signals.size(); i++) {
		cvec in = tx_signals(i);
		cvec out = H * in + sigma*randn_c(num_rx_ant);	// y = Hx + z
		rx_signals(i) = out;
	}

	return rx_signals;
}
	
// Generate new channel
void MimoMac::genH() {
	H = randn_c(num_rx_ant, num_user);
}

