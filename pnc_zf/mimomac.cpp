#include "mimomac.h"

using namespace itpp;
using namespace std

MimoMac::MimoMac(int nuser, int nrx) {

	num_user = nuser;
	num_rx_ant = nrx;
	
	H.set_size(nrx, nuser);
}

MimoMac::~MimoMac() { 
}



