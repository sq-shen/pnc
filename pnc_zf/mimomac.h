#ifndef MIMOMAC_H
#define MIMOMAC_H

#include <itpp/itcomm.h>

/**
 *	Mimo Multiple-Access Channel
 *	Each user has only one antenna
 */
class MimoMac {
public:
	MimoMac(int nuser, int nrx, double var_2d);
	virtual ~MimoMac();
	
	// channel
	itpp::Array<itpp::cvec> channel(const itpp::Array<itpp::cvec> &tx_signals);
	
	// Generate new channel
	void genH();
	
	itpp::cmat get_H();

protected:
	
	// Number of user
	int num_user;
	
	// Number of rx antenna of the receiver
	int num_rx_ant;
	
	// Noiser variance
	double N0;

	// mimo channel: num_rx_ant x num_user
	itpp::cmat H;

};

inline itpp::cmat MimoMac::get_H() {
	return H;
}


#endif
