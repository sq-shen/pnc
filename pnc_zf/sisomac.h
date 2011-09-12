#ifndef SISOMAC_H
#define SISOMAC_H

#include <itpp/itcomm.h>

/**
 *	SISO Multiple-Access Channel
 *	Each user has only one antenna
 */
class SisoMac {
public:
	SisoMac(int nuser);
	SisoMac(int nuser, double var_2d);
	virtual ~SisoMac();
	
	// channel
	itpp::Array<std::complex<double> > channel(const itpp::Array<itpp::cvec> &tx_signals);
	
	// Generate new channel
	void gen_h();
	
	void set_h(itpp::cvec &ch);
	
	void set_N0(double var_2d);

	itpp::cvec get_h();

protected:
	
	// Number of user
	int num_user;

	// Noiser variance
	double N0;

	// siso channel: num_rx_ant x num_user
	itpp::cvec h;

};

inline itpp::cvec SisoMac::get_h() {
	return h;
}

inline void SisoMac::set_N0(double var_2d) {
	N0 = var_2d;
}

inline void SisoMac::set_h(itpp::cvec &ch) {
	h = ch;
}


#endif
