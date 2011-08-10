#ifndef MIMOMAC_H
#define MIMOMAC_H

#include <itpp/itcomm.h>

/**
 *	Mimo Multiple-Access Channel
 *	Each user has only one antenna
 */
class MimoMac {
public:
	MimoMac(int nuser, int nrx);
	virtual ~MimoMac();
	
	// channel
	Array<cvec> channel(const Array<cvec> &tx_signals);

protected:
	
	// Number of user
	int num_user;
	
	// Number of rx antenna of the receiver
	int num_rx_ant;

	cmat H;
	


};



#endif
