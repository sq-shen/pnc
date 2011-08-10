#ifndef MIMOMAC_H
#define MIMOMAC_H

#include <itpp/itcomm.h>

/**
 *	Mimo Multiple-Access Channel
 */
class MimoMac {
public:
	MimoMac(int nuser, int ntx_ant, int rx_ant);
	virtual ~MimoMac();
	
protected:
	
	// Number of user
	int num_user;
	
	// Number of tx antenna of "each user"
	int num_tx_ant;
	
	// Number of rx antenna of the receiver
	int num_rx_ant;
	

};



#endif