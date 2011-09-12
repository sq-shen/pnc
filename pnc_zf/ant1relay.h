/*
 * ant1relay.h
 *
 *  Created on: Aug 13, 2011
 *      Author: ivan
 */

#ifndef ANT1RELAY_H_
#define ANT1RELAY_H_

#include <itpp/itcomm.h>


class Ant1Relay {
public:
	Ant1Relay();
	virtual ~Ant1Relay();

	// Set the channel matrix
	void set_channel(itpp::cvec &ch);

	itpp::ivec bpsk_chcord_demapping(itpp::Array<std::complex<double> > &miso_out);

	itpp::ivec bpsk_nc_ml_demapping(itpp::Array<std::complex<double> > &miso_out, itpp::cvec symbols);


	// MRC-L
	itpp::ivec bpsk_mrc_pnc_demapping(itpp::Array<std::complex<double> > &miso_out);

protected:
	// 2x1 channel vector
	itpp::cvec channel;

};

inline void Ant1Relay::set_channel(itpp::cvec &ch) {
	channel = ch;
}



#endif /* ANT1RELAY_H_ */
