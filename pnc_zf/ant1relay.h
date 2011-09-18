/*
 * ant1relay.h
 *
 *  Created on: Aug 13, 2011
 *      Author: ivan
 */

#ifndef ANT1RELAY_H_
#define ANT1RELAY_H_

#include <itpp/itcomm.h>

#include "sp_pt.h"

typedef struct {
	std::complex<double> pt;
	int label;
	int base_pt_label;
} proj_pt_info;


class Ant1Relay {
public:
	Ant1Relay();
	virtual ~Ant1Relay();

	void show_sp_constellation();

	// Set the channel matrix
	void set_channel(itpp::cvec &ch);

	void init_sp_const(itpp::cvec &m1, itpp::cvec &m2);


	itpp::ivec bpsk_chcord_demapping(itpp::Array<std::complex<double> > &miso_out);

	itpp::ivec bpsk_nc_ml_demapping(itpp::Array<std::complex<double> > &miso_out, itpp::cvec symbols);


	// MRC-L
	itpp::ivec bpsk_mrc_pnc_demapping(itpp::Array<std::complex<double> > &miso_out);
	
	// proj
	itpp::ivec bpsk_vecproj_pnc_demapping(itpp::Array<std::complex<double> > &miso_out);

	//-----------------------------------------
	// My proposed method

	// Initialize projected points
	itpp::Array<proj_pt_info> init_proj_pt(itpp::cvec &ch, itpp::cvec &symbols);


	// QPSK vector projection
	itpp::ivec qpsk_vecproj_pnc_demapping(itpp::Array<std::complex<double> > &miso_out, itpp::cvec &symbols);
	
	//-----------------------------------------
	// NC-ML
	itpp::ivec nc_ml_demapping(itpp::Array<std::complex<double> > &miso_out, itpp::cvec &symbols);

protected:
	// 2x1 channel vector
	itpp::cvec channel;

	itpp::Array<proj_pt_info> proj_pts;

	// Superimposed constellation points
	itpp::Array<sp_pt> sp_constellation;

};

inline void Ant1Relay::set_channel(itpp::cvec &ch) {
	channel = ch;
}



#endif /* ANT1RELAY_H_ */
