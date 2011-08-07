/*
 * CRelay.h
 *
 *  Created on: Aug 7, 2011
 *      Author: ivan
 */

#ifndef CRELAY_H_
#define CRELAY_H_

#include <iostream>
#include <complex>


typedef struct _label {
	std::complex<double> pt;	// point in the constellation
	int bits_dec;				// mapped bits in decimal
} label;


/*
 * De-mapping boundry
 */
typedef struct _dem_region_1d {
	double low_db;	// lower bound
	double up_bd;	// upper bound
	int idx;		// 1-d index of the mapping label
} dem_region_1d;


class CRelay {
public:
	CRelay();
	virtual ~CRelay();

	//
	void init_dem_region(int x_num, dem_region_1d *x_rgns,
						 int y_num, dem_region_1d *y_rgns);

	void set_rx_sym_len(int len);

	// PNC demapping
	void pnc_demapping(int len, std::complex<double> *input, int *res_dem);

protected:

	// received symbol length
	int rx_sym_len;

	// superimposed constellation
	//int num_const_pts;
	//std::complex<label> *const_pts;

	// pnc de-mapping
	int x_dem_rgn_size;
	int y_dem_rgn_size;
	dem_region_1d *x_dem_rgn;
	dem_region_1d *y_dem_rgn;


};


inline void CRelay::set_rx_sym_len(int len) {
	rx_sym_len = len;
}

#endif /* CRELAY_H_ */
