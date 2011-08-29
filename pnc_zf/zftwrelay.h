/*
 * zftwrelay.h
 *
 *  Created on: Aug 13, 2011
 *      Author: ivan
 */

#ifndef ZFTWRELAY_H_
#define ZFTWRELAY_H_

#include <itpp/itcomm.h>


typedef struct {
	// coordination of the point
	std::complex<double> pt; 
	
	// decimal labeling
	int label;
} sp_pt;

typedef struct {
	double min_x;
	double max_x;
	double min_y;
	double max_y;
	int label;
} region;



class ZfTwRelay {
public:
	ZfTwRelay();
	virtual ~ZfTwRelay();

	// Set the channel matrix
	void set_H(itpp::cmat &ch);

	// Calculate mmse detector
	itpp::cmat cal_mmseG(double N0);

	// Get the pseudo-inverse of H
	itpp::cmat get_pinvH();
		
	// Initialize superimposed constellation
	void init_sp_const(itpp::cvec &comb_coeff, itpp::cvec &m1, itpp::cvec &m2);

	// Initialize superimposed constellation for pnc-mmse
	void init_mmse_sp_const(itpp::cvec &comb_coeff, itpp::cvec &m1, itpp::cvec &m2);

	// Get Rayleigh quotient; lambda are in increasing order
	bool ralgh_quot(itpp::vec &lambda, itpp::Array<itpp::cvec> &eigvec);

	// Get Rayleigh quotient; lambda are in increasing order
	bool ralgh_quot(itpp::cmat &Q, itpp::vec &lambda, itpp::Array<itpp::cvec> &eigvec);

	// Get superimposed constellation
	itpp::Array<sp_pt> get_sp_constellation();
	
	// Calculate pseudo-inverse of H
	void cal_pinvH();

	 // Set the linear combination factor
	void set_lincoeff(itpp::vec &a);
	void set_lincoeff(itpp::cvec &a);
	
	// Calculate optimized a
	itpp::vec calc_opt_lincoeff();

	// Initialize demodulation region
	void init_dem_region(itpp::vec &a, itpp::cvec &m1, itpp::cvec &m2);

	// Get the index of the demodulation region
	int get_region_idx(std::complex<double> symbol);

	//
	itpp::ivec pnc_demapping(itpp::Array<itpp::cvec> &mimo_out);
	
	/*
	 *	PNC-MMSE functions
	 */
	// Get the optimized factor of pnc-mmse linear combination
	itpp::cvec pnc_mmse_a(double N0);

	// Calculate R_y
	itpp::cmat calc_Ry(double N0);

	// Get PNC-MMSE detector
	itpp::cvec pnc_mmse_detector(itpp::cvec &a, itpp::cmat Ry);
	
	/*
	 *	PNC-ML
	 */
	//
	itpp::ivec pnc_ml_demapping(itpp::cvec &a, itpp::Array<itpp::cvec> &mimo_out);

	
	/*
	 *	NC functions
	 */
	//
	itpp::ivec nc_zf_demapping(itpp::Array<itpp::cvec> &mimo_out, itpp::QAM &qam);

	//
	itpp::ivec nc_mmse_demapping(itpp::Array<itpp::cvec> &mimo_out, double N0, itpp::QAM &qam);

	//
	itpp::ivec nc_ml_demapping(itpp::Array<itpp::cvec> &mimo_out, itpp::QAM &qam);

	/*
	 *	Testing functions
	 */
public:
	void show_sp_constellation();
	void show_dem_regions();
		

protected:
	// Estimated channel matrix
	itpp::cmat H;
	
	// Pseudo-inverse of H
	itpp::cmat pinvH;

	// Linear combination factor of the two transmitted signals
	itpp::cvec lincoeff;

	// Superimposed constellation points
	itpp::Array<sp_pt> sp_constellation;

	// Demodulation region (9 regions)
	itpp::Array<region> dem_regions;

	itpp::vec xbndry;
	itpp::vec ybndry;

};

inline void ZfTwRelay::set_H(itpp::cmat &ch) {
	H = ch;
	cal_pinvH();
}

inline itpp::cmat ZfTwRelay::get_pinvH() {
	return pinvH;
}

inline void ZfTwRelay::set_lincoeff(itpp::cvec &a) {
	lincoeff = a;
}

inline void ZfTwRelay::set_lincoeff(itpp::vec &a) {
	lincoeff.set_size(a.size());
	for(int i=0; i<a.size(); i++) {
		lincoeff(i) = a(i);
	}
}

inline itpp::Array<sp_pt> ZfTwRelay::get_sp_constellation() {
	return sp_constellation;
}


#endif /* ZFTWRELAY_H_ */
