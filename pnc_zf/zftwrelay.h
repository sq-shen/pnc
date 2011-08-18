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

	// Set the linear combination factor
	void set_lincoeff(itpp::vec &a);
	void set_lincoeff(itpp::cvec &a);


	// Initialize demodulation region
	void init_dem_region(itpp::vec &a, itpp::cvec &m1, itpp::cvec &m2);

	// Get the index of the demodulation region
	int get_region_idx(std::complex<double> symbol);

	//
	itpp::ivec pnc_demapping(itpp::Array<itpp::cvec> &mimo_out);
	
	/*
	 *	Testing functions
	 */
public:
	void show_sp_constellation();
	void show_dem_regions();
		

protected:
	// Estimated channel matrix
	itpp::cmat H;

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


#endif /* ZFTWRELAY_H_ */
