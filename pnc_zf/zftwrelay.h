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



class ZfTwRelay {
public:
	ZfTwRelay();
	virtual ~ZfTwRelay();

	// Set the channel matrix
	void set_H(itpp::cmat &ch);

	// Set the linear combination factor
	void set_lincoeff(itpp::vec &a);


	// Initialize demodulation region
	void init_dem_region(itpp::vec &a, itpp::cvec &m1, itpp::cvec &m2);

	//  
		

protected:
	// Estimated channel matrix
	itpp::cmat H;

	// Linear combination factor of the two transmitted signals
	itpp::vec lincoeff;

	// Superimposed constellation points
	itpp::Array<sp_pt> sp_constellation;

};

inline void ZfTwRelay::set_H(itpp::cmat &ch) {
	H = ch;
}

inline void ZfTwRelay::set_lincoeff(itpp::vec &a) {
	lincoeff = a;
}


#endif /* ZFTWRELAY_H_ */
