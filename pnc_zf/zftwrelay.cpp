/*
 * zftwrelay.cpp
 *
 *  Created on: Aug 13, 2011
 *      Author: ivan
 */

#include "zftwrelay.h"

using namespace itpp;
using namespace std;


ZfTwRelay::ZfTwRelay() {
	sp_constellation.set_size(16);
}

ZfTwRelay::~ZfTwRelay() {
}


void ZfTwRelay::init_dem_region(vec &a, itpp::cvec &m1, itpp::cvec &m2) {
	
	// TODO: Assertion
	// m1.size()==2 && m2.size()==2
	
	// Set linear combination coefficient
	set_lincoeff(a);

	// Generate the superimposed points
	int k = 0;
	for(int i=0; i<m1.size(); i++) {
		for(int j=0; j<m2.size(); j++) {
			complex<double> pt = m1(i) + m2(j);
			int label = (i^j) & 0x03; // xor
			//cout<<i<<","<<j<<","<<label<<endl;
			
			sp_constellation(k++).pt    = pt;
			sp_constellation(k++).label = label;
		}
	}
	


}
