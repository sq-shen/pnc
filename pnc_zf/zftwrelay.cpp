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
	dem_regions.set_size(9);
	xbndry.set_size(2);
	ybndry.set_size(2);
}

ZfTwRelay::~ZfTwRelay() {
}


void ZfTwRelay::init_dem_region(vec &a, itpp::cvec &m1, itpp::cvec &m2) {
	
	// TODO: Assertion
	// a.size()==2 && m1.size()==2 && m2.size()==2
	
	// Set linear combination coefficient
	set_lincoeff(a);

	// Generate the superimposed points
	int k = 0;
	for(int i=0; i<m1.size(); i++) {
		for(int j=0; j<m2.size(); j++) {
			complex<double> pt = a(0)*m1(i) + a(1)*m2(j);
			int label = (i^j) & 0x03; // xor
			//cout<<i<<","<<j<<","<<label<<endl;
			
			sp_constellation(k).pt    = pt;
			sp_constellation(k).label = label;
			k++;
		}
	}
	
	/*
	 * Get the 9 demodulation regions
	 *
	 *  6 | 7 | 8
	 * ---+---+---
	 *  3 | 4 | 5
	 * ---+---+---
	 *  0 | 1 | 2
	 */
	double max_a = itpp::max(a);
	xbndry(0) = -max_a * sqrt(2);
	xbndry(1) =  max_a * sqrt(2);
	ybndry(0) = -max_a * sqrt(2);
	ybndry(1) =  max_a * sqrt(2);


	//
	//	Not necessary ?
	//
	for(int i=0; i<3; i++) {

		/*
		 *	min_x = numeric_limits<double>::min()
		 *	man_x = xbndry(0)
		 *
		 *	min_y = numeric_limits<double>::min()
		 *	max_y = ybndry(0)
		 */
		dem_regions(3*i).min_x = numeric_limits<double>::min();
		dem_regions(3*i).max_x = xbndry(0);
		dem_regions(i).min_y = numeric_limits<double>::min();
		dem_regions(i).max_y = ybndry(0);

		/*
		 *	min_x = xbndry(0)
		 *	man_x = xbndry(1)
		 *
		 *	min_y = ybndry(0)
		 *	max_y = ybndry(1)
		 */
		dem_regions(3*i+1).min_x = xbndry(0);
		dem_regions(3*i+1).max_x = xbndry(1);
		dem_regions(3+i).min_y = ybndry(0);
		dem_regions(3+i).max_y = ybndry(1);

		/*
		 *	min_x = xbndry(1)
		 *	man_x = numeric_limits<double>::max()
		 *
		 *	min_y = ybndry(1)
		 *	max_y = numeric_limits<double>::max()
		 */
		dem_regions(3*i+2).min_x = xbndry(1);
		dem_regions(3*i+2).max_x = numeric_limits<double>::max();
		dem_regions(6+i).min_y = ybndry(1);
		dem_regions(6+i).max_y = numeric_limits<double>::max();
	}

	// Get the labeling for each region
	for(int i=0; i<sp_constellation.size(); i++) {
		int idx = get_region_idx( sp_constellation(i).pt );
		dem_regions(idx).label = sp_constellation(i).label;
	}
}



int ZfTwRelay::get_region_idx(std::complex<double> symbol) {
	int x, y;

	if(symbol.real() <= xbndry(0))
		x = 0;
	else if(symbol.real() > xbndry(0) && symbol.real() < xbndry(1))
		x = 1;
	else
		x = 2;

	if(symbol.imag() <= ybndry(0))
		y = 0;
	else if(symbol.imag() > ybndry(0) && symbol.imag() < ybndry(1))
		y = 1;
	else
		y = 2;

	return 3*y + x;
}


ivec ZfTwRelay::pnc_demapping(itpp::cvec &recv_signal) {

	ivec result;
	result.set_size(recv_signal.size());

	for(int i=0; i<recv_signal.size(); i++) {
		int idx = get_region_idx(recv_signal(i));
		result(i) = dem_regions(idx).label;
	}

	return result;
}



