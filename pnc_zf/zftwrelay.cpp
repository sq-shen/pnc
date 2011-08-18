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
	double max_a = abs(itpp::max(a));
	xbndry(0) = -max_a * sqrt(2)/2;
	xbndry(1) =  max_a * sqrt(2)/2;
	ybndry(0) = -max_a * sqrt(2)/2;
	ybndry(1) =  max_a * sqrt(2)/2;


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
		//cout<<idx<<","<<sp_constellation(i).label<<endl;
	}
}



int ZfTwRelay::get_region_idx(std::complex<double> symbol) {
	int x, y;

	if(symbol.real() <= xbndry(0))
		x = 0;
	else if(symbol.real() > xbndry(0) && symbol.real() <= xbndry(1))
		x = 1;
	else
		x = 2;

	if(symbol.imag() <= ybndry(0))
		y = 0;
	else if(symbol.imag() > ybndry(0) && symbol.imag() <= ybndry(1))
		y = 1;
	else
		y = 2;
	//cout<<"symbol="<<symbol<<","<<x<<","<<y<<endl;
	return 3*y + x;
}


ivec ZfTwRelay::pnc_demapping(Array<cvec> &mimo_output) {

	ivec res_label;
	res_label.set_size(mimo_output.size());
	
	// Get pseudo-inverse of H
	cmat herm_H = hermitian_transpose(H);
	cmat herm_H_H = herm_H * H;
	cmat inv_herm_H_H = inv(herm_H_H);
	cmat pinv_H = inv_herm_H_H * herm_H;

	// cout<<"H="<<H<<endl;
	// cout<<"pinv_H="<<pinv_H<<endl;
	// cout<<"pinv_H_H"<<pinv_H * H<<endl;
	
	// 
	for(int i=0; i<mimo_output.size(); i++) {
		cvec in = mimo_output(i);
		complex<double> transformed_symbol = lincoeff * (pinv_H * in);
		
		int idx = get_region_idx(transformed_symbol);
		res_label(i) = dem_regions(idx).label;
	}
	
	
	return res_label;
}


void ZfTwRelay::show_sp_constellation() {
	cout<<"sp_constellation: "<<endl;
	for(int i=0; i<sp_constellation.size(); i++) {
		cout<<"pt="<<sp_constellation(i).pt<<",label="<<sp_constellation(i).label<<endl;
	}
}

void ZfTwRelay::show_dem_regions() {
	for(int i=0; i<dem_regions.size(); i++) {
		cout<<"label:"<<dem_regions(i).label<<", "
		    <<"("<<dem_regions(i).min_x<<","<<dem_regions(i).max_x<<"), "
			<<"("<<dem_regions(i).min_y<<","<<dem_regions(i).max_y<<")"<<endl;
	}
}



