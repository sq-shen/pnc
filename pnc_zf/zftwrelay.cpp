/*
 * zftwrelay.cpp
 *
 *  Created on: Aug 13, 2011
 *      Author: ivan
 */
#include <limits>
#include <itpp/itstat.h>

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

void ZfTwRelay::cal_pinvH() {
	// Get pseudo-inverse of H
	cmat herm_H = hermitian_transpose(H);
	cmat herm_H_H = herm_H * H;
	cmat inv_herm_H_H = inv(herm_H_H);
	pinvH = inv_herm_H_H * herm_H;
}

cmat ZfTwRelay::cal_mmseG(double N0) {

	cvec cv_ones = ones_c(2);
	cmat herm_H = hermitian_transpose(H);
	cmat herm_H_H = herm_H * H;
	cmat N0_herm_H_H = herm_H_H + N0*diag(cv_ones);
	cmat inv_N0_herm_H_H = inv(N0_herm_H_H);
	cmat G = inv_N0_herm_H_H * herm_H;

	return G;

}


vec ZfTwRelay::calc_opt_lincoeff() {
	
	// Get (H*H)^-1
	cmat herm_H = hermitian_transpose(H);
	cmat herm_H_H = herm_H * H;
	cmat inv_herm_H_H = inv(herm_H_H);
	
	// variables
	double p = inv_herm_H_H(0, 0).real();
	double r = inv_herm_H_H(1, 1).real();
	double s = 2 * inv_herm_H_H(0, 1).real();
	
	double min_val, min_b1, min_b2;
	double val, tmpb, b1, b2;
	
	//-----------------------------
	// Try b1 = 1
	//-----------------------------
	min_b1 = 1;
	tmpb = -s/(2*r);
	if(tmpb <= -1 || tmpb >= 1) {
		min_b2 = tmpb;
		min_val = (4*p*r - s*s) / (4*r);
	} else if(tmpb >= 0) {
		min_b2 = 1;
		min_val = p+r+s;
	} else {
		min_b2 = -1;
		min_val = p+r-s;
	}
	
	

	//-----------------------------
	// Try b1 = -1
	//-----------------------------
	b1 = -1;
	tmpb = s/(2*r);
	if(tmpb <= -1 || tmpb >= 1) {
		b2 = tmpb;
		val = (4*p*r - s*s) / (4*r);
	} else if(tmpb >= 0) {
		b2 = 1;
		val = p+r-s;
	} else {
		b2 = -1;
		val = p+r+s;
	}
	
	if(val<min_val) {
		min_val = val;
		min_b1 = b1;
		min_b2 = b2;
	}
	
	
	//-----------------------------
	// Try b2 = 1
	//-----------------------------
	b2 = 1;
	tmpb = -s/(2*p);
	if(tmpb <= -1 || tmpb >= 1) {
		b1 = tmpb;
		val = (4*p*r - s*s) / (4*p);
	} else if(tmpb >= 0) {
		b1 = 1;
		val = p+r+s;
	} else {
		b1 = -1;
		val = p+r-s;
	}
	
	if(val<min_val) {
		min_val = val;
		min_b1 = b1;
		min_b2 = b2;
	}
	
	
	//-----------------------------
	// Try b2 = -1
	//-----------------------------
	b2 = -1;
	tmpb = s/(2*p);
	if(tmpb <= -1 || tmpb >= 1) {
		b1 = tmpb;
		val = (4*p*r - s*s) / (4*p);
	} else if(tmpb >= 0) {
		b1 = 1;
		val = p+r-s;
	} else {
		b1 = -1;
		val = p+r+s;
	}
	
	if(val<min_val) {
		min_val = val;
		min_b1 = b1;
		min_b2 = b2;
	}
	

	vec a_hat;
	a_hat.set_size(2);
	a_hat(0) = min_b1;
	a_hat(1) = min_b2;
	return a_hat;
}

void ZfTwRelay::init_sp_const(itpp::cvec &comb_coeff, itpp::cvec &m1, itpp::cvec &m2) {

	// Generate the superimposed points
	int k = 0;
	cvec orin_pt = zeros_c(2);
	for(int i=0; i<m1.size(); i++) {
		orin_pt(0) = m1(i);
		for(int j=0; j<m2.size(); j++) {
			orin_pt(1) = m2(j);

			complex<double> pt = comb_coeff * ( H * orin_pt );
			int label = (i^j) & 0x03; // xor
			//cout<<i<<","<<j<<","<<label<<endl;

			sp_constellation(k).pt    = pt;
			sp_constellation(k).label = label;
			k++;
		}
	}

}

void ZfTwRelay::init_mmse_sp_const(itpp::cvec &comb_coeff, itpp::cvec &m1, itpp::cvec &m2) {

	// Generate the superimposed points
	int k = 0;
	cvec a_conj = conj(comb_coeff);
	cvec orin_pt = zeros_c(2);
	for(int i=0; i<m1.size(); i++) {
		orin_pt(0) = m1(i);
		for(int j=0; j<m2.size(); j++) {
			orin_pt(1) = m2(j);

			complex<double> pt = a_conj * orin_pt;
			int label = (i^j) & 0x03; // xor
			//cout<<i<<","<<j<<","<<label<<endl;

			sp_constellation(k).pt    = pt;
			sp_constellation(k).label = label;
			k++;
		}
	}

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
	double max_a = max(abs(a));
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
		 *	max_x = xbndry(0)
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


bool ZfTwRelay::ralgh_quot(vec &lambda, Array<cvec> &eigvec) {

	// (H*H)^{-1}
	cmat herm_H = hermitian_transpose(H);
	cmat H_herm_H = H * herm_H;
	cmat inv_H_herm_H = inv(H_herm_H);

	vec eigval;
	cmat V;
    bool res = eig_sym(inv_H_herm_H , eigval, V);

    if(!res)
    	return false;

    // increasing order
    ivec idxs = sort_index(eigval);

    lambda.set_size(eigval.size());
    eigvec.set_size(V.cols());
    for(int i=0; i<idxs.size(); i++) {
    	lambda(i) = eigval(idxs(i));
    	eigvec(i) = V.get_col(idxs(i));
    }

    return true;
}


bool ZfTwRelay::ralgh_quot(cmat &Q, vec &lambda, Array<cvec> &eigvec) {

	vec eigval;
	cmat V;
    bool res = eig_sym(Q , eigval, V);

    if(!res)
    	return false;

    // increasing order
    ivec idxs = sort_index(eigval);

    lambda.set_size(eigval.size());
    eigvec.set_size(V.cols());
    for(int i=0; i<idxs.size(); i++) {
    	lambda(i) = eigval(idxs(i));
    	eigvec(i) = V.get_col(idxs(i));
    }

    return true;
}

cvec ZfTwRelay::pnc_mmse_a(double N0) {

	// (H^*H + N_0 I)^{-1}
	int ncol = H.cols();
	cmat herm_H = hermitian_transpose(H);
	cmat P = herm_H * H;
	P += N0 * diag(ones_c(ncol));
	cmat invP = inv(P);	

	cmat Q = N0 * invP;
	vec lambda;
	Array<cvec> eigvec;
	bool res = ralgh_quot(Q, lambda, eigvec);
	if(!res)
		return false;

	cvec a;
	//a = sqrt(2) * eigvec(0);
	a = eigvec(0);
	return a;
}


/**
 * Calculate Ry once per channel change
 */
cmat ZfTwRelay::calc_Ry(double N0) {

	int ncol = H.cols();
	cmat herm_H = hermitian_transpose(H);
	cmat Ry = H * herm_H;
	Ry += N0 * diag(ones_c(ncol));

	return Ry;
}

cvec ZfTwRelay::pnc_mmse_detector(cvec &a, cmat Ry) {
	cmat iRy = inv(Ry);
	cvec w = iRy * H * a;
	cvec w_conj = conj(w);
	return w_conj;
}



ivec ZfTwRelay::pnc_ml_demapping(cvec &a, Array<itpp::cvec> &mimo_out) {

	ivec res_label;
	res_label.set_size(mimo_out.size());
	
	for(int i=0; i<mimo_out.size(); i++) {
		vec norm2 = zeros(4);
		complex<double> aYr = a * mimo_out(i);
		for(int j=0; j<sp_constellation.size(); j++) {
			complex<double> pt = sp_constellation(j).pt;
			complex<double> diff = aYr - pt;
			double calc_norm2 = diff.real() * diff.real() + diff.imag() * diff.imag();
			norm2(sp_constellation(j).label) += calc_norm2;
		}
		
		int min_pt = 0;
		double min_norm2 = min(norm2, min_pt);
		res_label(i) = min_pt;
	}

	return res_label;
}

ivec ZfTwRelay::pnc_demapping(Array<cvec> &mimo_output) {

	ivec res_label;
	res_label.set_size(mimo_output.size());
	
	// Get pseudo-inverse of H
	// cmat herm_H = hermitian_transpose(H);
	// cmat herm_H_H = herm_H * H;
	// cmat inv_herm_H_H = inv(herm_H_H);
	// cmat pinvH = inv_herm_H_H * herm_H;

	// cout<<"H="<<H<<endl;
	// cout<<"pinvH="<<pinvH<<endl;
	// cout<<"pinvH_H"<<pinvH * H<<endl;
	
	// 
	for(int i=0; i<mimo_output.size(); i++) {
		cvec in = mimo_output(i);
		complex<double> transformed_symbol = lincoeff * (pinvH * in);
		
		int idx = get_region_idx(transformed_symbol);
		res_label(i) = dem_regions(idx).label;
	}
	
	
	return res_label;
}



///////////////////////////////////////////////////////////////////////////////
//
// Traditional Network Coding
//
///////////////////////////////////////////////////////////////////////////////

ivec ZfTwRelay::nc_zf_demapping(Array<itpp::cvec> &mimo_out, QAM &qam) {

	ivec res_label;
	res_label.set_size(mimo_out.size());

	//
	for(int i=0; i<mimo_out.size(); i++) {
		cvec in = mimo_out(i);
		cvec zfed_sig = pinvH * in;
		ivec dem_sym = qam.demodulate(zfed_sig);

		res_label(i) = dem_sym(0) ^ dem_sym(1);
	}


	return res_label;
}

ivec ZfTwRelay::nc_mmse_demapping(Array<itpp::cvec> &mimo_out, double N0, QAM &qam) {

	ivec res_label;
	res_label.set_size(mimo_out.size());

	cmat G = cal_mmseG(N0);

	//
	for(int i=0; i<mimo_out.size(); i++) {
		cvec in = mimo_out(i);
		cvec zfed_sig = G * in;
		ivec dem_sym = qam.demodulate(zfed_sig);
		res_label(i) = dem_sym(0) ^ dem_sym(1);
	}


	return res_label;
}

ivec ZfTwRelay::nc_ml_demapping(itpp::Array<itpp::cvec> &mimo_out, itpp::QAM &qam) {

	ivec res_label;
	res_label.set_size(mimo_out.size());

	cvec symbols = qam.get_symbols();

	for(int i=0; i<mimo_out.size(); i++) {

		cvec Yr = mimo_out(i);
		cvec X = zeros_c(2);

		double min_norm = numeric_limits<double>::max();
		int min_pt = -1;


		for(int q=0; q<symbols.size(); q++) {
			X(0) = symbols(q);
			for(int r=0; r<symbols.size(); r++) {
				X(1) = symbols(r);
				cvec HX = H * X;
				double calc_norm = norm(Yr-HX);
				if(calc_norm < min_norm) {
					min_norm = calc_norm;
					min_pt = (q^r) & 0x03; // xor
				}
			}
		}

		res_label(i) = min_pt;
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



