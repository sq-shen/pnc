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
	
	equiv_noise.set_size(2);
	
	mapper = NULL;
}

ZfTwRelay::~ZfTwRelay() {
}

void ZfTwRelay::set_H(itpp::cmat &ch) {
	H = ch;
	cal_pinvH();

	// MIMO-PNC
	calc_mimopnc_H();
}

void ZfTwRelay::cal_pinvH() {
	// Get pseudo-inverse of H
	cmat herm_H = hermitian_transpose(H);
	cmat herm_H_H = herm_H * H;
	cmat inv_herm_H_H = inv(herm_H_H);
	pinvH = inv_herm_H_H * herm_H;
}

void ZfTwRelay::cal_pinvH(cmat &pH) {
	// Get pseudo-inverse of H
	cmat herm_H = hermitian_transpose(H);
	cmat herm_H_H = herm_H * H;
	cmat inv_herm_H_H = inv(herm_H_H);
	pH = inv_herm_H_H * herm_H;
}

cmat ZfTwRelay::cal_mmseG(double N0) {

	cvec cv_ones = ones_c(2);
	cmat herm_H = hermitian_transpose(H);
	cmat herm_H_H = herm_H * H;
	cmat N0_herm_H_H = herm_H_H + N0*diag(cv_ones);
	cmat inv_N0_herm_H_H = inv(N0_herm_H_H);
	cmat mmseG = inv_N0_herm_H_H * herm_H;
	return mmseG;
}


vec ZfTwRelay::calc_opt_lincoeff(int type, double N0) {
	
	// Get G
	//	ZF:   G = (H^*H)^-1
	//	MMSE: G = (H^*H+ N0I)^-1
	if(type==0) {
		cal_pinvH(G);
	} else {
		G = cal_mmseG(N0);
	}
	
	cmat G2 = G * hermitian_transpose(G);
	
	// variables
	double p = G2(0, 0).real();
	double r = G2(1, 1).real();
	double s = 2 * G2(0, 1).real();
	
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
	if(!res) {
		return ones_c(2);
	}

	cvec a;
	a = 0.5 * eigvec(0);
	//a = eigvec(0);
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
		double min_norm2 = numeric_limits<double>::max();
		int min_pt = 0;
		
		complex<double> aYr = a * mimo_out(i);
		for(int j=0; j<sp_constellation.size(); j++) {
			complex<double> pt = sp_constellation(j).pt;
			complex<double> diff = aYr - pt;
			double calc_norm2 = diff.real() * diff.real() + diff.imag() * diff.imag();
			
			if(calc_norm2<min_norm2) {
				min_norm2 = calc_norm2;
				min_pt = sp_constellation(j).label;
			}
		}
		
		res_label(i) = min_pt;
	}

	return res_label;
}

ivec ZfTwRelay::pnc_demapping(int type, Array<cvec> &mimo_output, double N0) {

	ivec res_label;
	res_label.set_size(mimo_output.size());
	
	for(int i=0; i<mimo_output.size(); i++) {
		cvec in = mimo_output(i);
		complex<double> transformed_symbol = lincoeff * (G * in);
		
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

/*
 *	MIMO-PNC
 */
void ZfTwRelay::calc_mimopnc_H() {
	cmat invD = "0.5 0.5; 0.5 -0.5";
	mimopnc_H = H * invD;
}

//
//
//
void ZfTwRelay::calc_mimopnc_G(int type, double N0) {
	
	cmat herm_H = hermitian_transpose(mimopnc_H);
	cmat herm_H_H = herm_H * mimopnc_H;
	
	if(type==0) {
		// Get pseudo-inverse of H	
		cmat inv_herm_H_H = inv(herm_H_H);
		mimopnc_G = inv_herm_H_H * herm_H;
	} else {
		cvec cv_ones = ones_c(2);
		cmat N0_herm_H_H = herm_H_H + N0*diag(cv_ones);
		cmat inv_N0_herm_H_H = inv(N0_herm_H_H);
		mimopnc_G = inv_N0_herm_H_H * herm_H;	
	}
	
	cmat hermG = hermitian_transpose(mimopnc_G);

	//cmat hermG_G = hermG * mimopnc_G;	// paper ??
	cmat hermG_G = mimopnc_G * hermG;

	equiv_noise(0) = hermG_G(0,0).real() * N0;
	equiv_noise(1) = hermG_G(1,1).real() * N0;
	//cout<<hermG_G(0,0)<<", "<<hermG_G(1,1)<<endl;
	//cout<<"equiv_noise="<<equiv_noise<<endl;

}


ivec ZfTwRelay::mimo_pnc_demapping(itpp::Array<itpp::cvec> &mimo_out) {
	
	ivec res_label = zeros_i(mimo_out.size());
	
	int bits_per_symbol = mapper->bits_per_symbol();
	cvec symbols = mapper->get_symbols();
	double upbd = 2 * abs(symbols(0).real());
	double lobd = -upbd;
	
	vec y1 = zeros(2);
	vec y2 = zeros(2);
	
	double n1 = equiv_noise(0)/2;
	double n2 = equiv_noise(1)/2;
	
	for(int i=0; i<mimo_out.size(); i++) {
		cvec r  = mimo_out(i);
		cvec Gr = mimopnc_G * r;
		
		y1(0) = Gr(0).real();
		y1(1) = Gr(0).imag();
		
		y2(0) = Gr(1).real();
		y2(1) = Gr(1).imag();
			
		
		ivec lbl = zeros_i(bits_per_symbol);		
		for(int k=0; k<bits_per_symbol; k++) {
			
			double numerator, denominator, LR;
			
			numerator   = ( exp( -pow(y1(k)-upbd, 2) / (2*n1) ) + exp( -pow(y1(k)-lobd, 2) / (2*n1) ) )* 
						  exp( -pow(y2(k), 2) / (2*n2) );
			
			denominator = exp(-pow(y1(k), 2) / (2*n1)) *
						  ( exp(-pow(y2(k)-upbd, 2) / (2*n2)) + exp(-pow(y2(k)-lobd, 2) / (2*n2)) );
			
			LR = numerator / denominator;
			lbl(k) = (LR>=1 ? 0 : 1);
			
			res_label(i) += lbl(k) << k;
		}
		//res_label(i) = lbl;
	}

	return res_label;
}

// ivec ZfTwRelay::mimo_pnc_demapping(itpp::Array<itpp::cvec> &mimo_out) {
	
	// ivec res_label;
	// res_label.set_size(mimo_out.size());
	// double var1 = equiv_noise(0)/2;
	// double var2 = equiv_noise(1)/2;
	
	// for(int i=0; i<mimo_out.size(); i++) {
		// cvec r  = mimo_out(i);
		// cvec Gr = mimopnc_G * r;
		// complex<double> y1 = Gr(0);
		// complex<double> y2 = Gr(1);
		
		// ivec lbl = "0 0";
		// //int lbl = 0 ;
		
		// // Real part
		// double L = exp(2/var2-2/var1) * cosh(2*y1.real()/var1) / cosh(2*y2.real()/var2);
		// lbl(0) = (L>=1 ? 0 : 1);
		// lbl = (L>=1 ? 0 : 1);
		
		// // Imaginary part 
		// L = exp(2/var2-2/var1) * cosh(2*y1.imag()/var1) / cosh(2*y2.imag()/var2);
		// lbl(1) = (L>=1 ? 0 : 1);
		
		
		// res_label(i) = 2*lbl(1) + lbl(0);
		// //res_label(i) = lbl;
	// }

	// return res_label;
// }
