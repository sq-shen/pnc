/*
 * ant1relay.cpp
 *
 *  Created on: Sep 12, 2011
 *      Author: ivan
 */
 
 #include "ant1relay.h"

using namespace itpp;
using namespace std;

#define DBGSHOW

#ifdef DBGSHOW
#define DBGCMD(cmd) \
	cmd;

#define SHOW1D(ary, size, type) \
{	for(int i=0; i<size; i++) \
		cout<<(type)(ary[i])<<" "; \
	cout<<endl; }

#define SHOW2D(ary, row, col, type) \
{	for(int i=0; i<row; i++){ \
		for(int j=0; j<col; j++) \
			cout<<(type)(ary[i][j])<<" "; \
		cout<<endl; \
	}}
#else
#define DBGCMD(cmd)
#define SHOW1D(ary, size, type)
#define SHOW2D(ary, row, col, type)
#endif




Ant1Relay::Ant1Relay() {
	channel.set_size(2);
}


Ant1Relay::~Ant1Relay() {

}


ivec Ant1Relay::bpsk_chcord_demapping(Array<complex<double> > &miso_out) {

	ivec res_label;
	res_label.set_size(miso_out.size());

	complex<double> h_sum = channel(0) + channel(1);
	complex<double> h_dif = channel(0) - channel(1);

	mat A = zeros(2,2);
	A(0,0) = h_sum.real();
	A(1,0) = h_sum.imag();
	A(0,1) = h_dif.real();
	A(1,1) = h_dif.imag();

	mat invA = inv(A);

	vec in_vec = zeros(2);
	for(int i=0; i<miso_out.size(); i++) {
		complex<double> in = miso_out(i);
		in_vec(0) = in.real();
		in_vec(1) = in.imag();

		vec cord = invA * in_vec;
		res_label(i) = ( abs(cord(0))>abs(cord(1)) ? 0 : 1 );
	}

	return res_label;
}

ivec Ant1Relay::bpsk_nc_ml_demapping(Array<complex<double> > &miso_out, cvec symbols) {

	ivec res_label;
	res_label.set_size(miso_out.size());

	for(int i=0; i<miso_out.size(); i++) {

		complex<double> Y = miso_out(i);
		cvec X = zeros_c(2);

		double min_norm = numeric_limits<double>::max();
		int min_pt = -1;


		for(int q=0; q<symbols.size(); q++) {
			X(0) = symbols(q);
			for(int r=0; r<symbols.size(); r++) {
				X(1) = symbols(r);
				complex<double> hx = channel * X;
				double calc_norm = norm(Y-hx);
				if(calc_norm < min_norm) {
					min_norm = calc_norm;
					min_pt = (q^r) & 0x01; // xor
				}
			}
		}

		res_label(i) = min_pt;
	}

	return res_label;
}


ivec Ant1Relay::bpsk_mrc_pnc_demapping(Array<complex<double> > &miso_output) {

	ivec res_label;
	res_label.set_size(miso_output.size());

	for(int i=0; i<miso_output.size(); i++) {

		complex<double> in = miso_output(i);

		complex<double> Ryh1 = conj(in) * channel(0);
		complex<double> Ryh2 = conj(in) * channel(1);
		complex<double> Rh1h2 = conj(channel(0)) * channel(1);

		double lhs = abs((Ryh1 + Ryh2).real());
		double rhs = abs((Ryh1 - Ryh2).real()) + 2 * Rh1h2.real();

		res_label(i) = (lhs>=rhs ? 0 : 1);
	}

	return res_label;
}


ivec Ant1Relay::bpsk_vecproj_pnc_demapping(Array<complex<double> > &miso_output) {

	ivec res_label;
	res_label.set_size(miso_output.size());
	
	int idx = 0;
	max(abs(channel), idx);
	complex<double> hmax = channel(idx);
	complex<double> hmin = channel((idx==0 ? 1 : 0));
	
	for(int i=0; i<miso_output.size(); i++) {
		
		int b1 = 0, b2 = 0;
		complex<double> in = miso_output(i);
				
		b1 = ((conj(in)*hmax).real() >=0 ? 0 : 1);
		
		complex<double> v2 = in - (1-2*b1)*hmax;
		b2 = ((conj(v2)*hmin).real() >=0 ? 0 : 1);
		
		res_label(i) = (b1 + b2) % 2;
	}
	
	
	return res_label;
}


ivec Ant1Relay::qpsk_vecproj_pnc_demapping(Array<complex<double> > &miso_output, cvec &symbols) {

	ivec res_label;
	res_label.set_size(miso_output.size());
	
	int idx = 0;
	max(abs(channel), idx);
	complex<double> hmax 	 = channel(idx);
	complex<double> hmax_r90 = hmax * polar(1.0, pi/2);
		
	complex<double> hmin  	 = channel((idx==0 ? 1 : 0));
	complex<double> hmin_r90 = hmin * polar(1.0, pi/2);
	
	DBGCMD(cout<<"hmax="<<hmax<<endl);
	DBGCMD(cout<<"hmax_r90="<<hmax_r90<<endl);
	DBGCMD(cout<<"hmin="<<hmin<<endl);
	DBGCMD(cout<<"hmin_r90="<<hmin_r90<<endl);

	for(int i=0; i<miso_output.size(); i++) {
		
		int xdir, ydir;
		int b1, b2;
		int msb, lsb;
		complex<double> r1, r2;
		complex<double> in = miso_output(i);
		
		DBGCMD(cout<<"in="<<in<<endl);

		xdir = ((conj(in)*hmax_r90).real() >=0 ? 0 : 1);
		ydir = ((conj(in)*hmax).real() >=0 ? 0 : 1);
		
		DBGCMD(cout<<"xdir="<<xdir<<", ydir="<<ydir<<endl);

		double r1_hmin, r2_hmin;
		double r1_hmin_r90, r2_hmin_r90;

		/*
		 *	first bit
		 */
		b1 = xdir;
		r1 = in - symbols(2*b1) * hmax;
		r2 = in - symbols(2*b1+1) * hmax;
		r1_hmin_r90 = (conj(r1)*hmin_r90).real();
		r2_hmin_r90 = (conj(r2)*hmin_r90).real();
		if(r1_hmin_r90>=0 && r2_hmin_r90>=0) {
			b2 = 0;
		} else if(r1_hmin_r90<0 && r2_hmin_r90<0) {
			b1 = 1;
		} else {
			r1_hmin = abs(conj(r1)*hmin);
			r2_hmin = abs(conj(r2)*hmin);
			if(r1_hmin<=r2_hmin) {
				b2 = (r1_hmin_r90>=0 ? 0 : 1);
			} else {
				b2 = (r2_hmin_r90>=0 ? 0 : 1);
			}
		}
		msb = (b1==b2 ? 0 : 1);

		DBGCMD(cout<<"r1="<<r1<<", r2="<<r2<<endl);

		
		/*
		 *	second bit
		 */
	    b1 = ydir;
		r1 = in - symbols(b1) * hmax;
		r2 = in - symbols(2+b1) * hmax;
		r1_hmin = (conj(r1)*hmin).real();
		r2_hmin = (conj(r2)*hmin).real();
		if(r1_hmin>=0 && r2_hmin>=0) {
			b2 = 0;
		} else if(r1_hmin<0 && r2_hmin<0) {
			b1 = 1;
		} else {
			r1_hmin_r90 = abs(conj(r1)*hmin_r90);
			r2_hmin_r90 = abs(conj(r2)*hmin_r90);
			if(r1_hmin_r90<=r2_hmin_r90) {
				b2 = (r1_hmin>=0 ? 0 : 1);
			} else {
				b2 = (r2_hmin>=0 ? 0 : 1);
			}
		}
		lsb = (b1==b2 ? 0 : 1);

		res_label(i) = msb*2 + lsb;

		DBGCMD(cout<<"[b1 b2]=["<<msb<<" "<<lsb<<"]"<<endl);
	}
	
	return res_label;
}


ivec Ant1Relay::nc_ml_demapping(Array<complex<double> > &miso_output, cvec &symbols) {

	ivec res_label;
	res_label.set_size(miso_output.size());

	for(int i=0; i<miso_output.size(); i++) {

		complex<double> Yr = miso_output(i);
		cvec X = zeros_c(2);

		double min_norm = numeric_limits<double>::max();
		int min_pt = -1;


		for(int q=0; q<symbols.size(); q++) {
			X(0) = symbols(q);
			for(int r=0; r<symbols.size(); r++) {
				X(1) = symbols(r);
				complex<double> HX = channel * X;
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



