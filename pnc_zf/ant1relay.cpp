/*
 * ant1relay.cpp
 *
 *  Created on: Sep 12, 2011
 *      Author: ivan
 */

#include "ant1relay.h"

using namespace itpp;
using namespace std;

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



