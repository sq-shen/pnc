/**
 *	Simulation of PNC-ZF using IT++
 *	Two antennas in the relay
 *
 */

#include <itpp/itstat.h>
#include <itpp/itcomm.h>

#include <limits>

using namespace itpp;
using namespace std;

cmat cal_pinvH(cmat &H) {
	// Get pseudo-inverse of H
	cmat herm_H = hermitian_transpose(H);
	cmat herm_H_H = herm_H * H;
	cmat inv_herm_H_H = inv(herm_H_H);
	cmat pinvH = inv_herm_H_H * herm_H;
	return pinvH;
}

double norm2_a_pinvH(vec a, cmat pinvH) {
	cvec ca = to_cvec(a);
	cmat pinvHt = transpose(pinvH);
	cvec res = pinvHt * ca;
	return energy(res);
}

int main(int argc, char *argv[]) {

	cmat read_H;
	if(argc>=2) {
		it_file ff;
		ff.open(argv[1]);
		ff>>Name("H")>>read_H;
		ff.close();	
		cout<<"Read H: "<<read_H<<endl;
	} else {
		cout<<"usage: search_a [H file]"<<endl;
		return 1;
	}
	
	// log file
	ofstream log;
	log.open("log_search_a.txt", ios::app|ios::out);
	log<<"==================================================================================================="<<endl;
	log<<"H="<<read_H<<endl;
	
	cmat pinvH = cal_pinvH(read_H);
	
	vec a1 = linspace(-20,20,41);
	vec a2 = linspace(-20,20,41);
	
	vec a, a_hat;
	a.set_size(2); a_hat.set_size(2);
	double max_val = numeric_limits<double>::min();
	double max_val2 = numeric_limits<double>::min();
	for(int i=0; i<a1.size(); i++) {
		if(a1(i)==0) continue;	
		for(int j=0; j<a2.size(); j++) {
			if(a2(j)==0) continue;
			a(0) = a1(i);
			a(1) = a2(j);
			
			double factor = norm2_a_pinvH(a, pinvH);
			double min_a = min(abs(a));
			double max_a = max(abs(a));
			
			double val = pow(min_a, 2) / factor;
			double val2  = pow(2*max_a-min_a, 2) / factor;
			
			cout<<"a="<<a<<", val="<<val<<", val2="<<val2<<endl;
			log<<"a="<<a<<", val="<<val<<", val2="<<val2<<endl;
			if(val>max_val) {
				max_val = val;
				max_val2 = val2;
				a_hat = a;
			} else if(val==max_val && val2>max_val2) {
				max_val = val;
				max_val2 = val2;
				a_hat = a;
			}
		}
	}
	
	cout<<"============================================"<<endl;
	cout<<"a_hat="<<a_hat<<", max_val="<<max_val<<", max_val2="<<max_val2<<endl;
	
	log<<"\na_hat="<<a_hat<<", max_val="<<max_val<<", max_val2="<<max_val2<<endl;
	
	log<<"===================================================================================================\n\n"<<endl;
	log.close();
	
	return 0;
}