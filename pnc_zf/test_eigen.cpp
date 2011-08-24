#include <itpp/itstat.h>
#include <itpp/itcomm.h>

#include "mimomac.h"
#include "zftwrelay.h"

using namespace itpp;
using namespace std;

int main(int argc, char *argv[]) 
{
	// channel matrix from file
	cmat H;
	if(argc>=2) {
		it_file ff;
		ff.open(argv[1]);
		ff>>Name("H")>>H;
		ff.close();	
		cout<<"H: "<<H<<endl;
	} else {
		cout<<"test_eigen [H file]"<<endl;
		return 1;
	}
	
	
	// (H*H)^{-1}
	cmat herm_H = hermitian_transpose(H);
	cmat herm_H_H = herm_H * H;
	cmat inv_herm_H_H = inv(herm_H_H);
	
	cmat V;
    vec d;
    bool res = eig_sym(inv_herm_H_H, d, V);
	
	cout<<"res="<<res<<endl;
	cout<<"d="<<d<<endl;
	cout<<"V="<<V<<endl;
	
	double eigval = d(0);
	cvec eigvec = V.get_col(0);
	cout<<"Av="<<round_to_zero(inv_herm_H_H*eigvec)<<endl;
	cout<<"lambda v ="<<round_to_zero(eigval*eigvec)<<endl;
	
	
	
	

	return 0;
}