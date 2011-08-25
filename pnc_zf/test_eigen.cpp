#include <itpp/itstat.h>
#include <itpp/itcomm.h>

#include "mimomac.h"
#include "zftwrelay.h"

using namespace itpp;
using namespace std;

bool ralgh_quot(cmat &H, vec &lambda, Array<cvec> &eigvec) {

	// (H*H)^{-1}
	cmat herm_H = hermitian_transpose(H);
	cmat herm_H_H = herm_H * H;
	cmat inv_herm_H_H = inv(herm_H_H);

	vec eigval;
	cmat V;
    bool res = eig_sym(inv_herm_H_H, eigval, V);

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
	
	vec lambda;
	Array<cvec> eigvec;
	bool res = ralgh_quot(H, lambda, eigvec);
	
	cout<<"res="<<res<<endl;
	cout<<"lambda="<<lambda<<endl;
	cout<<"eigvec="<<eigvec<<endl;


	return 0;
}
