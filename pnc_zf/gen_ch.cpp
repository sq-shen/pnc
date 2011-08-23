/**
 *	Generate a Nr x Nt channel matrix
 *
 */

#include <itpp/itcomm.h>

using namespace itpp;
using namespace std;

int main(int argc, char *argv[]) {

	int Nt = 2, 	// # of transmit antenna
	    Nr = 2;		// # of receive antenna
		
	char filename[256];
	if(argc>=2) {
		sprintf(filename, "%s", argv[1]);
	} else {
		sprintf(filename, "H_2by2.it");
	}

	RNG_randomize();
	
    it_file ff;
    ff.open(filename);
	cmat H = randn_c(Nr, Nt);
	ff<<Name("H")<<H;
	ff.flush();
	ff.close();
	
	// read back
	cmat read_H;
    ff.open(filename);
	ff>>Name("H")>>read_H;
	ff.close();
	
	cout<<"Read back:\n"<<read_H<<endl;
	
	return 0;
}