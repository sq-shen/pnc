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
		
	char *filename = NULL;
	if(argc>=2) {
		filename = argv[1];
	} else {
		filename = "H_2by2.it";
	}

	RNG_randomize();
	
    it_file ff;
    ff.open(filename);
	cmat H = randn_c(Nr, Nt);
	ff<<Name("H")<<H;
	ff.flush();
	ff.close();
	
	return 0;
}