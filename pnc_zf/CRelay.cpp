/*
 * CRelay.cpp
 *
 *  Created on: Aug 7, 2011
 *      Author: ivan
 */


#include "CRelay.h"

using namespace std;

CRelay::CRelay() {

	x_dem_rgn_size = 0;
	y_dem_rgn_size = 0;
	x_dem_rgn = NULL;
	y_dem_rgn = NULL;
	dem_label = NULL;

}

CRelay::~CRelay() {

	if(x_dem_rgn!=NULL) delete [] x_dem_rgn;
	if(y_dem_rgn!=NULL) delete [] y_dem_rgn;

	if(dem_label!=NULL) {
		for(int i=0; i<x_dem_rgn_size; i++) {
			delete [] dem_label[i];
		}
		delete [] dem_label;
	}	

}


void CRelay::init_dem_region(int x_size, dem_region_1d *x_rgns, 
							 int y_size, dem_region_1d *y_rgns,
							 int **dem_lbl) {

	if(x_dem_rgn != NULL) delete [] x_dem_rgn;
	if(y_dem_rgn != NULL) delete [] y_dem_rgn;
	if(dem_label != NULL) {
		for(int i=0; i<x_dem_rgn_size; i++) {
			delete [] dem_label[i];
		}
		delete [] dem_label;
	}	

	/*
	 * X demodulation boundary
	 */
	x_dem_rgn_size = x_size;
	x_dem_rgn = new dem_region_1d[x_size];
	for(int i=0; i<x_dem_rgn_size; i++) {
		x_dem_rgn[i] = x_rgns[i];
	}

	/*
	 * Y demodulation boundary
	 */
	y_dem_rgn_size = y_size;
	y_dem_rgn = new dem_region_1d[y_size];
	for(int i=0; i<y_dem_rgn_size; i++) {
		y_dem_rgn[i] = y_rgns[i];
	}


	/*
	 *	Label
	 */
	dem_label = new int*[x_dem_rgn_size];
	for(int i=0; i<x_dem_rgn_size; i++)
		dem_label[i] = new int[y_dem_rgn_size];

	for(int i=0; i<x_dem_rgn_size; i++) {
		for(int j=0; j<y_dem_rgn_size; j++) {
			dem_label[i][j] = dem_lbl[i][j];
		}
	}
}



void CRelay::pnc_demapping(int len, complex<double> *input, int *res_dem) {

	
	for(int i=0; i<len; i++) {

		int x_rgn_idx = 0, y_rgn_idx = 0;
		
		/*
	 	 *	Determine X region
	 	 */
		for(int k=0; k<x_dem_rgn_size; k++) {
			if(input[i].real()>=x_dem_rgn[k].low_bd && 
			   input[i].real()<x_dem_rgn[k].up_bd) {
					x_rgn_idx = x_dem_rgn[k].idx;	
					break;
			}
		}



		/*
	 	 *	Determine Y region
	 	 */
		for(int k=0; k<y_dem_rgn_size; k++) {
			if(input[i].imag()>=y_dem_rgn[k].low_bd && 
			   input[i].imag()<y_dem_rgn[k].up_bd) {
					y_rgn_idx = y_dem_rgn[k].idx;	
					break;
			}
		}

		res_dem[i] = dem_label[x_rgn_idx][y_rgn_idx];

		//cout<<x_rgn_idx<<","<<y_rgn_idx<<"=>"<<res_dem[i]<<endl;
	}

}



