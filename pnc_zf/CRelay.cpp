/*
 * CRelay.cpp
 *
 *  Created on: Aug 7, 2011
 *      Author: ivan
 */

#include "CRelay.h"

CRelay::CRelay() {

	x_dem_rgn_size = 0;
	y_dem_rgn_size = 0;
	x_dem_rgn = NULL;
	y_dem_rgn = NULL;

}

CRelay::~CRelay() {

	if(x_dem_rgn!=NULL) delete [] x_dem_rgn;
	if(y_dem_rgn!=NULL) delete [] y_dem_rgn;

}


void CRelay::init_dem_region(
		int x_size, dem_region_1d *x_rgns, int y_size, dem_region_1d *y_rgns) {

	/*
	 * X demodulation boundary
	 */
	if(x_dem_rgn_size != x_size) {
		if(x_dem_rgn != NULL) delete [] x_dem_rgn;
		x_dem_rgn_size = x_size;
		x_dem_rgn = new dem_region_1d[x_size];
	}

	for(int i=0; i<x_dem_rgn_size; i++) {
		x_dem_rgn[i] = x_rgns[i];
	}

	/*
	 * Y demodulation boundary
	 */
	if(y_dem_rgn_size != y_size) {
		if(y_dem_rgn != NULL) delete [] y_dem_rgn;
		y_dem_rgn_size = y_size;
		y_dem_rgn = new dem_region_1d[y_size];
	}

	for(int i=0; i<y_dem_rgn_size; i++) {
		y_dem_rgn[i] = y_rgns[i];
	}

}


