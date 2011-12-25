/*
 * sp_pt.h
 *
 *  Created on: Sep 17, 2011
 *      Author: ivan
 */

#ifndef SP_PT_H_
#define SP_PT_H_

typedef struct {
	// coordination of the point
	std::complex<double> pt;

	// decimal labeling
	int label;
} sp_pt;


#endif /* SP_PT_H_ */
