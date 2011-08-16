/*
 * zftwrelay.h
 *
 *  Created on: Aug 13, 2011
 *      Author: ivan
 */

#ifndef ZFTWRELAY_H_
#define ZFTWRELAY_H_

#include <itpp/itcomm.h>

class ZfTwRelay {
public:
	ZfTwRelay();
	virtual ~ZfTwRelay();

	//
	void setH(itpp::cmat &ch);


private:
	void init_dem_region();


protected:
	// channel
	itpp::cmat cm_H;

	// coefficient of two
	itpp::ivec iv_a;

};

inline void ZfTwRelay::setH(itpp::cmat &ch) {
	cm_H = ch;
}

#endif /* ZFTWRELAY_H_ */
