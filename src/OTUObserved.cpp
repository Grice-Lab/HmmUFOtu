/*
 * ObservedOTU.cpp
 *
 *  Created on: Jul 11, 2017
 *      Author: zhengqi
 */

#include "OTUObserved.h"

namespace EGriceLab {

int OTUObserved::numObservedSites() const {
	return ((freq.colwise().sum() + gap).array() > 0).count();
}

int OTUObserved::numSymSites() const {
	return (freq.colwise().sum().array() > 0).count();
}

} /* namespace EGriceLab */
