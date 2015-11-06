/*
 * BandedHMMP7Bg.cpp
 *
 *  Created on: May 11, 2015
 *      Author: zhengqi
 */

#include <cassert>
#include "BandedHMMP7Bg.h"

namespace EGriceLab {

} /* namespace EGriceLab */

void EGriceLab::BandedHMMP7Bg::init_transPr() {
	if(K > 0)
		transGG = static_cast<float>(K) / (K + 1);
	else
		transGG = 0;
}

void EGriceLab::BandedHMMP7Bg::setSize(int size) {
	K = size;
	init_transPr(); // re-init the transition probs
}

void EGriceLab::BandedHMMP7Bg::init_bgFreq() {
	bgFreq = VectorXf::Ones(nuclAbc->getSize()); // all equal frequencies
	bgFreq /= bgFreq.sum();
}

void EGriceLab::BandedHMMP7Bg::setBgFreq(const VectorXf& q) {
	assert(q.size() == bgFreq.size() && q.minCoeff() >= 0 && q.sum() > 0);
	bgFreq = q / q.sum(); // re-normalize, even if already normalized
}
