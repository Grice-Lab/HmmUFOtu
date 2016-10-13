/*
 * BandedHMMP7Bg.cpp
 *
 *  Created on: May 11, 2015
 *      Author: zhengqi
 */

#include <cassert>
#include "BandedHMMP7Bg.h"

namespace EGriceLab {

void BandedHMMP7Bg::init_transPr() {
	p1 = K >= MIN_BG_K ? K / (K + 1.0) : MIN_BG_K / (MIN_BG_K + 1.0);
}

void BandedHMMP7Bg::setSize(int size) {
	K = size;
	init_transPr(); // re-init the transition probs
}

void BandedHMMP7Bg::init_bgFreq() {
	bgFreq = Vector4d::Ones() / 4.0;
}

void BandedHMMP7Bg::setBgFreq(const Vector4d& q) {
	assert((q.array() >= 0).all());
	if(q.sum() > 0)
		bgFreq = q / q.sum(); /* re-normalize */
	else
		bgFreq = Eigen::Vector4d::Ones() / 4.0; /* use all equal frequencies */
}

} /* namespace EGriceLab */
