/*******************************************************************************
 * This file is part of HmmUFOtu, an HMM and Phylogenetic placement
 * based tool for Ultra-fast taxonomy assignment and OTU organization
 * of microbiome sequencing data with species level accuracy.
 * Copyright (C) 2017  Qi Zheng
 *
 * HmmUFOtu is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HmmUFOtu is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
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
