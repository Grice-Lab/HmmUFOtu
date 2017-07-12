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
 * ObservedOTU.h
 *  An observed OTU type, that is a PTUNodePtr with additional observed data
 *  Created on: Jul 11, 2017
 *      Author: zhengqi
 */

#ifndef SRC_OTUOBSERVED_H_
#define SRC_OTUOBSERVED_H_

#include <string>
#include <Eigen/Dense>

namespace EGriceLab {
using std::string;
using Eigen::Matrix4Xd;
using Eigen::RowVectorXd;

struct OTUObserved {
	/** constructors */
	/** default constructor, do nothing */
	OTUObserved() {  }

	/** construct an OTUObserved with given information */
	OTUObserved(const string& id, const string& taxon, int csLen, int N) :
		id(id), taxon(taxon), csLen(csLen), N(N), freq(4, csLen), gap(csLen), count(N)
	{
		freq.setZero();
		gap.setZero();
		count.setZero();
	}

	virtual ~OTUObserved() {  }

	/** member methods */
	/** get observed number of reads */
	double numReads() const {
		return count.sum();
	}

	/** get observed number of samples */
	int numSamples() const {
		return (count.array() > 0).count();
	}

	/** get observed site number across all samples */
	int numObservedSites() const;

	/** get observed site fraction */
	double fracObservedSites() const {
		return numObservedSites() / static_cast<double>(csLen);
	}

	/** get non-gap symbol site number */
	int numSymSites() const;

	/** get non-gap symbol site fraction */
	double fracSymSites() const {
		return numSymSites() / static_cast<double> (csLen);
	}

	string id; /* id for this OTU */
	string taxon; /* taxon for this OTU */
	int csLen;  /* consensus sequence length */
	int N;      /* number of total samples */
	Matrix4Xd freq;  /* observed aggregate base frequency over all samples */
	RowVectorXd gap; /* observed aggregate gap over all samples */
	RowVectorXd count;  /* observed sequence count for each sample separately */
};

} /* namespace EGriceLab */

#endif /* SRC_OTUOBSERVED_H_ */
