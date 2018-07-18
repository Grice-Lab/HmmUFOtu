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
 * OTUTable.h
 *  An OTU Table represent the summary of OTU abundance over multiple samples
 *  Created on: Jul 11, 2017
 *      Author: zhengqi
 */

#ifndef SRC_OTUTABLE_H_
#define SRC_OTUTABLE_H_

#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <boost/algorithm/string.hpp> /* for boost::split */
#include <boost/random/mersenne_twister.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/random/discrete_distribution.hpp>
#include "ProgLog.h"
#include "OTUObserved.h"

namespace EGriceLab {
namespace HmmUFOtu {

using std::string;
using std::vector;
using std::map;
using std::invalid_argument;
using std::find;
using std::istream;
using std::ostream;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

class OTUTable {
public:
	/* typedefs */
	typedef map<string, string> otuMap;
	typedef boost::random::mt11213b RNG; /* preferred random number generator type */
	typedef boost::random::discrete_distribution<size_t> ReadDistrib; /* base (nucleotide) distribution */
	typedef ReadDistrib::param_type ReadParam; /* read distribution parameters */

	/** constructors */
	/** default constructor */
	OTUTable() {  }

	/** construct an OTUTable with given samples and OTU list */
	OTUTable(const vector<string>& samples, const vector<string>& otus, const otuMap& otu2Taxon, const MatrixXd& otuMetric) :
		samples(samples), otus(otus), otu2Taxon(otu2Taxon), otuMetric(otuMetric)
	{  }

	/** construct an OTUTable with initial samples only */
	explicit OTUTable(const vector<string>& samples) :
			samples(samples), otuMetric(0, samples.size())
	{  }

	/** destructor, do nothing */
	virtual ~OTUTable() { }

	/** member methods */
	/** test whether this OTU table is empty */
	bool isEmpty() const {
		return numOTUs() == 0 && numSamples() == 0;
	}

	/** get number of samples */
	size_t numSamples() const {
		return samples.size();
	}

	/** get number of OTUs */
	size_t numOTUs() const {
		return otus.size();
	}

	/** get all sample names */
	const vector<string>& getSamples() const {
		return samples;
	}

	/** get all OTU names */
	const vector<string>& getOTUs() const {
		return otus;
	}

	/** get sample at given position */
	const string& getSample(string::size_type j) const {
		return samples[j];
	}

	/** get OTU at given position */
	const string& getOTU(string::size_type i) const {
		return otus[i];
	}

	/** get entire OTU metric count */
	const MatrixXd& getMetric() const {
		return otuMetric;
	}

	/**
	 * test whether this OTUTable is relative abundance
	 * @return  true only if all metrics are no greater than 1
	 */
	bool isRelative() const {
		return (otuMetric.array() <= 1.0).all();
	}

	/**
	 * test whether this OTUTable contains a specific sample
	 */
	bool hasSample(const string& sampleName) const {
		return find(samples.begin(), samples.end(), sampleName) != samples.end();
	}

	/**
	 * test whether this OTUTable contains a specific OTU
	 */
	bool hasOTU(const string& otuID) const {
		return find(otus.begin(), otus.end(), otuID) != otus.end();
	}

	/**
	 * get the index of a given sample
	 * @return  0..numSamples-1 if found, or numSamples if not
	 */
	size_t getSampleIndex(const string& sampleName) const {
		return find(samples.begin(), samples.end(), sampleName) - samples.begin();
	}

	/**
	 * get the index of a given OTU
	 * @return  0..numOTUs-1 if found, or numOTUs if not
	 */
	size_t getOTUIndex(const string& otuID) const {
		return find(otus.begin(), otus.end(), otuID) - otus.begin();
	}

	/**
	 * get taxon of given otuID
	 */
	const string& getTaxon(const string& otuID) const {
		return otu2Taxon.at(otuID);
	}

	/** get taxon of given OTU index */
	const string& getTaxon(size_t i) const {
		return getTaxon(getOTU(i));
	}

	/** get count of given sample index */
	VectorXd numSampleMetric(size_t j) const {
		return otuMetric.col(j);
	}

	/** get count of given sample name */
	VectorXd numSampleMetric(const string& sn) const {
		return numSampleMetric(getSampleIndex(sn));
	}

	/** get count of given OTU index */
	RowVectorXd numOTUMetric(size_t i) const {
		return otuMetric.row(i);
	}

	/** get count of given OTU */
	RowVectorXd numOTUMetric(const string& otuID) const {
		return numOTUMetric(getOTUIndex(otuID));
	}

	/** get count of given OTU and sample idx */
	double numOTUMetric(size_t i, size_t j) const {
		return otuMetric(i, j);
	}

	/** get count of given OTU and sample */
	double numOTUMetric(const string& otuID, const string& sn) const {
		return numOTUMetric(getOTUIndex(otuID), getSampleIndex(sn));
	}

	/**
	 * get number of reads in sample j
	 */
	double sumSampleMetric(size_t j) const {
		return otuMetric.col(j).sum();
	}

	/**
	 * get number of reads in OTU i
	 */
	double sumOTUMetric(size_t i) const {
		return otuMetric.row(i).sum();
	}

	/**
	 * get total count of a given sample across all OTUs
	 * @return  column-sum of given sample, or undefined behavior if not found
	 */
	double sumSampleMetric(const string& sn) const {
		return sumSampleMetric(getSampleIndex(sn));
	}

	/**
	 * get total count of an OTU across all samples
	 * @return  row-sum of given otuID, or undefined behavior if not found
	 */
	double sumOTUMetric(const string& otuID) const {
		return sumOTUMetric(getOTUIndex(otuID));
	}

	/**
	 * add a new sample into this OTUTable, ignored if already exists
	 * @param sampleName  new sample name
	 * @return  true if this sampleName is not already exist
	 */
	bool addSample(const string& sampleName);

	/** delete an existing sample from this OTUTable, ignored if not exists
	 * @param sampleName  existing sample name
	 * @return  true only if this sampleName exists
	 */
	bool removeSample(const string& sampleName);


	/** delete an existing sample from this OTUTable using its index
	 * @param j  sample index
	 * @return  true only if the index is in range
	 */
	bool removeSample(size_t j);

	/**
	 * add a new OTU into this OTUTable, ignored if already exists
	 * @param otuID  new OTU ID
	 * @param taxon  taxon for new OTU
	 * @param count  count row vector for new OTU
	 * @return  true if this OTUID is not already exist
	 */
	bool addOTU(const string& otuID, const string& taxon, const RowVectorXd& count);

	/**
	 * add a new OTU into this OTUTable, ignored if already exists
	 * @param otuID  new OTU ID
	 * @param taxon  taxon for new OTU
	 * @return  true if this OTUID is not already exist
	 */
	bool addOTU(const string& otuID, const string& taxon = "") {
		return addOTU(otuID, taxon, RowVectorXd::Zero(numSamples()));
	}

	/**
	 * add an OTUObserved into this OTUTable, ignore if already exists
	 * @param otu  new OTUObserved
	 * @return  true if this otu is not already exist
	 */
	bool addOTU(const OTUObserved& otu) {
		return addOTU(otu.id, otu.taxon, otu.count);
	}

	/** delete an existing otuID from this OTUTable at index i, ignored if outside range
	 * @param i  OTU index
	 * @return  true only if this OTU index is in range
	 */
	bool removeOTU(size_t i);

	/** delete an existing otuID from this OTUTable, ignored if not exists
	 * @param otuID  existing OTU ID
	 * @return  true only if this otuID exists
	 */
	bool removeOTU(const string& otuID);

	/**
	 * clear this entire OTUTable to empty
	 */
	void clear();

	/**
	 * prune bad samples with less than min reads, usually after calling subset
	 */
	void pruneSamples(size_t min = 0);

	/**
	 * prune bad OTUs with less than min reads, usually after calling subset
	 */
	void pruneOTUs(size_t min = 0);

	/**
	 * normalize the metric with constant method
	 * @param Z  normalization constant
	 */
	void normalizeConst(double Z = 0);

	/**
	 * normalize the metric
	 * @param Z  normalization constant
	 * @param method  normalization method
	 */
	void normalize(double Z = 0, const string& method = "constant") {
		if(method == "constant")
			normalizeConst(Z);
		else
			throw invalid_argument("Unsupported subsetting method '" + method + "'");
	}

	/**
	 * set seed for subset functions
	 */
	void seed(unsigned newSeed) {
		rng.seed(newSeed);
	}

	/**
	 * subset this OTU table to a minimum read count using given method
	 * samples that have less than min reads will be removed
	 * samples that have more than min reads will be subsampled
	 * @param min  min read requirement
	 * @param method  sampleing method
	 * @throw  invalid_argument if the sampling method is not supported
	 */
	void subset(size_t min, const string& method) {
		if(method == "uniform")
			subsetUniform(min);
		else if(method == "multinomial")
			subsetMultinom(min);
		else
			throw invalid_argument("Unsupported subsetting method '" + method + "'");
	}

	/**
	 * subset this OTU table to a minimum read count using uniform sampling
	 */
	void subsetUniform(size_t min);

	/**
	 * subset this OTU table to a minimum read count using Multinomial sampling
	 */
	void subsetMultinom(size_t min);

	/**
	 * load raw table object from input in given format
	 */
	istream& load(istream& in, const string& format = "table");

	/**
	 * save this table to an output stream in given format
	 */
	ostream& save(ostream& out, const string& format = "table") const;

	/**
	 * load raw table object from input in table format
	 */
	istream& loadTable(istream& in);

	/**
	 * save this table to an output stream in table format
	 */
	ostream& saveTable(ostream& out) const;

	/**
	 * load a table from an input stream in BIOM hdf5 format
	 */
	istream& loadHdf5(istream& in);

	/**
	 * save this table to an output stream in BIOM hdf5 format
	 */
	ostream& saveHdf5(ostream& out) const;

	/** merge this OTUTable with another */
	OTUTable& operator+=(const OTUTable& other);

	/* non-member functions */
	friend OTUTable operator+(const OTUTable& lhs, const OTUTable& rhs);

private:
	/** member fields */
	vector<string> samples; /* 0..N sample names */
	vector<string> otus;    /* 0..M OTUs */
	otuMap otu2Taxon;
	MatrixXd otuMetric; /* M * N matrix of OTU (relative) abundance metric */

	/** static fields */
	static const Eigen::IOFormat dblTabFmt;
	static const Eigen::IOFormat fltTabFmt;
	static RNG rng;
};

inline std::istream& OTUTable::load(istream& in, const string& format) {
	if(format == "table")
		return loadTable(in);
	else {
		errorLog << "Cannot load OTUTable, unsupported format '" << format << "'" << endl;
		in.setstate(std::ios_base::failbit);
		return in;
	}
}

inline std::ostream& OTUTable::save(ostream& out, const string& format) const {
	if(format == "table")
		return saveTable(out);
	else {
		errorLog << "Cannot save OTUTable, unsupported format '" << format << "'" << endl;
		out.setstate(std::ios_base::failbit);
		return out;
	}
}

inline bool OTUTable::removeSample(const string& sampleName) {
	return removeSample(std::find(samples.begin(), samples.end(), sampleName) - samples.begin());
}

inline bool OTUTable::removeOTU(const string& otuID) {
	return removeOTU(std::find(otus.begin(), otus.end(), otuID) - otus.begin());
}

inline OTUTable operator+(const OTUTable& lhs, const OTUTable& rhs) {
	OTUTable otuMerged(lhs); /* make a local copy */
	return otuMerged += rhs;
}


} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

#endif /* SRC_OTUTABLE_H_ */
