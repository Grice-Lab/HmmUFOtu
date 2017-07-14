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
#include "ProgLog.h"
#include "OTUObserved.h"

namespace EGriceLab {

using std::string;
using std::vector;
using std::map;
using std::invalid_argument;
using std::find;
using std::istream;
using std::ostream;
using Eigen::MatrixXd;
using Eigen::RowVectorXd;

class OTUTable {
	/* typedefs */
	typedef map<string, string> otuMap;

public:
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
		return otuMetric.size() == 0;
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
	 * get the index of a given sample, or -1 if not found
	 */
	size_t getSampleIndex(const string& sampleName) const {
		vector<string>::const_iterator result = find(samples.begin(), samples.end(), sampleName);
		return result != samples.end() ? result - samples.begin() : -1;
	}

	/**
	 * get the index of a given OTU, or -1 if not found
	 */
	size_t getOTUIndex(const string& otuID) const {
		vector<string>::const_iterator result = find(otus.begin(), otus.end(), otuID);
		return result != otus.end() ? result - otus.begin() : -1;
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

	/**
	 * get total count of a given sample across all OTUs, or 0 if not found
	 */
	double sampleTotal(const string& sampleName) const;

	/**
	 * get total count of an OTU across all samples, or 0 if not found
	 */
	double otuTotal(const string& otuTotal) const;

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
	 * normalize the metric using given normalization constant
	 * @param Z  normalization constant
	 * @return  the modifled OTUTable
	 */
	void normalize(double Z = 0);

	/**
	 * subset this OTU table to a minimum read count using given method
	 * samples that have less than min reads will be removed
	 * samples that have more than min reads will be subsampled
	 * @param min  min read requirement
	 * @param method  sampleing method
	 * @throw  invalid_argument if the sampling method is not supported
	 */
	void subset(double min, const string& method) {
		if(method == "uniform")
			subsetUniform(min);
		else if(method == "multinom")
			subsetMultinom(min);
		else
			throw invalid_argument("Unsupported subsetting method '" + method + "'");
	}

	/**
	 * subset this OTU table to a minimum read count using uniform sampling
	 */
	void subsetUniform(double min);

	/**
	 * subset this OTU table to a minimum read count using Multinomial sampling
	 */
	void subsetMultinom(double min);

	/**
	 * load a table from an input stream
	 */
	istream& load(istream& in, const string& format = "table");

	/**
	 * save this table to an output stream
	 */
	ostream& save(ostream& out, const string& format = "table") const;

	/**
	 * load a table from an input stream in OTU table format
	 */
	istream& loadTable(istream& in);

	/**
	 * save this table to an output stream in OTU table format
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

	/** non-member functions */
	/** formatted input, equivalent to the default behavior of load() */
	friend istream& operator>>(istream& in, OTUTable& otu);

	/** formatted output, equivalent to the default behavior of save() */
	friend ostream& operator<<(ostream& out, const OTUTable& otu);

private:
	/** member fields */
	vector<string> samples; /* 0..N sample names */
	vector<string> otus;    /* 0..M OTUs */
	otuMap otu2Taxon;
	MatrixXd otuMetric; /* M * N matrix of OTU (relative) abundance metric */

	/** static fields */
	static const Eigen::IOFormat dblTabFmt;
	static const Eigen::IOFormat fltTabFmt;
};

inline void OTUTable::clear() {
	samples.clear();
	otus.clear();
	otuMetric.resize(0, 0);
	otu2Taxon.clear();
}

inline double OTUTable::sampleTotal(const string& sampleName) const {
	size_t j = getSampleIndex(sampleName);
	return j != -1 ? otuMetric.col(j).sum() : 0;
}

inline double OTUTable::otuTotal(const string& otuID) const {
	size_t i = getOTUIndex(otuID);
	return i != -1 ? otuMetric.row(i).sum() : 0;
}

inline std::istream& operator>>(std::istream& in, OTUTable& otu) {
	return otu.load(in);
}

inline std::ostream& operator<<(std::ostream& out, OTUTable& otu) {
	return otu.save(out);
}

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

} /* namespace EGriceLab */

#endif /* SRC_OTUTABLE_H_ */
