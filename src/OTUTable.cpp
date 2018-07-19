/*
 * OTUTable.cpp
 *
 *  Created on: Jul 11, 2017
 *      Author: zhengqi
 */

#include <ctime>
#include <cassert>
#include <algorithm>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/random/random_number_generator.hpp> /* an adapter between randome_number_generator and uniform_random_number_generator */
#include <Eigen/Dense>
#include "StringUtils.h"
#include "HmmUFOtuConst.h"
#include "OTUTable.h"

namespace EGriceLab {
namespace HmmUFOtu {

using namespace std;
using namespace Eigen;

const IOFormat OTUTable::dblTabFmt(FullPrecision, DontAlignCols, "\t", "\n", "", "", "");
const IOFormat OTUTable::fltTabFmt(FullPrecision, DontAlignCols, "\t", "\n", "", "", "");
OTUTable::RNG OTUTable::rng(time(NULL)); /* initiate RNG with random seed */

void OTUTable::clear() {
	samples.clear();
	otus.clear();
	metric.resize(0, 0);
	otu2Taxon.clear();
}

size_t OTUTable::addSample(const string& sampleName) {
	const size_t N = getSampleIndex(sampleName);
	if(N < numSamples()) /* already exists */
		return N;

	samples.push_back(sampleName);
	metric.conservativeResize(Eigen::NoChange, N + 1);
	metric.col(N).setZero();

	return N;
}

void OTUTable::removeSample(size_t j) {
	const size_t N = metric.cols();
	samples.erase(samples.begin() + j);
	MatrixXd oldMetric(metric); /* copy the old values */
	metric.resize(Eigen::NoChange, N - 1);
	for(MatrixXd::Index n = 0, k = 0; n < N; ++n) {
		if(n != j) { /* not the deleted column */
			metric.col(k) = oldMetric.col(n);
			k++;
		}
	}
}

size_t OTUTable::addOTU(const string& otuID, const string& taxon, const RowVectorXd& count) {
	const size_t M = getOTUIndex(otuID);
	if(M < numOTUs()) /* already exists */
		return M;

	otus.push_back(otuID);
	otu2Taxon[otuID] = taxon;
	metric.conservativeResize(M + 1, Eigen::NoChange);
	metric.row(M) = count;

	return M;
}

void OTUTable::removeOTU(size_t i) {
	const size_t M = metric.rows();
	otu2Taxon.erase(otus[i]); /* remove taxon first */
	otus.erase(otus.begin() + i); /* remove the actual OTU */
	MatrixXd oldMetric(metric); /* copy the old values */
	metric.resize(M - 1, Eigen::NoChange);
	for(MatrixXd::Index m = 0, k = 0; m < M; ++m) {
		if(m != i) { /* not the deleted row */
			metric.row(k) = oldMetric.row(m);
			k++;
		}
	}
}

void OTUTable::pruneSamples(size_t min) {
	if(min == 0)
		return;

	const size_t N = numSamples();
	/* remove samples backwards */
	for(size_t j = N; j > 0; --j) {
		if(sumSampleMetric(j - 1) < min)
			removeSample(j - 1);
	}
}

void OTUTable::pruneOTUs(size_t min) {
	const size_t M = numOTUs();
	/* remove samples backwards */
	for(size_t i = M; i > 0; --i) {
		double nRead = sumOTUMetric(i - 1);
		if(min > 0 && nRead < min || min == 0 && nRead == 0)
			removeOTU(i - 1);
	}
}

void OTUTable::normalizeConst(double Z) {
	assert(Z >= 0);
	if(empty() || (metric.array() == 0).all()) /* empty or all zero metric */
		return;
	if(Z == 0)
		Z = metric.colwise().sum().maxCoeff(); /* use max column sum as constant */

	const size_t N = metric.cols();
	RowVectorXd norm = metric.colwise().sum() / Z;
	for(MatrixXd::Index j = 0; j < N; ++j)
		metric.col(j) /= norm(j);
}

istream& OTUTable::loadTable(istream& in) {
	clear(); /* clear old data */
	/* input header */
	string line;
	size_t N = 0;
	while(std::getline(in, line)) {
		if(StringUtils::startsWith(line, "otuID")) { /* header line */
			vector<string> headers;
			boost::split(headers, line, boost::is_any_of("\t"));
			N = headers.size() - 2;
			/* update samples */
			samples.resize(N);
			std::copy(headers.begin() + 1, headers.end() - 1, samples.begin());
			metric.resize(0, N);
		}
		else { /* value line */
			string otuID, taxon;
			istringstream lineIn(line);
			std::getline(lineIn, otuID, '\t');
			RowVectorXd counts(N);
			for(size_t j = 0; j < N; ++j)
				lineIn >> counts(j);
			lineIn.ignore(1, '\t');
			std::getline(lineIn, taxon);
			addOTU(otuID, taxon, counts); /* add a new OTU */
		}
	}

	return in;
}

ostream& OTUTable::saveTable(ostream& out) const {
	/* output header */
	out << "otuID\t" << boost::join(samples, "\t") << "\ttaxonomy" << endl;

	/* output each OTU */
	const size_t M = numOTUs();
	for(size_t i = 0; i < M; ++i)
		out << otus[i] << "\t" << metric.row(i).format(fltTabFmt) << "\t" << otu2Taxon.at(otus[i]) << endl;

	return out;
}

void OTUTable::subsetUniform(size_t min) {
	for(int j = 0; j < numSamples(); ++j) {
		double sampleTotal = sumSampleMetric(j);
		if(sampleTotal <= min) /* not enough reads to subset */
			continue;

		/* generate an sampling index with length M */
		std::vector<bool> otuIdx(static_cast<size_t> (sampleTotal), false); /* use the efficient std::vector<bool>, default all false */
		fill_n(otuIdx.begin(), min, true);
		boost::random_number_generator<RNG, size_t> gen(rng);
		boost::random_shuffle(otuIdx, gen);
		/* subset reads in OTUs without replacement using the random index */
		for(size_t i = 0, k = 0; i < numOTUs(); ++i) { /* k is the start index of current OTU */
			size_t N = static_cast<size_t> (metric(i, j));
			metric(i, j) = std::count(otuIdx.begin() + k, otuIdx.begin() + k + N, true);
			assert(metric(i, j) <= N);
			k += N;
		}
	}
}

void OTUTable::subsetMultinom(size_t min) {
	const size_t M = numOTUs();
	double *otuPr = new double[M]; /* raw read sample probabilities */
	Map<VectorXd> otuPrMap(otuPr, M); /* use a map to access indirectly */
	otuPrMap.setOnes(); /* use all equal probs by default */
	ReadDistrib rdist(otuPr, otuPr + M); /* construct the discrete distribution */

	for(int j = 0; j < numSamples(); ++j) {
		double sampleTotal = sumSampleMetric(j);
		if(sampleTotal <= min) /* not enough reads to subset */
			continue;

		/** reset rdist probabilities according to current counts */
		otuPrMap = metric.col(j);
		rdist.param(ReadParam(otuPr, otuPr + M));
		/* sample min reads */
		VectorXd sampled = VectorXd::Zero(M);
		for(size_t m = 0; m < min; ++m)
			sampled(rdist(rng))++;
		metric.col(j) = sampled;
	}
	delete[] otuPr;
}

OTUTable& OTUTable::operator+=(const OTUTable& other) {
	if(empty()) {
		*this = other;
		return *this;
	}
	if(other.empty())
		return *this;

	/* add non-existing samples */
	for(size_t j = 0; j < other.numSamples(); ++j)
		addSample(other.getSample(j));

	/* add non-existing OTUs */
	for(size_t i = 0; i < other.numOTUs(); ++i) {
		string otuID = other.getOTU(i);
		addOTU(otuID, other.getTaxon(otuID));
	}

	/** merge counts */
	for(size_t i = 0; i < other.numOTUs(); ++i) {
		string otuID = other.getOTU(i);
		int i0 = getOTUIndex(otuID);
		for(size_t j = 0; j < other.numSamples(); ++j) {
			int j0 = getSampleIndex(other.getSample(j));
			metric(i0, j0) += other.numMetric(i, j);
		}
	}

	return *this;
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */
