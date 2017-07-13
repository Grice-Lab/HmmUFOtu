/*
 * OTUTable.cpp
 *
 *  Created on: Jul 11, 2017
 *      Author: zhengqi
 */

#include <cassert>
#include <boost/lexical_cast.hpp>
#include "HmmUFOtuConst.h"
#include "OTUTable.h"

namespace EGriceLab {
using namespace std;
using namespace Eigen;

const IOFormat OTUTable::dblTabFmt(FullPrecision, DontAlignCols, "\t", "\n", "", "", "");
const IOFormat OTUTable::fltTabFmt(FullPrecision, DontAlignCols, "\t", "\n", "", "", "");

bool OTUTable::addSample(const string& sampleName) {
	if(hasSample(sampleName))
		return false;

	const size_t N = otuMetric.cols();
	samples.push_back(sampleName);
	otuMetric.conservativeResize(Eigen::NoChange, N + 1);
	otuMetric.col(N).setZero();

	return true;
}

bool OTUTable::removeSample(const string& sampleName) {
	vector<string>::iterator result = find(samples.begin(), samples.end(), sampleName);
	if(result == samples.end())
		return false;

	const size_t N = otuMetric.cols();
	samples.erase(result);
	MatrixXd oldMetric(otuMetric); /* copy the old values */
	otuMetric.resize(Eigen::NoChange, N - 1);
	for(MatrixXd::Index j = 0, k = 0; j < N; ++j) {
		if(j != result - samples.begin()) { /* not the deleted column */
			otuMetric.col(k) = oldMetric.col(j);
			k++;
		}
	}

	return true;
}

bool OTUTable::addOTU(const string& otuID, const string& taxon, const RowVectorXd& count) {
	if(hasOTU(otuID))
		return false;

	const size_t M = otuMetric.rows();
	otus.push_back(otuID);
	otu2Taxon[otuID] = taxon;
	otuMetric.conservativeResize(M + 1, Eigen::NoChange);
	otuMetric.row(M) = count;

	return true;
}

bool OTUTable::removeOTU(const string& otuID) {
	vector<string>::iterator result = find(otus.begin(), otus.end(), otuID);
	if(result == otus.end())
		return false;

	const size_t M = otuMetric.rows();
	otus.erase(result);
	otu2Taxon.erase(otuID);
	MatrixXd oldMetric(otuMetric); /* copy the old values */
	otuMetric.resize(M - 1, Eigen::NoChange);
	for(MatrixXd::Index i = 0, k = 0; i < M; ++i) {
		if(i != result - otus.begin()) { /* not the deleted row */
			otuMetric.row(k) = oldMetric.row(k);
			k++;
		}
	}

	return true;
}

void OTUTable::normalize(double Z) {
	assert(Z > 0);
	const size_t N = otuMetric.cols();
	RowVectorXd norm = otuMetric.colwise().sum() / Z;
	for(MatrixXd::Index j = 0; j < N; ++j)
		otuMetric.col(j) /= norm(j);
}

istream& OTUTable::loadTable(istream& in) {
	clear(); /* clear old data */
	/* input header */
	string line;
	bool isHeader = true;
	size_t M = 0;
	size_t N = 0;
	while(std::getline(in, line)) {
		if(line.front() == '#') { /* header line */
			char pname[1024];
			char pver[1024];
			sscanf("# OTU table generated by %s %s", pname, pver);
			if(pname != progName) {
				errorLog << "Not a OTU Table file" << endl;
				in.setstate(ios_base::failbit);
				return in;
			}
			if(cmpVersion(progVersion, pver) < 0) {
				errorLog << "You are trying using an older version " << (progName + progVersion) <<
						" to read a newer OTU Table file that was build by " << (string(pname) + string(pver)) << endl;
				in.setstate(ios_base::failbit);
				return in;
			}
		}
		vector<string> fields;
		boost::split(fields, line, boost::is_any_of("\t"));
		if(isHeader) {
			N = fields.size() - 2;
			/* update samples */
			samples.resize(N);
			std::copy(fields.begin() + 1, fields.end() - 1, samples.begin());
			otuMetric.resize(0, N);
			isHeader = false;
		}
		else { /* value line */
			if(fields.size() != N + 2)
				continue;
			bool isNew = addOTU(fields.front(), fields.back()); /* add a new OTU */
			if(!isNew)
				continue;
			for(size_t j = 0; j < fields.size(); ++j)
				otuMetric(M, j) = boost::lexical_cast<double> (fields[j + 1]);
			M++;
		}
	}

	return in;
}

ostream& OTUTable::saveTable(ostream& out) const {
	/* output header */
	out << "# OTU table generated by " << progName << " " << progVersion << endl;
	out << "otuID\t" << boost::join(samples, "\t") << "\ttaxonomy" << endl;

	/* output each OTU */
	const size_t M = numOTUs();
	for(size_t i = 0; i < M; ++i)
		out << otus[i] << "\t" << otuMetric.row(i).format(fltTabFmt) << "\t" << otu2Taxon.at(otus[i]) << endl;

	return out;
}

} /* namespace EGriceLab */
