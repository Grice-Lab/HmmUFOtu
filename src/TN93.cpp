/*
 * TN93.cpp
 *
 *  Created on: Mar 7, 2017
 *      Author: zhengqi
 */

#include "TN93.h"
#include "ProgLog.h"

namespace EGriceLab {
using namespace std;
using namespace Eigen;

const string TN93::name = "TN93";

istream& TN93::read(istream& in) {
	string line, tag, value;
	while(in >> tag) {
		if(tag[0] == '#') { /* comment or header */
			std::getline(in, line); /* ignore the entire line */
			continue;
		}
		if(tag == "Type:") {
			in >> value; // read in model type
			if(value != modelType()) {
				errorLog << "Unmatched Model Type!" << endl;
				errorLog << "Trying to read in a " << value << " model into a " << modelType() << " object" << endl;
				in.setstate(ios_base::badbit);
				return in;
			}
		}
		else if(tag == "pi:") {
			for(Vector4d::Index i = 0; i != pi.rows(); ++i)
				in >> pi(i);
		}
		else if(tag == "kr:")
			in >> kr;
		else if(tag == "ky:")
			in >> ky;
		else if(tag == "beta:") {
			in >> beta;
			std::getline(in, line); /* ignore the entire line */
			break;
		}
		else {
			errorLog << "Un-recognized line found in TN93 Model input: tag: " << tag << endl << line << endl;
			in.setstate(ios_base::badbit);
			return in;
		}
	}

	return in;
}

ostream& TN93::write(ostream& out) const {
	out << "# DNA Substitution Model" << endl;
	out << "Type: " << modelType() << endl;
	out << "pi: " << pi.transpose().format(FULL_FORMAT) << endl;
	out << "kr: " << kr << " ky: " << ky << " beta: " << beta << endl;

	return out;
}

void TN93::trainParams(const vector<Matrix4d>& Pv, const Vector4d& f) {
	/* estimate pi using mean f */
	pi = f / f.sum();
	/* estimate kr and ky */
	for(vector<Matrix4d>::const_iterator P = Pv.begin(); P != Pv.end(); ++P) {
		double Tr = (*P)(A, G) + (*P)(G, A);
		double Ty = (*P)(C, T) + (*P)(T, C);
		double Tv = (*P)(A, C) + (*P)(A, T) + (*P)(C, A) + (*P)(C, G) + (*P)(G, C) + (*P)(G, T) + (*P)(T, A) + (*P)(T, G);
		kr += Tr / Tv;
		ky += Ty / Tv;
	}
	kr /= Pv.size();
	ky /= Pv.size();
	/* estimate beta */
	setBeta();
}

} /* namespace EGriceLab */
