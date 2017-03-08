/*
 * JC69.cpp
 *
 *  Created on: Mar 7, 2017
 *      Author: zhengqi
 */

#include "JC69.h"
#include "ProgLog.h"

namespace EGriceLab {
using namespace std;
using namespace Eigen;

const string JC69::name = "JC69";
const Vector4d JC69::pi = Vector4d::Constant(1.0 / 4);

istream& JC69::read(istream& in) {
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
			std::getline(in, line); /* ignore the reset of line */
			break;
		}
		else {
			errorLog << "Un-recognized line found in JC69 Model input: tag: " << tag << endl << line << endl;
			in.setstate(ios_base::badbit);
			return in;
		}
	}

	return in;
}

ostream& JC69::write(ostream& out) const {
	out << "# DNA Substitution Model" << endl;
	out << "Type: " << modelType() << endl;

	return out;
}

} /* namespace EGriceLab */
