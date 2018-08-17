/*
 * GFF.cpp
 *
 *  Created on: Aug 7, 2018
 *      Author: zhengqi
 */

#include <math.h> // requires C99
#include <limits>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "StringUtils.h"
#include "GFF.h"

namespace EGriceLab {
namespace UCSC {

using namespace std;
const string GFF::GFF_SUFFIX = ".gff";
const string GFF::GTF_SUFFIX = ".gtf";
const string GFF::GFF3_SUFFIX = ".gff3";

const double GFF::INVALID_SCORE = NAN;
const string GFF::INVALID_TOKEN = ".";

void GFF::setAttr(const string& name, const string& val) {
	if(attrValues.count(name) == 0) /* not exists yet */
		attrNames.push_back(name);
	attrValues[name] = val;
}

istream& GFF::load(istream& in) {
	StringUtils::loadString(seqname, in);
	StringUtils::loadString(source, in);
	StringUtils::loadString(type, in);
	in.read((char *) &start, sizeof(long));
	in.read((char *) &end, sizeof(long));
	in.read((char *) &score, sizeof(double));
	in.read(&strand, 1);
	in.read((char *) &frame, sizeof(int));
	size_t nAttr = 0;
	string name, val;
	in.read((char *) &nAttr, sizeof(size_t));
	for(size_t i = 0; i < nAttr; ++i) {
		StringUtils::loadString(name, in);
		StringUtils::loadString(val, in);
		setAttr(name, val);
	}
	return in;
}

ostream& GFF::save(ostream& out) const {
	StringUtils::saveString(seqname, out);
	StringUtils::saveString(source, out);
	StringUtils::saveString(type, out);
	out.write((const char *) &start, sizeof(long));
	out.write((const char *) &end, sizeof(long));
	out.write((const char *) &score, sizeof(double));
	out.write(&strand, 1);
	out.write((const char *) &frame, sizeof(int));
	size_t nAttr = numAttrs();
	out.write((const char *) &nAttr, sizeof(size_t));
	for(attr_map::const_iterator pair = attrValues.begin(); pair != attrValues.end(); ++pair) {
		StringUtils::saveString(pair->first, out);
		StringUtils::saveString(pair->second, out);
	}
	return out;
}

istream& operator>>(istream& in, GFF& record) {
	string token;
	std::getline(in, record.seqname, GFF::SEP);
	std::getline(in, record.source, GFF::SEP);
	std::getline(in, record.type, GFF::SEP);
	in >> record.start >> record.end;
	/* safe-guard invalid fields */
	in >> token;
	record.score = token != GFF::INVALID_TOKEN ? boost::lexical_cast<double>(token) : GFF::INVALID_SCORE;

	in >> record.strand;

	in >> token;
	record.frame = token != GFF::INVALID_TOKEN ? boost::lexical_cast<int>(token) : GFF::INVALID_FRAME;

	in.ignore(std::numeric_limits<streamsize>::max(), GFF::SEP);
	string attrStr;
	std::getline(in, attrStr);
	record.readAttributes(attrStr);
	return in;
}

ostream& operator<<(ostream& out, const GFF& record) {
	out << record.seqname << GFF::SEP << record.source << GFF::SEP << record.type << GFF::SEP
		<< record.start << GFF::SEP << record.end << GFF::SEP;
	if(::isnan(record.score))
		out << GFF::INVALID_FLAG;
	else
		out << record.score;
	out << GFF::SEP << record.strand << GFF::SEP;
	if(record.frame == GFF::INVALID_FRAME)
		out << GFF::INVALID_FLAG;
	else
		out << record.frame;
	out << GFF::SEP << record.writeAttributes();

	return out;
}

GFF::Version GFF::guessVersion(const string& fn) {
	if(StringUtils::endsWith(fn, GTF_SUFFIX))
		return GTF;
	else if(StringUtils::endsWith(fn, GFF_SUFFIX) || StringUtils::endsWith(fn, GFF3_SUFFIX))
		return GFF3;
	else return UNK;
}

void GFF::readGTFAttributes(const string& attrStr) {
	vector<string> attrs;
	boost::split(attrs, attrStr, boost::is_any_of(" \";"), boost::token_compress_on);
	for(vector<string>::size_type i = 0; i < attrs.size(); i += 2)
		if(!attrs[i].empty())
			setAttr(attrs[i], attrs[i+1]);
}

string GFF::writeGTFAttributes() const {
	string attrStr;
	const vector<string>& attrNames = getAttrNames();
	for(vector<string>::const_iterator name = attrNames.begin(); name != attrNames.end(); ++name) {
		if(!attrStr.empty()) /* non-first */
			attrStr += "; ";
		attrStr += (*name) + " \"" + getAttr(*name) + "\"";
	}
	return attrStr;
}

void GFF::readGFF3Attributes(const string& attrStr) {
	vector<string> attrs;
	boost::split(attrs, attrStr, boost::is_any_of("=;"), boost::token_compress_on);
	for(vector<string>::size_type i = 0; i < attrs.size(); i += 2)
		if(!attrs[i].empty())
			setAttr(attrs[i], attrs[i+1]);
}

string GFF::writeGFF3Attributes() const {
	string attrStr;
	const vector<string>& attrNames = getAttrNames();
	for(vector<string>::const_iterator name = attrNames.begin(); name != attrNames.end(); ++name) {
		if(!attrStr.empty()) /* non-first */
			attrStr += ";";
		attrStr += (*name) + "=" + getAttr(*name);
	}
	return attrStr;
}

} /* namespace UCSC */
} /* namespace EGriceLab */

