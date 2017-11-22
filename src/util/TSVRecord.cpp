/*
 * TSVRecord.cpp
 *
 *  Created on: Jul 12, 2017
 *      Author: zhengqi
 */

#include <algorithm>
#include <cctype>
#include "TSVRecord.h"
#include "StringUtils.h"

namespace EGriceLab {

const string TSVRecord::DEFAULT_SEP = "\t";

void TSVRecord::parse(const string& line, const string& sep, char quote) {
	boost::split(fields, line, boost::is_any_of(sep));
	if(quote != '\0') { /* quote requested */
		for(vector<string>::iterator val = fields.begin(); val != fields.end(); ++val)
			*val = StringUtils::stripQuotes(*val, quote);
	}
}

string TSVRecord::toString(const string& sep, char quote) const {
	if(!::isprint(quote)) /* non-printable quote charaster */
		return boost::join(fields, sep);

	vector<string> fieldsQuoted;
	fieldsQuoted.reserve(numFields());
	for(vector<string>::const_iterator val = fields.begin(); val != fields.end(); ++val)
		fieldsQuoted.push_back(quote + *val + quote);

	return boost::join(fieldsQuoted, sep);
}

/** add a new header field */
void TSVRecord::TSVHeader::addHeader(const string& name) {
	if(hasHeader(name))
		return;
	names.push_back(name);
	size_t n = index.size();
	index[name] = n;
}

/** remove a header */
void TSVRecord::TSVHeader::removeHeader(const string& name) {
	if(!hasHeader(name))
		return;
	index.erase(name);
	names.erase(std::remove(names.begin(), names.end(), name), names.end());
}

void TSVRecord::TSVHeader::setHeaderIndex() {
	for(vector<string>::size_type i = 0; i < names.size(); ++i)
		index[names[i]] = i;
}

} /* namespace EGriceLab */
