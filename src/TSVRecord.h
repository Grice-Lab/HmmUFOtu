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
 * TSVRecord.h
 *  Class designed for Tab-Separated-Values format plain-text tables
 *  Created on: Jul 12, 2017
 *      Author: zhengqi
 */

#ifndef SRC_TSVRECORD_H_
#define SRC_TSVRECORD_H_
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

namespace EGriceLab {
using std::string;
using std::vector;
using std::map;
using std::istream;
using std::ostream;

class TSVRecord {
public:
	/** nested types and enums */
	struct TSVHeader {
		/** constructors */
		/** default constructor */
		TSVHeader() {  }

		/** construct a header with a list of names */
		explicit TSVHeader(const vector<string>& headerNames): names(headerNames) {
			setHeaderIndex();
		}

		/** construct a header with a header line */
		explicit TSVHeader(const string& line, const string& sep = DEFAULT_SEP) {
			boost::split(names, line, boost::is_any_of(sep));
			setHeaderIndex();
		}

		/** test weither this header is empty */
		bool empty() const {
			return names.empty();
		}

		/** get number of header fields */
		size_t numHeader() const {
			return names.size();
		}

		/** test wether a header name exists */
		bool hasHeader(const string& name) const {
			return index.count(name) > 0;
		}

		/** get header name by index */
		const string& getHeaderName(size_t i) const {
			return names[i];
		}

		/** get header index by name */
		size_t getHeaderIndex(const string& name) const {
			return index.at(name);
		}

		/** add a new header field */
		void addHeader(const string& name);

		/** remove a header */
		void removeHeader(const string& name);

		/** convert header to a line */
		string toString(const string& sep = DEFAULT_SEP, char quote = DEFAULT_QUOTE) const;

		/** non-member functions */
		friend istream& operator>>(istream& in, TSVHeader& header);

		friend ostream& operator<<(ostream& out, const TSVHeader& header);

	private:
		/** internal method */
		void setHeaderIndex();

		/** member fields */
		vector<string> names;
		map<string, size_t> index;

	};

	typedef boost::shared_ptr<TSVHeader> TSVHeaderPtr;

	/** constructors */
	/** default constructor */
	TSVRecord() {  }

	/** construct a record with given fields */
	explicit TSVRecord(const vector<string>& fields, const TSVHeaderPtr& headerIdx = NULL) :
			fields(fields), header(headerIdx)
	{  }

	/** construct a record with a input line */
	explicit TSVRecord(const string& line, const TSVHeaderPtr& headerIdx = NULL,
			string sep = DEFAULT_SEP, char quote = DEFAULT_QUOTE) :
					header(headerIdx)
	{
		parse(line, sep, quote);
	}

	/** destructor, do nothing */
	virtual ~TSVRecord() {  }

	/** member methods*/
	/** getters and setters */
	const vector<string>& getFields() const {
		return fields;
	}

	/** get number of fields */
	size_t numFields() const {
		return fields.size();
	}

	/** get field by index */
	const string& getField(size_t i) const {
		return fields[i];
	}

	/** set field at index */
	void setField(size_t i, const string& val) {
		fields[i] = val;
	}

	/** test whether this TSVRecord has an associated header */
	bool hasHeader() const {
		return header != NULL && !header->empty();
	}

	/** get field by name */
	const string& getFieldByName(const string& name) const {
		return getField(header->getHeaderIndex(name));
	}

	/** set field by name */
	void setFieldByName(const string& name, const string& val) {
		setField(header->getHeaderIndex(name), val);
	}

	/** convert record to a line */
	string toString(const string& sep = DEFAULT_SEP, char quote = DEFAULT_QUOTE) const;

	/** non-member functions */
	friend istream& operator>>(istream& in, TSVRecord& record);

	friend ostream& operator<<(ostream& out, const TSVRecord& record);

private:
	/** member fields */
	vector<string> fields;
	TSVHeaderPtr header;

	/** static fields */
public:
	static const string DEFAULT_SEP;
	static const char DEFAULT_QUOTE = '\0';

	/** internal member methods */
private:
	void parse(const string& line, const string& sep, char quote);

};

inline istream& operator>>(istream& in, TSVRecord& record) {
	string line;
	std::getline(in, line);
	record.parse(line, TSVRecord::DEFAULT_SEP, TSVRecord::DEFAULT_QUOTE);
	return in;
}

inline ostream& operator<<(ostream& out, const TSVRecord& record) {
	out << record.toString();
	return out;
}

inline istream& operator>>(istream& in, TSVRecord::TSVHeader& header) {
	string line;
	std::getline(in, line);
	boost::split(header.names, line, boost::is_any_of(TSVRecord::DEFAULT_SEP));
	header.setHeaderIndex();
	return in;
}

inline ostream& operator<<(ostream& out, const TSVRecord::TSVHeader& header) {
	out << header.toString();
	return out;
}

} /* namespace EGriceLab */

#endif /* SRC_TSVRECORD_H_ */
