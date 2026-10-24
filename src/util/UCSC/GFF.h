/*
 * GFF.h
 *
 *  Created on: Aug 7, 2018
 *      Author: zhengqi
 */

#ifndef GFF_H_
#define GFF_H_
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <cstdint>

namespace EGriceLab {
namespace UCSC {

using std::string;
using std::vector;
using std::map;
using std::istream;
using std::ostream;
using std::cerr;
using std::endl;

/**
 * an abstract class represent a UCSC GFF (GTF or GFF3) record
 */
class GFF {
public:
	typedef vector<string> attr_list;
	typedef map<string, string> attr_map;

	enum Version { UNK = 1, GTF, GFF3 };

	/* constructors */
	/** default constructor */
	GFF() = default;

	/** construct a GFF with given attr_map */
	GFF(const string& seqname, const string& source, const string& type,
			long start, long end, double score, char strand, int frame, const attr_map& attrs)
	: seqname(seqname), source(source), type(type),
	  start(start), end(end), score(score), strand(strand), frame(frame), attrValues(attrs) {
		attrNames.reserve(attrValues.size());
		for(const attr_map::value_type& entry : attrValues) // use natual order
			attrNames.push_back(entry.first);
	}

	/** construct a GFF record with given info */
	GFF(GFF::Version ver, const string& seqname, const string& source, const string& type,
			long start, long end, double score, char strand, int frame, const string& attrStr = "")
	: seqname(seqname), source(source), type(type),
	  start(start), end(end), score(score), strand(strand), frame(frame) {
		readAttributes(attrStr, ver);
	}

	const vector<string>& getAttrNames() const {
		return attrNames;
	}

	const attr_map& getAttrValues() const {
		return attrValues;
	}

	long getEnd() const {
		return end;
	}

	void setEnd(long end) {
		this->end = end;
	}

	int getFrame() const {
		return frame;
	}

	void setFrame(int frame) {
		this->frame = frame;
	}

	double getScore() const {
		return score;
	}

	void setScore(double score) {
		this->score = score;
	}

	const string& getSeqname() const {
		return seqname;
	}

	void setSeqname(const string& seqname) {
		this->seqname = seqname;
	}

	const string& getSource() const {
		return source;
	}

	void setSource(const string& source) {
		this->source = source;
	}

	long getStart() const {
		return start;
	}

	void setStart(long start) {
		this->start = start;
	}

	char getStrand() const {
		return strand;
	}

	void setStrand(char strand) {
		this->strand = strand;
	}

	const string& getType() const {
		return type;
	}

	void setType(const string& type) {
		this->type = type;
	}

	/* member methods */
	/** test whether this GFF is valid */
	bool isValid() const {
		return start > 0 && start <= end && attrNames.size() == attrValues.size();
	}

	/** get length of this GFF record */
	long length() const {
		return end - start + 1;
	}

	/** get number of attrs */
	size_t numAttrs() const {
		return attrNames.size();
	}

	/** clear all content and reset this GFF to defaul state */
	void clear() {
		seqname.clear();
		source.clear();
		type.clear();
		start = 0;
		end = 0;
		score = INVALID_SCORE;
		strand = DEFAULT_STRAND;
		frame = INVALID_FRAME;
		attrNames.clear();
		attrValues.clear();
	}

	/** get the attr value given its name */
	string getAttr(const string& name) const {
		return attrValues.at(name);
	}

	/** set the attr value with given name, create if not already exists */
	void setAttr(const string& name, const string& val);

	/** test whether this attr name exists */
	bool hasAttr(const string& name) const {
		return attrValues.count(name) > 0;
	}

	/** shift this GFF record of given offset */
	GFF& shift(long offset) {
		start += offset;
		end += offset;
		return *this;
	}

	/** load object from binary input */
	istream& load(istream& in);

	/** save object to binary output */
	ostream& save(ostream& out) const;

	/** formated write to text output */
	ostream& write(ostream& out, GFF::Version ver = GFF::GFF3) const;

	/** formated read from text input */
	istream& read(istream& in, GFF::Version ver = GFF::GFF3);

	/** read attrs from formated text */
	void readAttributes(const string& attrStr, GFF::Version ver);

	/** write attrs to formated text */
	string writeAttributes(GFF::Version ver) const;

	/** read attrs from formated text in GTF format */
	void readGTFAttributes(const string& attrStr);

	/** write attrs to formated text in GTF format */
	string writeGTFAttributes() const;

	/** read attrs from formated text in GFF3 format */
	void readGFF3Attributes(const string& attrStr);

	/** write attrs to formated text in GFF3 format */
	string writeGFF3Attributes() const;

	/* non-member operators */
	/* formatted read */
	friend istream& operator>>(istream& in, GFF& record);

	/* formatted write */
	friend ostream& operator<<(ostream& out, const GFF& record);

	/** relational operators */
	friend bool operator==(const GFF& lhs, const GFF& rhs);
private:
	/* member fields */
	string seqname;
	string source;
	string type;
	int64_t start = 0; /* 1-based */
	int64_t end = 0;   /* 1-based */
	double score = INVALID_SCORE;
	char strand = DEFAULT_STRAND;
	int frame = INVALID_FRAME;
	vector<string> attrNames; /* attr names in original order */
	attr_map attrValues; /* attribute name->value map */

public:
	/* class constants */
	static const char COMMENT_CHAR = '#';
	static const char SEP = '\t';
	static const char GTF_QUOTE = '"';
	static const char GTF_ATTR_SEP = ';';
	static const char GFF3_ATTR_SEP = ';';
	static const char GTF_TAG_SEP = ' ';
	static const char GFF3_TAG_SEP = '=';

	static const string GFF_SUFFIX;
	static const string GTF_SUFFIX;
	static const string GFF3_SUFFIX;
	static const double INVALID_SCORE;
	static const char DEFAULT_STRAND = '.';
	static const int INVALID_FRAME = -1;
	static const char INVALID_FLAG = '.';
	static const string INVALID_TOKEN;

	/* static methods */
	static Version guessVersion(const string& fn);
};

inline void GFF::readAttributes(const string& attrStr, GFF::Version ver) {
	switch(ver) {
	case GTF:
		return readGTFAttributes(attrStr);
	case GFF3:
		return readGFF3Attributes(attrStr);
	}
}

inline string GFF::writeAttributes(GFF::Version ver) const {
	switch(ver) {
	case GTF:
		return writeGTFAttributes();
	case GFF3:
		return writeGFF3Attributes();
	default:
		return "";
	}
}

inline istream& operator>>(istream& in, GFF& record) {
	return record.read(in);
}

inline ostream& operator<<(ostream& out, const GFF& record) {
	return record.write(out);
}

inline bool operator==(const GFF& lhs, const GFF& rhs) {
	return lhs.seqname == rhs.seqname && lhs.source == rhs.source && lhs.type == rhs.type &&
			lhs.start == rhs.start && lhs.end == rhs.end &&
			(std::isnan(lhs.score) && std::isnan(rhs.score) || lhs.score == rhs.score) &&
			lhs.strand == rhs.strand && lhs.frame == rhs.frame;
}

inline bool operator!=(const GFF& lhs, const GFF& rhs) {
	return !(lhs == rhs);
}

} /* namespace UCSC */
} /* namespace EGriceLab */

#endif /* GFF_H_ */
