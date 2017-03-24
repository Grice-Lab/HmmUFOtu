/*
 * DigitalSeq.h
 *
 *  Created on: May 5, 2015
 *      Author: zhengqi
 */

#ifndef DIGITALSEQ_H_
#define DIGITALSEQ_H_

#include <string>
#include <bits/basic_string.h>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include "DegenAlphabet.h"
#include "PrimarySeq.h"

namespace EGriceLab {
using std::string;
using std::istream;
using std::ostream;

/**
 * A Digital representation of a sequence, so characters in the Alphabet will be represented as 0,1,...,N-1
 * Note that A DigitalSeq is always case-insensitive, and everything is stored in upper case internally
 * Also note that a DigitalSeq's life is dependent on the life-span of the underlying alphabet; no automatic memory
 * management is carried out
 */
class DigitalSeq: public std::basic_string<int8_t> {
public:
	/* constructors */
	/** default constructor, do nothing */
	DigitalSeq() : abc(NULL) { }

	/** Construct a DigitalSeq with given alphabet, name and string, invalid chars ignored
	 * @param dgAbc  A DegenAlphabet
	 * @param name  name of this ds
	 * @param str  string of this ds
	 */
	explicit DigitalSeq(const DegenAlphabet* abc, const string& name = "", const string& str = "");

	/**
	 Construct a DigitalSeq from a PrimrarySeq
	 */
	explicit DigitalSeq(const PrimarySeq& seq);

	/* virtual destructor */
	virtual ~DigitalSeq() { }

	/* Getters and Setters */
	const DegenAlphabet* getAbc() const {
		return abc;
	}

	void setAbc(const DegenAlphabet* abc) {
		this->abc = abc;
	}

	const string& getName() const {
		return name;
	}

	void setName(const string& name) {
		this->name = name;
	}

	/* utility member methods */

	/** Return the non-gap length of this seq */
	DigitalSeq::size_type nonGapLength() const {
		return length() - std::count(begin(), end(), DegenAlphabet::GAP_BASE);
	}

	/**
	 * Return the string representation of this DigitalSeq
	 */
	string toString() const;

	/**
	 * Generate the reverse complement copy of this DigitalSeq
	 * return a new copy in reverse complement version
	 * or throw an exception if the Alphabet desn't not support complement
	 */
	DigitalSeq revcom() const;

	/**
	 * Get a joint digit representation of current DigitalSeq, such as 0,1,2,1,3
	 * @param sep  the separator character or string, default is ","
	 */
	string join(const string& sep = ",");

	/**
	 * Get the decoded character at given position
	 * @param i  position within this object
	 * @return  the decoded character according to the underlying DegenAlphabet
	 */
	char decodeAt(DigitalSeq::size_type i) const {
		return abc->decode((*this)[i]);
	}

	/**
	 * Alias of toString method
	 */
	string decode() const {
		return toString();
	}

	/**
	 * test whether the encoded value position i is a symbol
	 * param i  position within this object
	 * @return  true if ith code is a symbol
	 */
	bool isSymbol(DigitalSeq::size_type i) const {
		return operator[](i) >= 0;
	}

	/**
	 * test whether the encoded value position i is a gap
	 * param i  position within this object
	 * @return  true if ith code is a gap
	 */
	bool isGap(DigitalSeq::size_type i) const {
		return operator[](i) == DegenAlphabet::GAP_BASE;
	}

	/**
	 * Append a new string to this DigitalSeq
	 * return the modified *this
	 */
	DigitalSeq& append(const string& str);

	/**
	 * Re-introduce all base class append methods
	 */
	using basic_string<int8_t>::append;

	/**
	 * load data from input
	 */
	istream& load(istream& in);

	/**
	 * save this seq to output in binary format
	 */
	ostream& save(ostream& out) const;

	/* static methods */

private:
	const DegenAlphabet* abc;
	string name;

	/* non-member operators */
	friend DigitalSeq operator+(const DigitalSeq& lhs, const DigitalSeq& rhs);
	friend bool operator==(const DigitalSeq& lhs, const DigitalSeq& rhs);
	friend bool operator<(const DigitalSeq& lhs, const DigitalSeq& rhs);
	friend ostream& operator<<(ostream& os, const DigitalSeq& dSeq);
};

/* non-member operator implementations */
/*
 * operator+ implemented based on operator+=
 */
inline DigitalSeq operator+(const DigitalSeq& lhs, const DigitalSeq& rhs) {
	DigitalSeq ds(lhs); // make a copy
	ds += rhs;
	return ds;
}

/*
 * operator!= implemented based on operator==
 */
inline bool operator!=(const DigitalSeq& lhs, const DigitalSeq& rhs) {
	return !(lhs == rhs);
}

/*
 * operator<= implemented based on operator== and operator<
 */
inline bool operator<=(const DigitalSeq& lhs, const DigitalSeq& rhs) {
	return lhs < rhs || lhs == rhs;
}

/*
 * operator> implemented based on operator<=
 */
inline bool operator>(const DigitalSeq& lhs, const DigitalSeq& rhs) {
	return !(lhs <= rhs);
}

/*
 * operator>= implemented based on operator<
 */
inline bool operator>=(const DigitalSeq& lhs, const DigitalSeq& rhs) {
	return !(lhs < rhs);
}

} /* namespace EGriceLab */

#endif /* DIGITALSEQ_H_ */
