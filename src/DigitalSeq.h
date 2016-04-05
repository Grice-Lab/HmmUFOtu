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
class DigitalSeq: public std::basic_string<int> {
public:
	/* constructors */
	/** Construct a DigitalSeq with given alphabet, name and string
	 * @param dgAbc  A DegenAlphabet
	 * @param name  name of this ds
	 * @param str  string of this ds, non-symbol characters will be discarded; synomynous will be resolved randomly
	 * (you need to call srand() in the main function)
	 * ignore any non-symbol, non-synoymous characters in the string
	 */
	DigitalSeq(const DegenAlphabet* abc, const string& name, const string& str = "");

	/**
	 Construct a DigitalSeq from a PrimrarySeq
	 */
	DigitalSeq(const PrimarySeq& seq);

	/* virtual destructor */
	virtual ~DigitalSeq() { }

	/* Getters and Setters */
	const DegenAlphabet* getAbc() const {
		return abc;
	}

	const string& getName() const {
		return name;
	}

	void setName(const string& name) {
		this->name = name;
	}

	/* utility member methods */
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
	 * Append a new string to this DigitalSeq
	 * return the modified *this
	 */
	DigitalSeq& append(const string& str);

private:
	const DegenAlphabet* const abc;
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
