/*
 * PrimarySeq.h
 *
 *  Created on: Jun 26, 2015
 *      Author: Qi Zheng
 */

#ifndef PRIMARYSEQ_H_
#define PRIMARYSEQ_H_

#include <string>
#include <iostream>
#include <stdexcept>
#include <climits>

#include "AlphabetFactory.h"
#include "StringUtils.h"

namespace EGriceLab {
using std::string;
using std::invalid_argument;
/**
 * A base class stands for a Biological sequence, similar to that of a the Bio::PrimarySeq class in BioPerl
 * Note the life-span of a PrimarySeq depends on the underlying alphabet object
 * @version v1.1
 * @since v1.1
 */
class PrimarySeq {

public:
	/* constructors */
	/**
	 * Construct a PrimarySeq with given alphabet, id, seq and optionally description
	 * @param alphabet  name of the alphabet
	 * @param id  display id
	 * @param seq  sequence
	 * @param desc  brief description
	 * @throw std::invalid_argument exception if the {@param seq} contains invalid alphabet characters
	 */
	PrimarySeq(const string& alphabet, const string& id, const string& seq,
			const string& desc = "", const string& qual = "") :
	abc(AlphabetFactory::getAlphabetByName(alphabet)), id(id), seq(seq),
	desc(desc), qual(qual), phredShift(DEFAULT_PHRED_SHIFT) {
		if(!isValidate())
			throw invalid_argument("Your sequence '" + seq + " ' contains invalid alphabet characters");
		if(!qual.empty() && qual.length() != seq.length())
			throw invalid_argument("qual length must be the same as seq length");
	}

	/**
	 * Construct a PrimarySeq with given alphabet pointer, id, seq and optionally description
	 * @param abc  pointer to an alphabet
	 * @param id  display id
	 * @param seq  sequence
	 * @param desc  brief description
	 * @throw std::invalid_argument exception if the {@param seq} contains invalid alphabet characters
	 */
	PrimarySeq(const DegenAlphabet* abc, const string& id, const string& seq,
			const string& desc = "", const string& qual = "") :
	abc(abc), id(id), seq(seq),
	desc(desc), qual(qual), phredShift(DEFAULT_PHRED_SHIFT) {
		if(!isValidate())
			throw invalid_argument("Your sequence '" + seq + " ' contains invalid alphabet characters");
		if(!qual.empty() && qual.length() != seq.length())
			throw invalid_argument("qual length must be the same as seq length");
	}

	/**
	 * destructor, do nothing
	 */
	virtual ~PrimarySeq() {  }

	/* Getters and Setters */
	const DegenAlphabet* getAbc() const {
		return abc;
	}

	const string& getDesc() const {
		return desc;
	}

	void setDesc(const string& desc) {
		this->desc = desc;
	}

	const string& getId() const {
		return id;
	}

	void setId(const string& id) {
		this->id = id;
	}

	const string& getSeq() const {
		return seq;
	}

	void setSeq(const string& seq) {
		this->seq = seq;
		if(!isValidate())
			throw invalid_argument(string("Cannot set seq to '") + seq + " ' that contains invalid alphabet characters");
	}

	string getQual() const {
		if(!qual.empty())
			return qual;
		else
			return string(length(), DEFAULT_QUAL + phredShift);
	}

	void setQual(const string& qual) {
		if(!qual.empty() && qual.length() != seq.length())
			throw invalid_argument("qual length must be the same as the seq length");
		this->qual = qual;
	}

	int getPhredShift() const {
		return phredShift;
	}

	void setPhredShift(int phredShift) {
		if(!(phredShift >= 0 && phredShift <= UCHAR_MAX))
			throw invalid_argument("phredShift not in valid range");
		this->phredShift = phredShift;
	}

	/* member functions */
	/**
	 * get the length of this PrimarySeq
	 * @return  the length of the underlying string
	 */
	string::size_type length() const {
		return seq.length();
	}

	/**
	 * test whether this seq is empty
	 * @return true if the underlying string is empty
	 */
	bool empty() const {
		return seq.empty();
	}

	/**
	 * get the total gaps of this PrimarySeq
	 * @return  number of gaps in this seq
	 */
	string::size_type numGap() const;

	/**
	 * get the non-gap length of this PrimarySeq
	 * @return  non-gap bases in this seq
	 */
	string::size_type nonGapLength() const {
		return length() - numGap();
	}

	/**
	 * validate this PrimarySeq
	 */
	bool isValidate() const;

	/**
	 * Modify this PrimarySeq internal seq to all upper-case
	 * @return the modified object
	 */
	PrimarySeq& toUpper() {
		StringUtils::toUpper(seq);
		return *this;
	}

	/**
	 * Modify this PrimarySeq internal seq to all lower-case
	 * @return the modified object
	 */
	PrimarySeq& toLower() {
		StringUtils::toLower(seq);
		return *this;
	}

	/**
	 * Remove gaps of this seq
	 * @return  modified this object with gaps removed
	 */
	PrimarySeq& removeGaps();

	/**
	 * Generate a reverse-complement copy of this seq
	 * @return a reverse-complement new seq
	 * @throw logic error if the alphabet doesn't support this operation
	 */
	PrimarySeq revcom() const;

	/**
	 * Get the subseq string of this PrimarySeq
	 * @return a subseq string
	 */
	string subseq(string::size_type pos, string::size_type len) const;

	/**
	 * Return a trucated copy of this PrimarySeq
	 * @return a new truncated copy
	 */
	PrimarySeq trunc(string::size_type pos, string::size_type len) const;

	/**
	 * Get the character at given pos
	 * @param pos  relative pos
	 * @return character at pos
	 * @throw out_of_range exception if pos is out of range
	 */
	char charAt(string::size_type pos) const;

	/**
	 * Get the alphabet encoded value of char at given pos
	 * @param pos  relative pos
	 * return 0..alphabet_size-1 of the char at pos, or negative value if not a valid character
	 * @throw out_of_range exception if pos is out of range
	 */
	int8_t encodeAt(string::size_type pos) const;

	/**
	 * Get the Phread Q-score at given pos
	 * @param pos  relative pos
	 * @return phred Q-score here
	 * @throw out_of_range exception if pos is out of range
	 */
	int qScoreAt(string::size_type pos) const;

	/* non-member friend operators */
	friend bool operator==(const PrimarySeq& lhs, const PrimarySeq& rhs);


private:
	const DegenAlphabet* abc; /* alphabet is not changeable after initialization */
	string id;
	string seq;
	string desc;
	string qual;
	int phredShift;

	/* static members */
	static const int DEFAULT_QUAL = 30;
	static const int DEFAULT_PHRED_SHIFT = 33;
};

inline string PrimarySeq::subseq(string::size_type pos,
		string::size_type len) const {
	return seq.substr(pos, len);
}

inline char PrimarySeq::charAt(string::size_type pos) const {
	return seq.at(pos);
}

inline int8_t PrimarySeq::encodeAt(string::size_type pos) const {
	return abc->encode(seq.at(pos));
}

inline int PrimarySeq::qScoreAt(string::size_type pos) const {
	return !qual.empty() ? qual.at(pos) - phredShift : DEFAULT_QUAL;
}

inline bool operator==(const PrimarySeq& lhs, const PrimarySeq& rhs) {
	return *lhs.abc == *rhs.abc && lhs.id == rhs.id
			&& StringUtils::toUpper(lhs.seq) == StringUtils::toUpper(rhs.seq);
}

inline bool operator!=(const PrimarySeq& lhs, const PrimarySeq& rhs) {
	return !(lhs == rhs);
}

inline PrimarySeq PrimarySeq::trunc(string::size_type pos, string::size_type len) const {
	return PrimarySeq(abc, id, seq.substr(pos, len), desc, qual);
}

} /* namespace EGriceLab */

#endif /* PRIMARYSEQ_H_ */
