/*
 * DegenAlphabet.h
 *
 *  Created on: May 5, 2015
 *      Author: zhengqi
 */

#ifndef DEGENALPHABET_H_
#define DEGENALPHABET_H_

#include <map>
#include <string>
#include <stdexcept>
#include "Alphabet.h"

namespace EGriceLab {
using std::string;
using std::map;

class DegenAlphabet: public Alphabet {
public:
	/* Constructors */
	/* customized constructors */
	/**
	 * Construct a DegenAlphabet with given name, symbol, expanded synonymous, and a map between sym and synon
	 * @param name  DegenAlphabet name
	 * @param sym_str  symbol string, redundant characters will be removed
	 * @param synon_str  additional synonymous characters, redundant will be removed
	 * @param my_map  degenerative mapping between synon char to multiple sym as a string
	 */
	DegenAlphabet(const string& name, const string& sym_str, const string& synon_str,
			const map<char, string>& my_map, const string& gap = "-._");

	/* virtual destructor, do nothing */
	virtual ~DegenAlphabet() { }

	/* Member methods */
	/* Getters and Setters */
	const string& getGap() const {
		return gap;
	}

	void setGap(const string& gap) {
		this->gap = gap;
	}

	/* test whether a char is a valid gap */
	bool isGap(char c) const {
		return gap_map[c] != 0;
	}

	/* Get alphabet size w/ gaps */
	int getSizeWithGap() const {
		return getSize() + gap.length();
	}

	/* Get size with all degenerative synonymous */
	int getDegenSize() const {
		return getSize() + synon.length();
	}

	/* Get size with all degenerative synonymous and gaps*/
	int getDegenSizeWithGap() const {
		return getDegenSize() + gap.length();
	}

	/* Get degenerative synonymous */
	string getSynonymous() const {
		return synon;
	}

	/* Get synonymous for a given symbol, or empty string if not exists */
	string getSynonymous(char c) const {
		if(isSynonymous(c))
			return degen_map.find(c)->second;
		return "";
	}

	/* encode a character to digital encoding
	 * synonymous is resolved randomly (you will need to call srand() in your main function once)
	 * @param c  given character
	 * @return value within 0..size-1, or a minus value if not a valid symbol
	 */
	virtual int8_t encode(char c) const {
		return !isSynonymous(c) ? Alphabet::encode(c) : synon_map[c];
	}

	/* test whether a character is a degenerative synonymous */
	bool isSynonymous(char c) const {
		return synon_map[c] != invalid_synon;
	}

	/*
	 *  test whether a character is a pure synonymous
	 */
	bool isPureSynonymous(char c) const {
		return !isSymbol(c) && isSynonymous(c);
	}

	/*
	 * test whether a character is a symbol or synonymous
	 */
	bool isSymbolOrSynonymous(char c) const {
		return isSymbol(c) || isSynonymous(c);
	}

	/* test whether a character is a valid symbol or synonymous */
	bool isValid(char c) const {
		return isSymbol(c) || isSynonymous(c) || isGap(c);
	}

	/* pure virtual member method to be overridden by subclass */
	virtual bool hasComplement() const = 0;

	virtual char getComplementSymbol(char c) const = 0;

private:
	string synon; /* Expanded symbols + synonymous */
	string gap; /* gap characters */
	int8_t synon_map[INT8_MAX + 1]; /* internal map for synonymous only */
	int8_t gap_map[INT8_MAX + 1]; /* internal map for gaps */
	map<char, string> degen_map; // map for degenerative synonymous, including original symbols' self-mapping

	static const int8_t invalid_synon = -2;
	//static const int8_t invalid_gap = -3;

	/* friend operators */
	friend bool operator==(const DegenAlphabet& lhs, const DegenAlphabet& rhs);
};

/* non-member operators */
inline bool operator!=(const DegenAlphabet& lhs, const DegenAlphabet& rhs) {
	return operator==(lhs, rhs);
}

} /* namespace EGriceLab */


#endif /* DEGENALPHABET_H_ */
