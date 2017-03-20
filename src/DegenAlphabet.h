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
#include <stdint.h>

namespace EGriceLab {
using std::string;
using std::map;

class DegenAlphabet {
public:
	/* Constructors */
	/* customized constructors */
	/**
	 * Construct a DegenAlphabet with given name, symbol, expanded synonymous, and a map between sym and synon
	 * ambigus characters are resolved arbitrary
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
	const string& getName() const {
		return name;
	}

	virtual string getAlias() const {
		return getName();
	}

	const string& getSymbol() const {
		return symbol;
	}

	string getSynonymous() const {
		return synon;
	}

	const string& getGap() const {
		return gap;
	}

	const int8_t* getInMap() const {
		return sym_map;
	}

	/* utility methods */
	/* test whether a char is a symbol or synonymous */
	bool isSymbol(char c) const {
		return sym_map[c] >= 0;
	}

	/* encode a character to digital encoding
	 * return an int within 0..size-1, or -2 if a gap, or -1 if other invalid symbol
	 */
	int8_t encode(char c) const {
		return sym_map[c];
	}

	/* decode a digital encoding to the original symbol
	 * return a char if within 0..length-1,
	 * or gapCh if is gap_sym,
	 * or undefined behavior if other invalid values
	 */
	char decode(int8_t i) const {
		return i == GAP_BASE ? gapCh : symbol[i];
	}

	/*
	 * test whether a character is a valid gap
	 * @return true if is a valid gap character
	 */
	bool isGap(char c) const {
		return sym_map[c] == GAP_BASE;
	}

	int getSize() const {
		return symbol.length();
	}

	/*
	 * Get alphabet size w/ gaps
	 */
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


	/* Get synonymous for a given symbol, or empty string if not exists */
	string getSynonymous(char c) const {
		if(isSynonymous(c))
			return degen_map.find(c)->second;
		return "";
	}

	/* test whether a character is a degenerative synonymous */
	bool isSynonymous(char c) const {
		return degen_map.find(c) != degen_map.end();
	}

	/* test whether a character is a valid symbol or gap */
	bool isValid(char c) const {
		return sym_map[c] != INVALID_BASE;
	}

	/* pure virtual member method to be overridden by subclass */
	virtual bool hasComplement() const = 0;

	virtual char getComplementSymbol(char c) const = 0;

private:
	string name;
	string symbol; /* symbols of this alphabet */
	string synon; /* Expanded synonymous */
	string gap; /* gap characters */
	char gapCh; /* representative gap char */
	int8_t sym_map[INT8_MAX + 1]; /* internal map for symbols */
	//int8_t gap_map[INT8_MAX + 1]; /* internal map for gaps */
	map<char, string> degen_map; // map for degenerative synonymous

public:
	static const int8_t INVALID_BASE = -1;
	static const int8_t GAP_BASE = -2; /* encoded gap symbol */
	static const char DEFAULT_GAP_CHAR = '-';

	/* friend operators */
	friend bool operator==(const DegenAlphabet& lhs, const DegenAlphabet& rhs);
};

/* non-member operators */
inline bool operator!=(const DegenAlphabet& lhs, const DegenAlphabet& rhs) {
	return operator==(lhs, rhs);
}

} /* namespace EGriceLab */


#endif /* DEGENALPHABET_H_ */
