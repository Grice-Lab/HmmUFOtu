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
	/* test whether a char is a symbol of this alphabet */
	bool isSymbol(char c) const {
		return sym_map[c] != invalid_sym;
	}

	/* encode a character to digital encoding
	 * return an int within 0..size-1, or -1 (string::nsize) if not a valid symbol
	 */
	int8_t encode(char c) const {
		return sym_map[c];
	}

	/* decode a digital encoding to the original symbol
	 * return a char if within 0..length-1, or undefined behavior if not
	 */
	char decode(int8_t i) const {
		return symbol[i];
	}

	/*
	 * test whether a character is a valid gap
	 * @return true if is a valid gap character
	 */
	bool isGap(char c) const {
		return gap_map[c] != 0;
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

	/* test whether a character is a valid symbol or synonymous */
	bool isValid(char c) const {
		return isSymbol(c) || isSynonymous(c) || isGap(c);
	}

	/* pure virtual member method to be overridden by subclass */
	virtual bool hasComplement() const = 0;

	virtual char getComplementSymbol(char c) const = 0;

private:
	string name;
	string symbol; /* symbols of this alphabet */
	string synon; /* Expanded synonymous */
	string gap; /* gap characters */
	int8_t sym_map[INT8_MAX + 1]; /* internal map for symbols */
	int8_t gap_map[INT8_MAX + 1]; /* internal map for gaps */
	map<char, string> degen_map; // map for degenerative synonymous, including original symbols' self-mapping

	static const int8_t invalid_sym = -1;
	//static const int8_t invalid_synon = -2;
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
