/*
 * Alphabet.h
 *
 *  Created on: May 4, 2015
 *      Author: zhengqi
 */

#ifndef ALPHABET_H_
#define ALPHABET_H_
#include <string>
#include <vector>
#include <climits>
#include <stdint.h>
//#include "HmmUFOtuDef.h"

namespace EGriceLab {
using std::string;

/**
 * An alphabet class using array as internal map.
 */
class Alphabet {
public:
	/* Constructors */
	/* Default constructor, do nothing */
	Alphabet() { };

	/** Construct an Alphabet with given name and symbol characters
	 * @param name  name of the Alphabet
	 * @param sym_str  string of alphabet characters, redundant characters will be removed
	 */
	Alphabet(const string& name, const string& sym_str);

	/* destructor */
	virtual ~Alphabet() { };

	/* utility methods */
	/* test whether a char is a symbol of this alphabet */
	bool isSymbol(char c) const {
		return sym_map[c] != invalid_sym;
	}

	/* encode a character to digital encoding
	 * return an int within 0..size-1, or -1 (string::nsize) if not a valid symbol
	 */
	virtual int8_t encode(char c) const {
		return sym_map[c];
	}

	/* decode a digital encoding to the original symbol
	 * return a char if within 0..length-1, or undefined behavior if not
	 */
	char decode(int8_t i) const {
		return sym[i];
	}

	/* Getters and Setters */
	const string& getName() const {
		return name;
	}
	int getSize() const {
		return sym.length();
	}
	string getSymbol() const {
		return sym;
	}

	const int8_t* getInmap() const {
		return sym_map;
	}

private:
	string name;
	string sym; // distinct symbols
	int8_t sym_map[INT8_MAX + 1]; // internal char to int8_t map for fast encoding and decoding
	static const int8_t invalid_sym = -1;

	/* friend operators */
	friend bool operator==(const Alphabet& lhs, const Alphabet& rhs);
};

/* non-member operators */
/**
 * Test whether two Alphabets are equal, only based on their symbols
 */
inline bool operator==(const Alphabet& lhs, const Alphabet& rhs) {
	return lhs.sym == rhs.sym;
}
/**
 * Test whether two Alphabets are not equal
 */
inline bool operator!=(const Alphabet& lhs, const Alphabet& rhs) {
	return !(lhs == rhs);
}
} /* namespace EGriceLab */


#endif /* ALPHABET_H_ */
