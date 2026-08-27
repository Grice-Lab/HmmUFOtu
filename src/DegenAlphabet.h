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
#include <cstdint>
#include <limits>
#include "HmmUFOtuDef.h"

namespace EGriceLab {
namespace HmmUFOtu {

using std::string;
using std::map;

class DegenAlphabet {
public:
	/* Constructors */
	/* explicitly disable default constructor */
	DegenAlphabet() = delete;
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

	/** test whether a base is a symbol/synonymous */
	bool isSymbol(int8_t b) const {
		return b >= 0;
	}

	/** test whether a char is a symbol/synonymous */
	bool isSymbol(char c) const {
		return isSymbol(encode(c));
	}

	/** test whether a base is a gap */
	bool isGap(int8_t b) const {
		return b == GAP_BASE;
	}

	/** test whether a char is a gap */
	bool isGap(char c) const {
		return isGap(encode(c));
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

	/** Get synonymous for a given symbol, or thorw exept if not exists */
	string getSynonymous(char c) const {
		return degen_map.at(c);
	}

	/** test whether a character is a synon */
	bool isSynonymous(char c) const {
		return degen_map.find(c) != degen_map.end();
	}

	/** test whether a base is valid */
	bool isValid(int8_t b) const {
		return b != INVALID_BASE;
	}

	/** test whether a char is valid */
	bool isValid(char c) const {
		return isValid(encode(c));
	}

	/* test whether two characters c1 and c2 is a match */
	bool isMatch(char c1, char c2) const;

	/* test whether a character is a match to a coded base */
	bool isMatch(char c, int8_t b) const;

	/* test whether a base is a match to a char */
	bool isMatch(int8_t b, char c) const {
		return isMatch(c, b);
	}

	/* pure virtual member method to be overridden by subclass */
	virtual bool hasComplement() const = 0;

	virtual char getComplementSymbol(char c) const = 0;

private:
	string name;
	string symbol; /* symbols of this alphabet */
	string synon; /* Expanded synonymous */
	string gap; /* gap characters */
	char gapCh {DEFAULT_GAP_CHAR}; /* representative gap char */
	int8_t sym_map[INT8_MAX + 1]; /* internal map for symbols */
	map<char, string> degen_map; // map for degenerative synonymous

public:
	static constexpr int8_t INVALID_BASE = -1;
	static constexpr int8_t GAP_BASE = -2; /* encoded gap symbol */
	static constexpr char DEFAULT_GAP_CHAR = '-';

	/* friend operators */
	friend bool operator==(const DegenAlphabet& lhs, const DegenAlphabet& rhs);
};

/* non-member operators */
inline bool operator!=(const DegenAlphabet& lhs, const DegenAlphabet& rhs) {
	return operator==(lhs, rhs);
}

} /* namespace HmHmmUFOtu */
} /* namespace EGriceLab */


#endif /* DEGENALPHABET_H_ */
