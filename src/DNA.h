/*
 * DNA.h
 *
 *  Created on: Oct 27, 2015
 *      Author: zhengqi
 */

#ifndef DNA_H_
#define DNA_H_

#include "DegenAlphabet.h"

namespace EGriceLab {

class DNA: public DegenAlphabet {
public:
	/* Constructors */
	/* default constructor */
	DNA();

	/* destructor, do nothing */
	virtual ~DNA() { };

	/* member methods */
	/* implementation of abstract superclass methods */
	/**
	 * always return true
	 */
	bool hasComplement() const {
		return true;
	}
	/**
	 * Get the complement char of given symbol
	 * @return the complement symbol, or '\0' if not a valid symbol
	 */
	char getComplementSymbol(char c) const {
		return compl_map[c];
	}

private:
	/* static initialization method */
	static map<char, string> init_DNA_map();
/*	static const map<char, string> degen_map;*/
	//map<char, char> compl_map;
	char compl_map[INT8_MAX + 1];
};

} /* namespace EGriceLab */

#endif /* DNA_H_ */
