/*
 * IUPACNucl.h
 *
 *  Created on: May 5, 2015
 *      Author: zhengqi
 */

#ifndef IUPACNUCL_H_
#define IUPACNUCL_H_

#include <map>
#include <string>
#include <stdexcept>
#include <iostream>
#include "DegenAlphabet.h"

namespace EGriceLab {
using std::string;
using std::map;

class IUPACNucl: public DegenAlphabet {
public:
	/* Constructors */
	/* default constructor */
	IUPACNucl();

	/* destructor, do nothing */
	virtual ~IUPACNucl() { };

	/* member methods */

	/**
	 * Get alias of this alphabet
	 * @override  base class method
	 */
	virtual string getAlias() const {
		return "DNA";
	}

	/**
	 * tell whether has complement symbols, always true
	 * @override  base class method
	 */
	virtual bool hasComplement() const {
		return true;
	}
	/**
	 * Get the complement char of given symbol
	 * @return the complement symbol, or leave unchanged if not defined
	 */
	virtual char getComplementSymbol(char c) const {
		return compl_map[c];
	}

private:
	/* static initialization method */
	static map<char, string> init_IUPAC_map();
/*	static const map<char, string> degen_map;*/
	int8_t compl_map[INT8_MAX + 1];
};

} /* namespace EGriceLab */

#endif /* IUPACNUCL_H_ */
