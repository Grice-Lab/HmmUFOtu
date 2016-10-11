/*
 * IUPACAmino.h
 *
 *  Created on: May 5, 2015
 *      Author: zhengqi
 */

#ifndef IUPACAMINO_H_
#define IUPACAMINO_H_

#include <map>
#include <string>
#include <stdexcept>
#include "DegenAlphabet.h"

namespace EGriceLab {
using std::string;
using std::map;

class IUPACAmino: public DegenAlphabet {
public:
	/* Constructors */
	/* default constructor */
	IUPACAmino() : DegenAlphabet("IUPACAmino",
			"ACDEFGHIKLMNPQRSTVWY", "BXZ", init_IUPAC_map()) {
	}

	/* destructor, do nothing */
	virtual ~IUPACAmino() { };

	/* member methods */
	/**
	 * Get alias of this alphabet
	 * @override  base class method
	 */
	virtual string getAlias() const {
		return "AMINO";
	}

	/**
	 * check whether has complement, always false
	 * @override  base class method
	 */
	virtual bool hasComplement() const {
		return false;
	}

	/**
	 * Get the complement char of given symbol
	 * @return unchanged amino acids don't have complementary symbols
	 */
	virtual char getComplementSymbol(char c) const {
		return c;
	}

private:
	/* static initialization method */
	static map<char, string> init_IUPAC_map();
};

} /* namespace EGriceLab */

#endif /* IUPACAMINO_H_ */
