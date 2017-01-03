/*
 * DNASubModelFactory.cpp
 *
 *  Created on: Dec 16, 2016
 *      Author: zhengqi
 */

#include "DNASubModelFactory.h"

namespace EGriceLab {

DNASubModel* DNASubModelFactory::createModel(const string& type) {
	if(type == "GTR")
		return new GTR();
	else
		throw std::invalid_argument("Unknown DNA Substitution Model type: " + type);
}

} /* namespace EGriceLab */


