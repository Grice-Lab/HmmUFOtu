/*
 * DNASubModelFactory.h
 *  A factory class with static methods creating new empty DNASubstitution models base on type names
 *  Created on: Dec 16, 2016
 *      Author: zhengqi
 */

#ifndef SRC_DNASUBMODELFACTORY_H_
#define SRC_DNASUBMODELFACTORY_H_

#include <stdexcept>
#include <string>
#include "DNASubModel.h"
#include "GTR.h"

namespace EGriceLab {

using std::string;

class DNASubModelFactory {
public:
	static DNASubModel* createModel(const string& type);
};

} /* namespace EGriceLab */

#endif /* SRC_DNASUBMODELFACTORY_H_ */
