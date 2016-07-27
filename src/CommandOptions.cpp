/*
 * CommandOptions.cpp
 *
 *  Created on: Jul 15, 2016
 *      Author: zhengqi
 */

#include "CommandOptions.h"
#include <iostream>

namespace EGriceLab {

CommandOptions::CommandOptions(int argc, char** argv) {
	/* parse options */
	for(int i = 1; i < argc; ++i) {
		if(*argv[i] == '-') { /* a tag name */
			if(i < argc - 1 && *argv[i+1] != '-') {/* a tag value */
				opts[argv[i]] = argv[i+1];
				i++;
			}
			else /* a flag tag */
				opts[argv[i]] = ""; /* use empty value */
		}
		else /* a main opt */
			mainOpts.push_back(argv[i]);
	}
}

} /* namespace EGriceLab */

