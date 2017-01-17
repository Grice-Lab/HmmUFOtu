/*
 * CommandOptions.h
 *
 *  Created on: Jul 15, 2016
 *      Author: zhengqi
 */

#ifndef SRC_COMMANDOPTIONS_H_
#define SRC_COMMANDOPTIONS_H_

#include <string>
#include <set>
#include <map>
#include <vector>

namespace EGriceLab {

using std::string;
using std::map;
using std::vector;

class CommandOptions {
public:
	CommandOptions(int argc, char** argv);

	bool hasOpt(const string& name) const {
		return opts.find(name) != opts.end();
	}

	string getOpt(const string& name) const {
		return hasOpt(name) ? opts.find(name)->second : "";
	}

	const char* getOptStr(const string& name) const {
		return hasOpt(name) ? opts.find(name)->second.c_str() : "";
	}

	int numMainOpts() const {
		return mainOpts.size();
	}

	string getMainOpt(int i) const {
		return mainOpts.at(i);
	}

private:
	vector<string> mainOpts; /* mandatory options not following any -tag-name */
	map<string, string> opts; /* optional named options in -tag-name [value] pairs */
};

} /* namespace EGriceLab */

#endif /* SRC_COMMANDOPTIONS_H_ */
