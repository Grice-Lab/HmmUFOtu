/*
 * VersionSequence.h
 *  A universal version sequence class for defining program and product version
 *  Created on: Oct 4, 2017
 *      Author: zhengqi
 */

#ifndef SRC_VERSIONSEQUENCE_H_
#define SRC_VERSIONSEQUENCE_H_

#include <string>
#include <iostream>
#include <strstream>
#include <stdexcept>

namespace EGriceLab {
using std::string;
using std::istream;
using std::ostream;
using std::invalid_argument;

class VersionSequence {
public:
	/* constructors */

	/** default constructor */
	VersionSequence() : major(0), minor(0), build(0) {  }

	/** construct from a given version string, if failed, using default values */
	VersionSequence(const string& str)
		: major(0), minor(0), build(0)
	{
		parseString(str, *this);
	}

	/** destructor, do nothing */
	virtual ~VersionSequence() {  }

	/* member methods */
	string toString() const {
		return string("v") + major + "." + minor + "." + build;
	}

	/** load ProgVer from input, with null termination */
	istream& load(istream& in);

	/** save ProgVer to output with null termination */
	ostream& save(ostream& out) const;

	/* non-member methods */
	/* write operator with non-null terminator */
	friend istream& operator>>(istream& in, VersionSequence& ver);

	/* read operator with non-null terminator */
	friend ostream& operator<<(ostream& out, const VersionSequence& ver);

	friend bool operator==(const VersionSequence& lhs, const VersionSequence& rhs);

	friend bool operator<(const VersionSequence& lhs, const VersionSequence& rhs);

	/* member fields */
private:
	int major;
	int minor;
	int build;

	/* static fields */
public:
	static const int MAX_LENGTH = 4096;

	/* static methods */
	static void parseString(const string& str, VersionSequence& ver);

};

inline ostream& operator<<(ostream& out, const VersionSequence& ver) {
	return out << ver.toString();
}

inline bool operator==(const VersionSequence& lhs, const VersionSequence& rhs) {
	return lhs.major == rhs.major && lhs.minor == rhs.minor && lhs.build == rhs.build;
}

inline bool operator<(const VersionSequence& lhs, const VersionSequence& rhs) {
	return lhs.major < rhs.major ||
			lhs.major == rhs.major && lhs.minor < rhs.minor ||
			lhs.major == rhs.major && lhs.minor == rhs.minor && lhs.build < rhs.build;
}

inline bool operator!=(const VersionSequence& lhs, const VersionSequence& rhs) {
	return !(lhs == rhs);
}

inline bool operator<=(const VersionSequence& lhs, const VersionSequence& rhs) {
	return lhs < rhs || lhs == rhs;
}

inline bool operator>=(const VersionSequence& lhs, const VersionSequence& rhs) {
	return !(lhs < rhs);
}

inline bool operator>(const VersionSequence& lhs, const VersionSequence& rhs) {
	return !(lhs <= rhs);
}

} /* namespace EGriceLab */

#endif /* SRC_VERSIONSEQUENCE_H_ */
