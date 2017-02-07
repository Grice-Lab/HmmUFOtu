/*
 * BandedHMMP7.cpp
 *
 *  Created on: May 13, 2015
 *      Author: zhengqi
 */

#include <math.h> /* using C99 */
#include <cstdlib>
#include <cassert>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <stdexcept>
#include <ctime>
#include <sstream>
#include <algorithm>
#include "BandedHMMP7.h"
#include "LinearAlgebraBasic.h"

namespace EGriceLab {
using namespace std;
using namespace Eigen;

/*const int BandedHMMP7::kMinProfile = 10000; // up-to 10K 16S rRNA profile*/
string HMM_TAG =
		"HMM\t\tA\tC\tG\tT\n\t\tm->m\tm->i\tm->d\ti->m\ti->i\td->m\td->d";

const double BandedHMMP7::kMinGapFrac = 0.2;
const double BandedHMMP7::CONS_THRESHOLD = 0.9;
const double BandedHMMP7::DEFAULT_ERE = 1;
const IOFormat tabFmt(StreamPrecision, DontAlignCols, "\t", "\n", "", "", "", "");

BandedHMMP7::BandedHMMP7() :
		hmmVersion(progName + "-" + progVersion), name("unnamed"), K(0), abc(NULL),
		hmmBg(0), nSeq(0), effN(0),
		transInit(false), emisInit(false),
		spInit(false), limitInit(false), indexInit(false), wingRetracted(false) {
	/* Assert IEE559 at construction time */
	assert(numeric_limits<double>::is_iec559);
}

BandedHMMP7::BandedHMMP7(const string& name, int K, const DegenAlphabet* abc) :
		hmmVersion(progName + "-" + progVersion), name(name), K(K), abc(abc),
		hmmBg(K), nSeq(0), effN(0),
		transInit(false), emisInit(false), spInit(false), limitInit(false), indexInit(false),
		wingRetracted(false) {
	if(!(abc->getName() == "IUPACNucl" && abc->getSize() == 4))
		throw invalid_argument("BandedHMMP7 only supports DNA alphabet");
	/* Assert IEE559 at construction time */
	assert(numeric_limits<double>::is_iec559);
	init_transition_params();
	init_emission_params();
	init_special_params();
	init_index();
	init_limits();
	enableProfileLocalMode(); // always in profile local alignment mode
	setSpEmissionFreq(hmmBg.getBgEmitPr()); // Use bg emission freq by default
}

/* non-member friend functions */
istream& operator>>(istream& in, BandedHMMP7& hmm) {
	string line;
	int k = 0; // pos on the profile
	while (getline(in, line)) {
		if (line == "//") {/* end of profile */
			hmm.resetProbByCost(); // set the cost matrices
			hmm.adjustProfileLocalMode();
			hmm.wingRetract();
			return in;
		}
		istringstream iss(line); // detail parse this line
		string tag; /* header tag names and values */
		string tmp;
		if (!isspace(line[0])) { /* header section starts with non-empty characters */
			iss >> tag;
			if (tag.substr(0, 6) == "HMMER3") { // do not override our version, check minor version
				if(tag.length() < 8 || tag[7] < 'f') {
					cerr << "Obsolete HMM file version: " << tag << ", must be HMMER3/f or higher" << endl;
					in.setstate(ios_base::badbit);
					return in;
				}
			}
			else if (tag == "NAME") {
				iss >> hmm.name;
			} else if (tag == "LENG") {
				int size;
				iss >> size;
				hmm.setProfileSize(size);
			} else if (tag == "ALPH") {
				string abc;
				iss >> abc;
				if (abc != "DNA")
					throw invalid_argument(
							"Not allowed alphabet '" + abc
									+ "' in the HMM input file! Must be DNA");
				// override the alphabet
				hmm.abc = AlphabetFactory::getAlphabetByName("DNA");
			} else if (tag == "STATS") {
				string mode;
				string distrib;
				iss >> mode >> distrib;
				tag += " " + mode + " " + distrib; // use STATS + mode + distribution as the new tag name
				string val;
				getline(iss, val);
				hmm.setOptTag(tag, BandedHMMP7::trim(val));

			} else if(tag == "HMM") { /* HMM TAG */
				string tmp;
				getline(in, tmp); /* ignore the next line too */
			}
			else { /* optional tags */
				string val;
				getline(iss, val); // get the entire remaining part of this line
				if(!tag.empty())
					hmm.setOptTag(tag, BandedHMMP7::trim(val)); // record this tag-value pair
				// check some optional tags
				if(tag == "NSEQ")
					hmm.nSeq = ::atoi(val.c_str());
				else if(tag == "EFFN")
					hmm.effN = ::atof(val.c_str());
				else
				{ /* do nothing */ }
			}
		} /* end of header section */
		else { /* Main body, starts with space */
			iss >> tag;
			if (tag == "COMPO" || BandedHMMP7::isInteger(tag)) { // A compo line can be treated as position 0
				assert((tag == "COMPO" && k == 0) || atoi(tag.c_str()) == k);
				/* process current emission line */
				Vector4d emitFreq;
				for (Vector4d::Index i = 0; i < 4; ++i)
					iss >> emitFreq(i);
				if (tag == "COMPO") { // COMPO line
					hmm.E_M_cost.col(0) = emitFreq;
					emitFreq = (-emitFreq).array().exp();
					hmm.setSpEmissionFreq(emitFreq);
					hmm.hmmBg.setBgFreq(emitFreq);
				} else {
					/* Mk emission line */
					hmm.E_M_cost.col(k) = emitFreq;
					/* Make sure the MAP tag is set */
					string val;
					if(hmm.getOptTag("MAP") != "yes") {
						cerr << "Error: HMM file must has the MAP flag set to 'yes'" << endl;
						in.setstate(ios_base::badbit);
						return in;
					}
					iss >> tmp;
					hmm.cs2ProfileIdx[atoi(tmp.c_str())] = k;
					hmm.profile2CSIdx[k] = atoi(tmp.c_str());
					hmm.setLocOptTag("MAP", tmp, k);
					/* read other optional tags */
					if(!hmm.getOptTag("CONS").empty()) { /* this tag is present, regarding yes or no */
						iss >> tmp;
						hmm.setLocOptTag("CONS", tmp, k);
					}
					if(!hmm.getOptTag("RF").empty()) { /* this tag is present, regarding yes or no */
						iss >> tmp;
						hmm.setLocOptTag("RF", tmp, k);
					}
					if(!hmm.getOptTag("MM").empty()) { /* this tag is present, regarding yes or no */
						iss >> tmp;
						hmm.setLocOptTag("MM", tmp, k);
					}
					if(!hmm.getOptTag("CS").empty()) { /* this tag is present, regarding yes or no */
						iss >> tmp;
						hmm.setLocOptTag("CS", tmp, k);
					}
				}
				/* process the following Ik emission line */
				for (MatrixXd::Index i = 0; i < hmm.E_I_cost.rows(); ++i)
					in >> hmm.E_I_cost(i, k);
				/* process the following state K transition line */
					in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::M, BandedHMMP7::M) = hmm.hmmValueOf(tmp);  // Mk -> Mk+1
					in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::M, BandedHMMP7::I) = hmm.hmmValueOf(tmp);  // Mk -> Ik
					in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::M, BandedHMMP7::D) = hmm.hmmValueOf(tmp);  // Mk -> Dk+1
					in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::I, BandedHMMP7::M) = hmm.hmmValueOf(tmp);  // Ik -> Mk+1
					in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::I, BandedHMMP7::I) = hmm.hmmValueOf(tmp);  // Ik -> Ik
					in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::D, BandedHMMP7::M) = hmm.hmmValueOf(tmp);  // Dk -> Mk+1
					in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::D, BandedHMMP7::D) = hmm.hmmValueOf(tmp);  // Dk -> Dk+1
			} /* combo line section or match state line section */
			else { // non-COMPO begin state line (M0)
				assert(k == 0);
				string tmp;
				/* process the BEGIN insert emission line */
				for (MatrixXd::Index i = 0; i < hmm.E_I_cost.rows(); ++i)
					in >> hmm.E_I_cost(i, k);
				/* process the B state K transition line */
				in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::M, BandedHMMP7::M) = hmm.hmmValueOf(tmp);  // Mk -> Mk+1
				in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::M, BandedHMMP7::I) = hmm.hmmValueOf(tmp);  // Mk -> Ik
				in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::M, BandedHMMP7::D) = hmm.hmmValueOf(tmp);  // Mk -> Dk+1
				in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::I, BandedHMMP7::M) = hmm.hmmValueOf(tmp);  // Ik -> Mk+1
				in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::I, BandedHMMP7::I) = hmm.hmmValueOf(tmp);  // Ik -> Ik
				in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::D, BandedHMMP7::M) = hmm.hmmValueOf(tmp);  // Dk -> Mk+1
				in >> tmp; hmm.Tmat_cost[k](BandedHMMP7::D, BandedHMMP7::D) = hmm.hmmValueOf(tmp);  // Dk -> Dk+1
			}
			k++;
		} /* end of main section */
	} /* end of each line */
	// somehow the hmm file reached end without '//'
	in.setstate(std::ios::failbit);
	return in;
}

void BandedHMMP7::scale(double r) {
	/* scale transitions */
	for(int k = 0; k <= K; ++k)
		Tmat[k] *= r;
	/* scale emissions */
	E_M *= r;
	E_I *= r;
	/* reset costs */
	resetCostByProb();
}

void BandedHMMP7::normalize() {
	for(int k = 0; k <= K; ++k) {
		/* normalize transitions */
		Tmat[k].row(M) /= Tmat[k].row(M).sum(); /* TMX */
		Tmat[k].row(I) /= Tmat[k].row(I).sum(); /* TIX */
		Tmat[k].row(D) /= Tmat[k].row(D).sum(); /* TDX */
		/* normalize emissions */
		E_M.col(k) /= E_M.col(k).sum(); /* EM */
		E_I.col(k) /= E_I.col(k).sum(); /* EI */
	}
	/* enforce the T[0] and T[K] specials */
	Tmat[0](D, M) = 1;
	Tmat[0](D, D) = 0;
	Tmat[K](M, D) = 0;
	Tmat[K](D, M) = 1;
	Tmat[K](D, D) = 0;

	/* reset costs */
	resetCostByProb();
}

void BandedHMMP7::estimateParams(const BandedHMMP7Prior& prior) {
	assert(abc->getSize() == prior.dmME.getK());

	/* normalize the COMPO Match emission, which is the B state emission */
//	E_M.col(0) /= E_M.col(0).sum();

	/* re-estimate parameters using the prior info */
	for(int k = 0; k <= K; ++k) {
		/* update transition parameters */
		/* TM */
		Tmat[k].row(M) = prior.dmMT.meanPostP(Tmat[k].row(M));
		/* TI */
		Tmat[k].row(I).segment(M, 2) = prior.dmIT.meanPostP(Tmat[k].row(I).segment(M, 2)); /* only use first two elements of the TI row */
		/* TD */
		VectorXd dt(2);
		dt(0) = Tmat[k](D, M);
		dt(1) = Tmat[k](D, D);
		dt = prior.dmDT.meanPostP(dt); /* replace observed frequency with meanPostP */
		Tmat[k](D, M) = dt(0);
		Tmat[k](D, D) = dt(1);

		/* update emission parameters */
		E_M.col(k) = prior.dmME.meanPostP(E_M.col(k));
		E_I.col(k) = prior.dmIE.meanPostP(E_I.col(k));
	}

	/* enforce the T[0] and T[K] specials */
	Tmat[0](D, M) = 1;
	Tmat[0](D, D) = 0;
	Tmat[K](M, D) = 0;
	Tmat[K](D, M) = 1;
	Tmat[K](D, D) = 0;

	/* reset costs */
	resetCostByProb();
}

double BandedHMMP7::meanRelativeEntropy() const {
	double ent = 0;
	for(int k = 1; k <= K; ++k)
		ent += Math::relative_entropy(E_M.col(k), hmmBg.getBgEmitPr());
	return ent / K;
}

ostream& operator<<(ostream& out, const BandedHMMP7& hmm) {
	/* write mandatory tags */
	out << "HMMER3/f\t" << hmm.hmmVersion << endl;
	out << "NAME\t" << hmm.name << endl;
	out << "LENG\t" << hmm.K << endl;
	out << "ALPH\t" << hmm.abc->getAlias() << endl;

	/* write optional tags */
	for(vector<string>::const_iterator it = hmm.optTagNames.begin(); it != hmm.optTagNames.end(); ++it)
		out << *it << "  " << hmm.getOptTag(*it) << endl;

	/* write optional HMM tags */
	out << HMM_TAG << endl;
	for(int k = 0; k <= hmm.K; ++k) {
		/* write M or background emission line */
		if(k == 0)
			out << "\tCOMPO\t" << hmm.E_M_cost.col(0).transpose().format(tabFmt) << endl;
		else {
			out << "\t" << k << "\t" << hmm.E_M_cost.col(k).transpose().format(tabFmt);
			/* write other optional tags, if present */
			if(!hmm.getOptTag("MAP").empty())
				out << "\t" << hmm.getLocOptTag("MAP", k);
			if(!hmm.getOptTag("CONS").empty())
				out << "\t" << hmm.getLocOptTag("CONS", k);
			if(!hmm.getOptTag("RF").empty())
				out << "\t" << hmm.getLocOptTag("RF", k);
			if(!hmm.getOptTag("MM").empty())
				out << "\t" << hmm.getLocOptTag("MM", k);
			if(!hmm.getOptTag("CS").empty())
				out << "\t" << hmm.getLocOptTag("CS", k);
			out << endl;
		}
		/* write insert emission line */
		double val;
		out << "\t";
		for(MatrixXd::Index i = 0; i != hmm.E_I_cost.rows(); ++i) {
			val = hmm.E_I_cost(i, k);
			hmmPrintValue(out << "\t", val);
		}
		out << endl;

		/* write state transition line */
		val = hmm.Tmat_cost[k](BandedHMMP7::M, BandedHMMP7::M); hmmPrintValue(out << "\t\t", val);
		val = hmm.Tmat_cost[k](BandedHMMP7::M, BandedHMMP7::I); hmmPrintValue(out << "\t", val);
		val = hmm.Tmat_cost[k](BandedHMMP7::M, BandedHMMP7::D); hmmPrintValue(out << "\t", val);
		val = hmm.Tmat_cost[k](BandedHMMP7::I, BandedHMMP7::M); hmmPrintValue(out << "\t", val);
		val = hmm.Tmat_cost[k](BandedHMMP7::I, BandedHMMP7::I); hmmPrintValue(out << "\t", val);
		val = hmm.Tmat_cost[k](BandedHMMP7::D, BandedHMMP7::M); hmmPrintValue(out << "\t", val);
		val = hmm.Tmat_cost[k](BandedHMMP7::D, BandedHMMP7::D); hmmPrintValue(out << "\t", val);

		out << endl;
	}
	out << "//" << endl;
	return out;
}

ostream& operator<<(ostream& os, const deque<BandedHMMP7::p7_state> path) {
	for(deque<BandedHMMP7::p7_state>::const_iterator it = path.begin(); it != path.end(); ++it)
		os << BandedHMMP7::decode(*it);
	return os;
}

BandedHMMP7 BandedHMMP7::build(const MSA& msa, double symfrac,
		const BandedHMMP7Prior& prior, const string& name) {
	if(msa.getMSALen() == 0)
		throw invalid_argument("Empty MSA encountered");
	if(!(symfrac > 0 && symfrac < 1))
		throw invalid_argument("symfrac must between 0 and 1");
	/* determine the bHMM size */
	const unsigned L = msa.getCSLen();
	const unsigned N = msa.getNumSeq();
	unsigned k = 0;
	int cs2ProfileIdx[kMaxProfile + 1]; // MAP index from consensus index -> profile index
	int profile2CSIdx[kMaxCS + 1]; // MAP index from profile index -> consensus index
	for(int i = 0; i <= kMaxProfile; ++i)
		cs2ProfileIdx[i] = 0;
	for(int i = 1; i <= kMaxCS; ++i)
		profile2CSIdx[i] = 0;

	for(unsigned j = 0; j < L; ++j) {
		if(msa.symWFrac(j) >= symfrac)
			profile2CSIdx[++k] = j + 1; /* all index are 1-based */
		cs2ProfileIdx[j+1] = k;
	}
	const int K = k;
	const int csStart = profile2CSIdx[1];
	const int csEnd = profile2CSIdx[K];

	BandedHMMP7 hmm(msa.getName(), K, msa.getAbc()); /* construct an empty hmm */

	/* reset transition and emisison matrices */
	hmm.reset_transition_params();
	hmm.reset_emission_params();
	/* copy the index */
	std::copy(cs2ProfileIdx, cs2ProfileIdx + kMaxProfile, hmm.cs2ProfileIdx);
	std::copy(profile2CSIdx, profile2CSIdx + kMaxCS, hmm.profile2CSIdx);

	/* train the hmm model using observed count, all index are 1-based */
	for(unsigned j = 1; j <= L; ++j) {
		unsigned k = cs2ProfileIdx[j];
		for(unsigned i = 1; i <= N; ++i) {
			int8_t b = msa.encodeAt(i - 1, j - 1);
			double w = msa.getSeqWeight(i - 1); /* use weighted count */
			p7_state sm = determineMatchingState(cs2ProfileIdx, j, b);
			if(sm == P)
				continue; // ignore this base
//			cerr << "j:" << j << " k:" << k << " i:" << i << " sm:" << sm << endl;
			/* update emission frequencies */
			if(sm == M) {
//				cerr << "i:" << i << " j:" << j << " b:" << (int) b << " db:" << hmm.abc->decode(b) << " k:" << k << endl;
				hmm.E_M(b, 0) += w; /* M0 as the COMPO freq */
				hmm.E_M(b, k) += w;
			}
			else if(sm == I) {
//				cerr << "i:" << i << " j:" << j << " b:" << (int) b << " db:" << hmm.abc->decode(b) << " k:" << k << endl;
				hmm.E_I(b, k) += w;
			}
			else { } // no emission

			/* update transition frequencies */
			unsigned jN;
			p7_state smN;
			/* find the next non P loc on this seq */
			for(jN = j + 1; jN <= L; ++jN) {
				int8_t bN = msa.encodeAt(i - 1, jN - 1);
				p7_state smN = determineMatchingState(cs2ProfileIdx, jN, bN);
				if(smN != P)
					break;
			}
			if(!(jN <= L && smN != P)) // no jN found
				continue;
			unsigned kN = cs2ProfileIdx[jN];
			if(sm == I && smN == D || sm == D && smN == I) // no I->D or D->I allowed
				continue;
//			if(sm == D && (j < msa->seqStart(i) + 1 || j > msa->seqEnd(i) + 1)) // 5' and 3' hanging gaps are ignored
//				continue;
			hmm.Tmat[k](sm, smN) += w;
		} // end each seq
	} // end each loc
	/* update B->M1/I0/D1 and MK/IK/DK->E frequencies */
	for(unsigned i = 0; i < N; ++i) {
		double w = msa.getSeqWeight(i);
		int start = msa.seqStart(i);
		int end = msa.seqEnd(i);
		int8_t bStart = msa.encodeAt(i, start);
		p7_state smStart = determineMatchingState(cs2ProfileIdx, start + 1, bStart);
		hmm.Tmat[0](M, smStart) += w;
		int8_t bEnd = msa.encodeAt(i, end);
		p7_state smEnd = determineMatchingState(cs2ProfileIdx, end + 1, bEnd);
		hmm.Tmat[K](smEnd, M) += w;
	}
	hmm.nSeq = msa.getNumSeq();
	hmm.effN = hmm.nSeq;

	/* tune the effN to target mean relative entropy */

	RelativeEntropyTargetFunc entFunc(DEFAULT_ERE, hmm, prior);
	Math::RootFinder rf(entFunc, 0, hmm.nSeq);
	hmm.effN = rf.rootBisection();
	if(::isnan(hmm.effN))
		hmm.effN = hmm.nSeq;
//	cerr << "Final HMM EFFN: " << hmm.effN << endl;
	hmm.scale(hmm.effN / hmm.nSeq);
	hmm.estimateParams(prior);

	/* set bgFreq */
	hmm.hmmBg.setBgFreq(hmm.E_M.col(0));
	hmm.setSpEmissionFreq(hmm.E_M.col(0));

	/* set optional tags */
	char value[128];
	sprintf(value, "%d", msa.getCSLen());
	hmm.setOptTag("MAXL", value);

	hmm.setOptTag("RF", "no");

	hmm.setOptTag("MM", "no");

	hmm.setOptTag("CONS", "yes");

	hmm.setOptTag("CS", "no");

	hmm.setOptTag("MAP", "yes");

	sprintf(value, "%d", hmm.nSeq);
	hmm.setOptTag("NSEQ", value);

	sprintf(value, "%g", hmm.effN);
	hmm.setOptTag("EFFN", value);

	/* set locOptTags */
	hmm.locOptTags["CONS"].resize(K + 1);
	hmm.locOptTags["MAP"].resize(K + 1);
	for(int k = 1; k <= K; ++k) {
		int map = profile2CSIdx[k];
		sprintf(value, "%d", map);
		hmm.setLocOptTag("MAP", value, k);
		char c = msa.CSBaseAt(map);
		int8_t b = msa.getAbc()->encode(c);
		if(msa.wIdentityAt(map) < CONS_THRESHOLD)
			c = ::tolower(c);
		hmm.setLocOptTag("CONS", string() + c, k);
	}

	/* set DATE tag after all done */
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(value, 128, "%c", timeinfo);
	hmm.setOptTag("DATE", value);

	return hmm;
}

string BandedHMMP7::trim(const string& str, const string& whitespace) {
	const string::size_type strBegin = str.find_first_not_of(whitespace);
	if(strBegin == string::npos) // no content
		return "";
	string::size_type strRange = str.find_last_not_of(whitespace) - strBegin + 1;
	return str.substr(strBegin, strRange);
}


/*void EGriceLab::BandedHMMP7::setSequenceMode(enum align_mode mode) {
 switch(mode) {
 case GLOBAL:
 }
 }*/

void BandedHMMP7::setProfileSize(int size) {
	K = size; // set self size
	hmmBg.setSize(size); // set bg size
	init_transition_params();
	init_emission_params();
	init_special_params();
	init_index();
	init_limits();
	enableProfileLocalMode(); // always in profile local alignment mode
	setSpEmissionFreq(hmmBg.getBgEmitPr()); // Use bg emission freq by default
}

void BandedHMMP7::setSequenceMode(enum align_mode mode) {
	switch (mode) {
	case GLOBAL:
		T_SP(N, N) = T_SP(C, C) = 0;
		break;
	case LOCAL:
		T_SP(N, N) = T_SP(C, C) = hmmBg.getBgTransPr();
		break;
	case NGCL:
		T_SP(N, N) = 0;
		T_SP(C, C) = hmmBg.getBgTransPr();
		break;
	case CGNL:
		T_SP(N, N) = hmmBg.getBgTransPr();
		T_SP(C, C) = 0;
		break;
	default:
		break; // do nothing
	}
	T_SP(N, B) = 1.0 - T_SP(N, N);
	T_SP(E, C) = 1.0; // always exit from E->C
	T_SP_cost = -T_SP.array().log(); // Eigen3 handle array to matrix assignment automatically
}

void BandedHMMP7::setSpEmissionFreq(const Vector4d& freq) {
	E_SP.col(N) = E_SP.col(C) = freq / freq.sum(); // re-do normalization, even if already done
	E_SP.col(B) = E_SP.col(E) = Vector4d::Zero(); // no emission for state B and E
	E_SP_cost = -E_SP.array().log();
}

void BandedHMMP7::init_transition_params() {
	if (transInit) // already initialized
		return;
	/*
	 * state 0 serves as the B state
	 */
	/* initiate transition matrices */
	Tmat = vector<Matrix3d>(K + 1);
	Tmat_cost = vector<Matrix3d>(K + 1);
	transInit = true;
}

void BandedHMMP7::init_emission_params() {
	if (emisInit) // already initialized
		return;
	/* state 0 serves as B state */
	E_M = E_M_cost = Matrix4Xd(4, K + 1);
	E_I = E_I_cost = Matrix4Xd(4, K + 1);
	emisInit = true;
}

void BandedHMMP7::init_special_params() {
	if (spInit) // already initialized
		return;
	entryPr = entryPr_cost = VectorXd(K + 1);
	exitPr = exitPr_cost = VectorXd(K + 1);
	E_SP = E_SP_cost = Matrix4Xd(4, kNS);
	T_SP = T_SP_cost = MatrixXd(kNS, kNS);
	spInit = true;
}

void BandedHMMP7::reset_transition_params() {
	/*
	 * state 0 serves as the B state
	 */
	for(int k = 0; k <= K; ++k) {
		Tmat[k].setZero();
		Tmat_cost[k].fill(inf);
	}
}

void BandedHMMP7::reset_emission_params() {
	/* state 0 serves as B state */
	E_M.setZero();
	E_I.setZero();
	E_M_cost.fill(inf);
	E_I_cost.fill(inf);
}

/*void BandedHMMP7::normalize_transition_params() {

	 * state 0 serves as the B state

	for(int k = 0; k <= K; ++k) {
		for(int i = 0; i < BandedHMMP7::kNM; ++i) {
			double C = Tmat[k].row(i).sum();
			double pseudoC = BandedHMMP7::pseudoCount(C);
			Tmat[k].row(i).array() += pseudoC / Tmat[k].cols();
			Tmat[k].row(i) /= C + pseudoC;
		}
		Tmat_cost[k] = Tmat[k].array().log();
	}

}*/

/*void BandedHMMP7::normalize_emission_params() {
	 state 0 serves as B state
	for(int k = 0; k <= K; ++k) {
		double emC = E_M.col(k).sum();
		double eiC = E_I.col(k).sum();
		if(emC > 0) {
			double emPseudo = BandedHMMP7::pseudoCount(emC);
			E_M.col(k).array() += emPseudo / E_M.rows();
			E_M.col(k) /= emC + emPseudo;
		}
		else {
			E_M.col(k).fill(1.0 / E_M.rows());  Nothing observed, use constants
		}

		if(eiC > 0) {
			double eiPseudo = BandedHMMP7::pseudoCount(eiC);
			E_I.col(k).array() += eiPseudo / E_I.rows();
			E_I.col(k) /= eiC + eiPseudo;
		}
		else {
			E_I.col(k).fill(1.0 / E_I.rows());  Nothing observed, use constants
		}
	}
	E_M_cost = -E_M.array().log();
	E_I_cost = -E_I.array().log();
}*/

void BandedHMMP7::init_limits() {
	if (limitInit) // already initialized
		return;
/*	insBeforeLimit = insAfterLimit = VectorXi::Constant(K + 1, static_cast<int> (kMinProfile));
	delBeforeLimit = delAfterLimit = VectorXi::Constant(K + 1, static_cast<int> (kMinProfile));*/

	gapBeforeLimit = gapAfterLimit = VectorXi(K + 1);
	//delBeforeLimit = delAfterLimit = VectorXi(K + 1);
	for(VectorXi::Index j = 1; j <= K; ++j) {
		gapBeforeLimit(j) = j * kMinGapFrac;
		gapAfterLimit(j) = (K - j) * kMinGapFrac;
	}
	limitInit = true;
}

void BandedHMMP7::init_index() {
	if(indexInit) // already initialized
		return;
	for(int i = 0; i <= kMaxProfile; ++i)
		cs2ProfileIdx[i] = 0;
	for(int i = 1; i <= kMaxCS; ++i)
		profile2CSIdx[i] = 0;
}

void BandedHMMP7::enableProfileLocalMode() {
	/* set entering costs */
	entryPr(0) = 0; // B->B not allowed
	entryPr.segment(1, K).setConstant(1 - hmmBg.getBgTransPr()); /* B->M1..MK equal cost */

	/* set exiting costs */
	exitPr(0) = 0; // B->E not allowed
	exitPr.segment(1, K).setConstant(1 - hmmBg.getBgTransPr()); /* M1..MK ->E equal cost */

	/* set log versions */
	entryPr_cost = -entryPr.array().log();
	exitPr_cost = -exitPr.array().log();
}

void BandedHMMP7::adjustProfileLocalMode() {
	/* adjust entering costs */
	entryPr(0) = 0; // B->B not allowed
	entryPr.segment(1, K).setConstant(Tmat[0](M, M)); /* B->M1..MK equal cost */

	/* set exiting costs */
	exitPr(0) = 0; // B->E not allowed
	exitPr.segment(1, K).setConstant(Tmat[K](M, M)); /* M1..MK ->E equal cost */

	/* set log versions */
	entryPr_cost = -entryPr.array().log();
	exitPr_cost = -exitPr.array().log();
}

BandedHMMP7::ViterbiScores EGriceLab::BandedHMMP7::initViterbiScores(
		const PrimarySeq& seq) const {
	assert(abc == seq.getDegenAlphabet());
	ViterbiScores vs(seq);

	return resetViterbiScores(vs, seq);
/*	resetForwardBackwardScores(bHmm);
	resetTraceScores(bHmm);*/
}

BandedHMMP7::ViterbiScores& EGriceLab::BandedHMMP7::resetViterbiScores(ViterbiScores& vs,
		const PrimarySeq& seq) const {
	vs.seq = &seq; // swap the pointer value but leave the pointing object untouched
	const int L = vs.seq->length();
	/* Use -Inf as initial values for Viterbi matrices */
/*	vs.DP_M = vs.DP_I = vs.DP_D = vs.S = MatrixXd::Constant(L + 1, K + 1,
			-numeric_limits<double>::infinity());*/
	/* set dimension */
	vs.DP_M.resize(L + 1, K + 1);
	vs.DP_I.resize(L + 1, K + 1);
	vs.DP_D.resize(L + 1, K + 1);

	/* set elements */
	vs.DP_M.setConstant(inf);
	vs.DP_I.setConstant(inf);
	vs.DP_D.setConstant(inf);
	return vs;
}

BandedHMMP7::ViterbiAlignPath EGriceLab::BandedHMMP7::initViterbiAlignPath(
		int L) const {
	ViterbiAlignPath vpath;
	vpath.K = K; // K is fixed between different seq
	vpath.N = 0;
	return resetViterbiAlignPath(vpath, L);
}

BandedHMMP7::ViterbiAlignPath& EGriceLab::BandedHMMP7::resetViterbiAlignPath(ViterbiAlignPath& vpath,
		int L) const {
	vpath.L = L;
	vpath.N = 0;
	/* reset values */
	vpath.start.clear();
	vpath.end.clear();
	vpath.from.clear();
	vpath.to.clear();
	vpath.numIns.clear();
	vpath.numDel.clear();
	vpath.alnStart = vpath.alnEnd = 0;
	vpath.alnFrom = vpath.alnTo = 0;
	//vpath.path.clear();
	vpath.alnPath.clear();
	/* vpath.TRACE.resize(vpath.L + 1, vpath.K + 1); */
	//vpath.TRACE.resize(L + 1, vpath.K + 1);
	//vpath.TRACE.setConstant(-1);
	return vpath;
}

void BandedHMMP7::calcViterbiScores(
		ViterbiScores& vs) const {
	assert(vs.seq != NULL);
	assert(wingRetracted); // make sure wing is retracted
	const int L = vs.seq->length();
	/* Initialize the M(,0), the B state */
	for (int i = 0; i <= L; ++i)
		vs.DP_M(i, 0) = i <= 1 ? 0 /* no N-N loop */ : T_SP_cost(N, N) * (i - 1); /* N->N loops */
	vs.DP_M.col(0).array() += T_SP_cost(N, B); /* N->B */

	/* calculate I0, a special case */
	for (int i = 1; i <= L; ++i)
		vs.DP_I(i, 0) = E_I_cost(vs.seq->encodeAt(i - 1), 0) + std::min(
					static_cast<double> (vs.DP_M(i - 1, 0) + Tmat_cost[0](M, I)), // B->I0
					static_cast<double> (vs.DP_I(i - 1, 0) + Tmat_cost[0](I, I))); // I0->I0
	/* Full Dynamic-Programming at row-first order */
	for (int j = 1; j <= K; ++j) {
		for (int i = 1; i <= L; ++i) {
			vs.DP_M(i, j) = E_M_cost(vs.seq->encodeAt(i-1), j) + EGriceLab::BandedHMMP7::min(
					static_cast<double> (vs.DP_M(i - 1, 0) + entryPr_cost(j)), // from the B state
					static_cast<double> (vs.DP_M(i - 1, j - 1) + Tmat_cost[j-1](M, M)), // from Mi-1,j-1
					static_cast<double> (vs.DP_I(i - 1, j - 1) + Tmat_cost[j-1](I, M)), // from Ii-1,j-1
					static_cast<double> (vs.DP_D(i - 1, j - 1) + Tmat_cost[j-1](D, M))); // from Di-1,j-1
			vs.DP_I(i, j) = E_I_cost(vs.seq->encodeAt(i - 1), j) + std::min(
							static_cast<double> (vs.DP_M(i - 1, j) + Tmat_cost[j](M, I)), // from Mi-1,j
							static_cast<double> (vs.DP_I(i - 1, j) + Tmat_cost[j](I, I))); // from Ii-1,j
			vs.DP_D(i, j) = std::min(
					static_cast<double> (vs.DP_M(i, j - 1) + Tmat_cost[j-1](M, D)), // from Mi,j-1
					static_cast<double> (vs.DP_D(i, j - 1) + Tmat_cost[j-1](D, D))); // from Di,j-1
		}
	}
	vs.S.resize(L + 1, K + 2); // column K+1 represent the exit from IK state
	vs.S.leftCols(K + 1) = vs.DP_M;; // 0..K columns copied from the calculated DP_M
	vs.S.col(K + 1) = vs.DP_I.col(K);
	vs.S.col(0).setConstant(inf); /* S(,0) is not useful */
	// add M-E exit costs
	vs.S.leftCols(K + 1).rowwise() += exitPr_cost.transpose();
	vs.S.col(K + 1).array() += Tmat_cost[K](I, M); // IK->E
	vs.S.array() += T_SP_cost(E, C); // add E->C transition
	for (int i = 1; i < L; ++i) // S(L,) doesn't have a C-> loop
		vs.S.row(i).array() += T_SP_cost(C, C) * (L - i); // add L-i C->C circles
}

void BandedHMMP7::calcViterbiScores(
		ViterbiScores& vs, const ViterbiAlignPath& vpath) const {
	assert(vs.seq != NULL);
	if(vpath.N == 0) // no known path provided, reduce to full Viterbi Dynamic-Programming
		return calcViterbiScores(vs);
	const int L = vs.seq->length();
	assert(wingRetracted);
	assert(K + 1 == vs.DP_M.cols() && K == vpath.K);
	assert(L == vpath.L);
	if (!isValidAlignPath(vpath)) {
		std::cerr << "Incorrect known seed alignment paths provided!" << endl;
		abort();
	}
/*	cerr << "path valid" << endl;*/

	/* vs.DP_M(0, 0) = 0; */
/*	cerr << "DP matrices initialized" << endl;*/
	/* Initialize the M(,0), the B state */
	for (int i = 1; i < L; i++)
		vs.DP_M(i, 0) = i == 1 ? 0 /* no N->N loop */ : T_SP_cost(N, N) * (i - 1); /* N->N loops */
	vs.DP_M.col(0).array() += T_SP_cost(N, B); /* N->B */
//	 cerr << "M(,0) initialized" << vs.DP_M.col(0).transpose() << endl;
	/* calculate I0, a special case */
	for (int i = 1; i <= L; ++i)
		vs.DP_I(i, 0) = E_I_cost(vs.seq->encodeAt(i - 1), 0) + std::min(
					static_cast<double> (vs.DP_M(i - 1, 0) + Tmat_cost[0](M, I)), // B->I0
					static_cast<double> (vs.DP_I(i - 1, 0) + Tmat_cost[0](I, I))); // I0->I0
	/* process each known path upstream and themselves */
	for(vector<int>::size_type n = 0; n < vpath.N; ++n) {
		/* Determine banded boundaries */
		int upQLen = n == 0 ? vpath.from[n] - 1 : vpath.from[n] - vpath.to[n - 1];
		int up_start = n == 0 ? vpath.start[n] - upQLen * (1 + kMinGapFrac) : vpath.end[n - 1];
		if (up_start < 1)
			up_start = 1;
		int up_from = n == 0 ? vpath.from[n] - upQLen * (1 + kMinGapFrac) : vpath.to[n - 1];
		if (up_from < 1)
			up_from = 1;
		//cerr << "upQLen:" << upQLen << endl;
//		cerr << "up_start:" << up_start << " up_end:" << vpath.start[n] << endl;
//		cerr << "up_from:" << up_from << " up_to:" << vpath.from[n] << endl;

		/* Dynamic programming of upstream of this known path at row-first order */
		for (int j = up_start; j <= vpath.start[n]; ++j) {
			for (int i = up_from; i <= vpath.from[n]; ++i) {
				vs.DP_M(i, j) = E_M_cost(vs.seq->encodeAt(i - 1), j)
						+ EGriceLab::BandedHMMP7::min(
								static_cast<double>(vs.DP_M(i, 0) + entryPr_cost(j)), // from B state
								static_cast<double>(vs.DP_M(i - 1, j - 1) + Tmat_cost[j-1](M, M)), // from Mi-1,j-1
								static_cast<double>(vs.DP_I(i - 1, j - 1) + Tmat_cost[j-1](I, M)), // from Ii-1,j-1
								static_cast<double>(vs.DP_D(i - 1, j - 1) + Tmat_cost[j-1](D, M))); // from Di-1,j-1
				vs.DP_I(i, j) = E_I_cost(vs.seq->encodeAt(i - 1), j)
								+ std::min(
										static_cast<double>(vs.DP_M(i - 1, j) + Tmat_cost[j](M, I)), // from Mi-1,j
										static_cast<double>(vs.DP_I(i - 1, j) + Tmat_cost[j](I, I))); // from Ii-1,j
				vs.DP_D(i, j) =	std::min(
						static_cast<double>(vs.DP_M(i, j - 1) + Tmat_cost[j-1](M, D)), // from Mi,j-1
						static_cast<double>(vs.DP_D(i, j - 1) + Tmat_cost[j-1](D, D))); // from Di,j-1
			}
		}
		/* Fill the score of the known alignment path */
		for (int j = vpath.start[n]; j <= vpath.end[n]; ++j) {
			for(int i = vpath.from[n]; i <= vpath.to[n]; ++i) {
				int dist = diagnalDist(i, j, vpath.from[n], vpath.start[n]);
				if(!(dist <= vpath.numIns[n] && dist >= -vpath.numDel[n]))
					continue;
				vs.DP_M(i, j) = E_M_cost(vs.seq->encodeAt(i - 1), j)
						+ EGriceLab::BandedHMMP7::min(
								static_cast<double>(vs.DP_M(i, 0) + entryPr_cost(j)), // from B state
								static_cast<double>(vs.DP_M(i - 1, j - 1) + Tmat_cost[j-1](M, M)), // from Mi-1,j-1
								static_cast<double>(vs.DP_I(i - 1, j - 1) + Tmat_cost[j-1](I, M)), // from Ii-1,j-1
								static_cast<double>(vs.DP_D(i - 1, j - 1) + Tmat_cost[j-1](D, M))); // from Di-1,j-1
				vs.DP_I(i, j) = E_I_cost(vs.seq->encodeAt(i - 1), j)
						+ std::min(
								static_cast<double>(vs.DP_M(i - 1, j) + Tmat_cost[j](M, I)), // from Mi-1,j
								static_cast<double>(vs.DP_I(i - 1, j) + Tmat_cost[j](I, I))); // from Ii-1,j
				vs.DP_D(i, j) = std::min(
						static_cast<double>(vs.DP_M(i, j - 1) + Tmat_cost[j-1](M, D)), // from Mi,j-1
						static_cast<double>(vs.DP_D(i, j - 1) + Tmat_cost[j-1](D, D))); // from Di,j-1
			}
		}
		// assert(i == vpath.to + 1 && j == vpath.end + 1);
	} /* end of each known path segment */
//	cerr << "known path aligned" << endl;
	/* Dynamic programming of the remaining downstream of the known paths, if any */
	int last_end = vpath.end[vpath.N - 1];
	int last_to = vpath.to[vpath.N - 1];
	int downQLen = L - last_to;
	int down_end = last_end + downQLen * (1 + kMinGapFrac);
	int down_to = last_to + downQLen * (1 + kMinGapFrac);
	if(down_end > K)
		down_end = K;
	if(down_to > L)
		down_to = L;

	for (int j = last_end; j <= down_end; ++j) {
		for (int i = last_to; i <= down_to; ++i) {
			vs.DP_M(i, j) = E_M_cost(vs.seq->encodeAt(i - 1), j) +
					EGriceLab::BandedHMMP7::min(
							// from Mi,0, the B state is not possible
							static_cast<double>(vs.DP_M(i - 1, j - 1) + Tmat_cost[j-1](M, M)), // from Mi-1,j-1
							static_cast<double>(vs.DP_I(i - 1, j - 1) + Tmat_cost[j-1](I, M)), // from Ii-1,j-1
							static_cast<double>(vs.DP_D(i - 1, j - 1) + Tmat_cost[j-1](D, M))); // from Di-1,j-1
			vs.DP_I(i, j) = E_I_cost(vs.seq->encodeAt(i - 1), j) +
					std::min(
							static_cast<double>(vs.DP_M(i - 1, j) + Tmat_cost[j](M, I)), // from Mi-1,j
							static_cast<double>(vs.DP_I(i - 1, j) + Tmat_cost[j](I, I))); // from Ii-1,j
			vs.DP_D(i, j) = std::min(
							static_cast<double>(vs.DP_M(i, j - 1) + Tmat_cost[j-1](M, D)), // from Mi,j-1
							static_cast<double>(vs.DP_D(i, j - 1) + Tmat_cost[j-1](D, D))); // from Di,j-1
		}
	}
//	cerr << "downstream done" << endl;
	vs.S.resize(L + 1, K + 2); // column K+1 represent the exit from IK state
	vs.S.leftCols(K + 1) = vs.DP_M;; // 0..K columns copied from the calculated DP_M
	vs.S.col(K + 1) = vs.DP_I.col(K);
	//vs.S.col(K + 1).setConstant(inf);
	vs.S.col(0).setConstant(inf);
	// add M-E exit costs
	vs.S.leftCols(K + 1).rowwise() += exitPr_cost.transpose();
	vs.S.col(K + 1).array() += Tmat_cost[K](I, M); // IK->E
	vs.S.array() += T_SP_cost(E, C); // add E->C transition
	for (int i = 1; i < L; ++i)
		vs.S.row(i).array() += T_SP_cost(C, C) * (L - i); // add L-i C->C circles
}

void BandedHMMP7::addKnownAlignPath(ViterbiAlignPath& vpath,
		const CSLoc& csLoc, int csFrom, int csTo) const {
	using std::string;
//	cerr << "csStart:" << csLoc.start << " csEnd:" << csLoc.end << " csFrom:" << csFrom << " csTo:" << csTo <<
//			" CS:" << csLoc.CS << endl;
	assert(csLoc.start >= 1 && csLoc.start <= csLoc.end
			&& csFrom >= 1 && csFrom <= csTo
			&& csLoc.CS.length() > csLoc.end - csLoc.start && csLoc.CS.length() > csTo - csFrom);
	/* calculate profile start, end and path */
	int start = 0;
	int end = 0;
	int from = 0;
	int to = 0;
	int numIns = 0;
	int numDel = 0;

	int i = csFrom;
	int j = csLoc.start;
	for(string::const_iterator it = csLoc.CS.begin(); it != csLoc.CS.end(); ++it) {
		int k = getProfileLoc(j); // position on profile
//		cerr << "i:" << i << " j:" << j << " k:" << k << endl;
//		cerr << "vpath.L:" << vpath.L << " vpath.K:" << vpath.K << endl;

		bool nonGap = abc->isSymbol(*it);

		if(from == 0 && nonGap)
			from = i;
		if(nonGap)
			to = i;
		if(k != 0) { // a non-D loc on profile
			if(start == 0) // first time a non-D loc on profile
				start = k;
			end = k; // keep updating
			if(!nonGap) // a deletion
				numDel++;
		}
		else { // a D loc on profile
			if(nonGap) // an insertion
				numIns++;
		}
		j++; // update j
		if(nonGap)
			i++; // update i
	}
//	cerr << "vpath.path:" << vpath.alnPath << endl;
//	cerr << "i:" << i << " j:" << j << " csTo:" << csTo << " csEnd:" << csLoc.end << endl;
	assert(i == csTo + 1 && j == csLoc.end + 1);
	if(from != 0 && to != 0 && start != 0 && end != 0) {
		vpath.start.push_back(start);
		vpath.end.push_back(end);
		vpath.from.push_back(from);
		vpath.to.push_back(to);
		vpath.numIns.push_back(numIns);
		vpath.numDel.push_back(numDel);
		vpath.N++;
	}
}

double BandedHMMP7::buildViterbiTrace(const ViterbiScores& vs, ViterbiAlignPath& vpath) const {
	assert(vs.seq != NULL);
	assert(vs.S.rows() == vpath.L + 1 && vs.S.cols() == K + 2);
	MatrixXd::Index minRow, minCol;
	double minScore = vs.S.minCoeff(&minRow, &minCol);
	if(minScore == inf)
		return minScore; // min score not found, do nothing
	/* do trace back in the vScore matrix */
	char s = minCol <= K ? 'M' : 'I'; // exiting state either M1..K or IK
	int i = minRow;
	int j = minCol <= K ? minCol : K;
	vpath.alnStart = minCol <= K ? minCol : K;
	vpath.alnEnd = minCol <= K ? minCol : K;
	vpath.alnFrom = vpath.alnTo = minRow;
	//cerr << "id:" << vs.seq->getId() << " minRow:" << minRow << " minCol:" << minCol << " minScore:" << minScore << " minEnd:" << getCSLoc(minCol) << endl;

	vpath.alnPath.push_back('E'); // ends with E
	while(i >= 0 && j >= 0) {
		//vpath.TRACE(i, j) = s;
		vpath.alnPath.push_back(s);
		if(j == 0) {
		//	cerr << "i:" << i << " j:" << j << " s:" << decode(B) << ":" << static_cast<double> (vs.DP_M(i, j)) << endl;
			break;
		}
		// update the status
		if(s == 'M') {
			//cerr << "i:" << i << " j:" << j << " s:" << s << ":" << static_cast<double> (vs.DP_M(i, j));
			vpath.alnStart--;
			vpath.alnFrom--;
			s = BandedHMMP7::whichMin(
					static_cast<double> (vs.DP_M(i, 0) + entryPr_cost(j)), /* from B-state */
					static_cast<double> (vs.DP_M(i - 1, j - 1) + Tmat_cost[j-1](M, M)), /* from M(i-1,j-1) */
					static_cast<double> (vs.DP_I(i - 1, j - 1) + Tmat_cost[j-1](I, M)), /* from I(i-1,j-1) */
					static_cast<double> (vs.DP_D(i - 1, j - 1) + Tmat_cost[j-1](D, M))); /* from D(i-1,j-1) */
			//cerr << "->" << s << endl;
			if(s == 'B')
				j = 0; // jump to B state
			else { // M, I or D state
				i--;
				j--;
			}
		}
		else if(s == 'I') {
			//cerr << "i:" << i << " j:" << j << " s:" << s << ":" << static_cast<double> (vs.DP_I(i, j));
			vpath.alnFrom--;
			s = BandedHMMP7::whichMin(
					static_cast<double> (vs.DP_M(i - 1, j) + Tmat_cost[j](M, I)), /* from M(i-1,j) */
					static_cast<double> (vs.DP_I(i - 1, j) + Tmat_cost[j](I, I)), /* from I(i-1,j) */
					"MI");
			//cerr << "->" << s << endl;
			i--;
		}
		else if(s == 'D') {
			vpath.alnStart--;
			//cerr << "i:" << i << " j:" << j << " s:" << s << ":" << static_cast<double> (vs.DP_D(i, j));
			s = BandedHMMP7::whichMin(
					static_cast<double> (vs.DP_M(i, j - 1) + Tmat_cost[j-1](M, D)), /* from M(i,j-1) */
					static_cast<double> (vs.DP_D(i, j - 1) + Tmat_cost[j-1](D, D)), /* from D(i,j-1) */
					"MD");
			//cerr << "->" << s << endl;
			j--;
		}
		else { // B state found
			//cerr << "i:" << i << " j:" << j << " s:" << s << ":" << static_cast<double> (vs.DP_M(i, j)) << endl;
			break;
		}
		//cerr << "vpath.alnStart:" << vpath.alnStart + 1 << " vpath.alnEnd:" << vpath.alnEnd << endl;
	} /* end of while */
	vpath.alnStart++; // 1-based
	vpath.alnFrom++; // 1-based
	if(vpath.alnPath[vpath.alnPath.length() - 1] != 'B')
		vpath.alnPath.push_back('B');
	reverse(vpath.alnPath.begin(), vpath.alnPath.end()); // reverse the alnPath string
	return minScore;
}

std::string BandedHMMP7::buildGlobalAlignSeq(const ViterbiScores& vs,
		const ViterbiAlignPath& vpath) const {
	assert(vs.seq != NULL);
	string aSeq;
	int profileNLen = vpath.alnStart - 1;
	int profileCLen = vpath.K - vpath.alnEnd;
	int seqNLen = vpath.alnFrom - 1;
	int seqCLen = vpath.L - vpath.alnTo;

	/* place N state residues to the beginning and ending of the unmatched up-stream of the profile */
	int NHalf = std::max(profileNLen / 2, seqNLen / 2);
	int CHalf = std::max(profileCLen / 2, seqCLen / 2);
	for(int i = 0; i < NHalf; ++i) /* put half the N residues at the beginning of the N' */
		aSeq.push_back(vs.seq->charAt(i - 1));
	if(profileNLen > 2 * NHalf)
		aSeq.append(profileNLen - 2 * NHalf, '-'); /* add the middle gaps, if any */
	for(int i = profileNLen - NHalf; i < profileNLen; ++i) /* put half the N residues at the beginning */
		aSeq.push_back(vs.seq->charAt(i - 1)); /* put half the N residues at the end of the N' */

	/* place aligned residues */
	int i = vpath.alnFrom;
	for(string::const_iterator it = vpath.alnPath.begin(); it != vpath.alnPath.end(); ++it) {
		switch(*it) {
		case 'M': case 'I': /* both case consumes residues on seq */
			aSeq.push_back(vs.seq->charAt(i - 1));
			i++;
			break;
		case 'D': /* gap relative to profile */
			aSeq.push_back('-');
			break;
		default:
			break; // do nothing
		}
	}
	for(int i = vpath.alnTo; i < vpath.alnTo + CHalf; ++i) /* put half the N residues at the beginning of the N' */
		aSeq.push_back(vs.seq->charAt(i - 1));
	aSeq.append(profileCLen - 2 * CHalf, '-'); /* add the middle gaps, if any */
	for(int i = vpath.K - NHalf; i < vpath.K; ++i) /* put half the N residues at the beginning */
		aSeq.push_back(vs.seq->charAt(i - 1)); /* put half the N residues at the end of the N' */
	return aSeq;
}

DigitalSeq BandedHMMP7::buildGlobalAlignDS(const ViterbiScores& vs,
		const ViterbiAlignPath& vpath) const {
	assert(vs.seq != NULL);
	DigitalSeq aSeq;
	int profileNLen = vpath.alnStart - 1;
	int profileCLen = vpath.K - vpath.alnEnd;
	int seqNLen = vpath.alnFrom - 1;
	int seqCLen = vpath.L - vpath.alnTo;

	/* place N state residues to the beginning and ending of the unmatched up-stream of the profile */
	int NHalf = std::max(profileNLen / 2, seqNLen / 2);
	int CHalf = std::max(profileCLen / 2, seqCLen / 2);
	for(int i = 0; i < NHalf; ++i) /* put half the N residues at the beginning of the N' */
		aSeq.push_back(vs.seq->encodeAt(i - 1));
	if(profileNLen > 2 * NHalf)
		aSeq.append(profileNLen - 2 * NHalf, DegenAlphabet::GAP_SYM); /* add the middle gaps, if any */
	for(int i = profileNLen - NHalf; i < profileNLen; ++i) /* put half the N residues at the beginning */
		aSeq.push_back(vs.seq->encodeAt(i - 1)); /* put half the N residues at the end of the N' */

	/* place aligned residues */
	int i = vpath.alnFrom;
	for(string::const_iterator it = vpath.alnPath.begin(); it != vpath.alnPath.end(); ++it) {
		switch(*it) {
		case 'M': case 'I': /* both case consumes residues on seq */
			aSeq.push_back(vs.seq->encodeAt(i - 1));
			i++;
			break;
		case 'D': /* gap relative to profile */
			aSeq.push_back(DegenAlphabet::GAP_SYM);
			break;
		default:
			break; // do nothing
		}
	}
	for(int i = vpath.alnTo; i < vpath.alnTo + CHalf; ++i) /* put half the N residues at the beginning of the N' */
		aSeq.push_back(vs.seq->encodeAt(i - 1));
	aSeq.append(profileCLen - 2 * CHalf, DegenAlphabet::GAP_SYM); /* add the middle gaps, if any */
	for(int i = vpath.K - NHalf; i < vpath.K; ++i) /* put half the N residues at the beginning */
		aSeq.push_back(vs.seq->encodeAt(i - 1)); /* put half the N residues at the end of the N' */
	return aSeq;
}

void BandedHMMP7::wingRetract() {
	if(wingRetracted) // already wing-retracted
		return;
	/* retract entering costs */
	/* increase the B->Mj cost by adding the chain B->D1->D2->...->Dj-1->Mj */
	for(MatrixXd::Index j = 2; j <= K; ++j) {
		double cost = 0; // additional retract cost in log-scale
		cost += Tmat_cost[0](M, D); // B->D1 (M0->D1)
		for(MatrixXd::Index i = 1; i < j - 1; ++i)
			cost += Tmat_cost[i](D, D); // Di->Di+1
		cost += Tmat_cost[j-1](D, M); // Dj-1->Mj
		assert(cost > 0);
		entryPr(j) += ::exp(-cost); // retract B->D1->D2...Dj-1->Mj to B->Mj
		if(entryPr(j) > 1)
			entryPr(j) = 1;
	}
	/* retract exiting costs */
	/* increase the Mj->E cost by adding the chain Mj->Dj+1->Dj+2->...->DK->E */
	for(MatrixXd::Index i = 1; i <= K - 1; ++i) {
		double cost = 0; // additional retract cost in log-scale
		cost += Tmat_cost[i](M, D); // Mj -> Di+1
		for(MatrixXd::Index j = i + 1; j < K; ++j)
			cost += Tmat_cost[j](D, D); // Dj->Dj+1
		cost += Tmat_cost[K](D, M); // DK -> E (DK->MK+1)
		assert(cost > 0);
		exitPr(i) += ::exp(-cost); // retract Mj->Dj+1->Dj+2...->DK->E to Mj->E
		if(exitPr(i) > 1)
			exitPr(i) = 1;
	}
	/* set transition matrices */
	/* reset log transition matrices */
//	cerr << "entry before retract: " << entryPr_cost.transpose() << endl;
	entryPr_cost = -entryPr.array().log();
	exitPr_cost = -exitPr.array().log();
//	cerr << "entry after retract: " << entryPr_cost.transpose() << endl;

	wingRetracted = true;
}

bool BandedHMMP7::isValidAlignPath(const ViterbiAlignPath& vpath) const {
	bool flag = true;
	for(int n = 0; n < vpath.N; ++n)
		if(!(vpath.start[n] >= 1 && vpath.start[n] <= vpath.end[n] && vpath.end[n] <= vpath.K
			&& vpath.from[n] >= 1 && vpath.from[n] <= vpath.to[n] && vpath.to[n] <= vpath.L) && // this path is valid
				(n == 0 || (vpath.start[n] > vpath.end[n - 1] && vpath.from[n] > vpath.to[n - 1]))) // path is properly ordered
			return false;
	return true;
}

double RelativeEntropyTargetFunc::operator()(double x) {
	BandedHMMP7 hmm2(hmm); // use a copy so original hmm won't be affected

	if(x > hmm2.effN) // do not scale up
		return 0;

	hmm2.effN = x;
	hmm2.scale(hmm2.effN / hmm2.nSeq);
	hmm2.estimateParams(prior);
	double relEnt = hmm2.meanRelativeEntropy();
//	cerr << "current effN: " << x << " ere: " << relEnt << endl;
//	return hmm.meanRelativeEntropy() - ere;
	return relEnt - ere;
}

} /* namespace EGriceLab */
