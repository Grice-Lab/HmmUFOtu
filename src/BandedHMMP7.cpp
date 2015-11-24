/*
 * BandedHMMP7.cpp
 *
 *  Created on: May 13, 2015
 *      Author: zhengqi
 */

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include "BandedHMMP7.h"

namespace EGriceLab {
using namespace std;
using namespace Eigen;

const int BandedHMMP7::kNM = 3; // number of matching states
const int BandedHMMP7::kNSP = 4; // number of special states
const int BandedHMMP7::kNS =  kNM + kNSP; // number of total states
/*const int BandedHMMP7::kMaxProfile = 10000; // up-to 10K 16S rRNA profile*/
string HMM_TAG =
		"HMM\t\tA\tC\tG\tT\n\t\tm->m\tm->i\tm->d\ti->m\ti->i\td->m\td->d";
const float BandedHMMP7::inf = std::numeric_limits<float>::infinity();

const float BandedHMMP7::infV = -inf;

const float BandedHMMP7::kMaxGapFrac = 0.2; // default at class level

ostream& operator<<(ostream& os, const deque<BandedHMMP7::p7_state> path) {
	for(deque<BandedHMMP7::p7_state>::const_iterator it = path.begin(); it != path.end(); ++it)
		os << BandedHMMP7::decode(*it);
	return os;
}

/* non-member friend functions */
istream& operator>>(istream& is, BandedHMMP7& hmm) {
	string line;
	int k = 0; // pos on the profile
	while (getline(is, line)) {
		if (line == "//") {/* end of profile */
			// reverse the signs
			hmm.E_I_log = -hmm.E_I_log;

			hmm.T_MM_log = -hmm.T_MM_log;
			hmm.T_MI_log = -hmm.T_MI_log;
			hmm.T_MD_log = -hmm.T_MD_log;
			hmm.T_IM_log = -hmm.T_IM_log;
			hmm.T_II_log = -hmm.T_II_log;
			hmm.T_DM_log = -hmm.T_DM_log;
			hmm.T_DD_log = -hmm.T_DD_log;

			// set the non-log matrices
			hmm.E_M = hmm.E_M_log.array().exp();
			hmm.E_I = hmm.E_I_log.array().exp();
			hmm.T_MM = hmm.T_MM_log.array().exp();
			hmm.T_MI = hmm.T_MI_log.array().exp();
			hmm.T_MD = hmm.T_MD_log.array().exp();
			hmm.T_IM = hmm.T_IM_log.array().exp();
			hmm.T_II = hmm.T_II_log.array().exp();
			hmm.T_DM = hmm.T_DM_log.array().exp();
			hmm.T_DD = hmm.T_DD_log.array().exp();
			hmm.wingRetract();
			return is;
		}
		istringstream iss(line); // detail parse this line
		string tag; /* header tag names and values */
		string tmp;
		if (!isspace(line[0])) { /* header section starts with non-empty characters */
			iss >> tag;
			if (tag.substr(0, 6) == "HMMER3") { // do not override our version, check minor version
				if(tag.length() < 8 || tag[7] < 'f') {
					cerr << "Obsolete HMM file version: " << tag << ", must be HMMER3/f [3.1] or higher" << endl;
					abort();
				}
			}
			else if (tag == "NAME") {
				string name;
				iss >> name;
				hmm.setName(name);
				hmm.setHeaderOpt(tag, name);
			} else if (tag == "LENG") {
				string size;
				iss >> size;
				hmm.setProfileSize(atoi(size.c_str()));
				hmm.setHeaderOpt(tag, size);
			} else if (tag == "ALPH") {
				string abc;
				iss >> abc;
				if (abc != "DNA")
					throw invalid_argument(
							"Not allowed alphabet '" + abc
									+ "' in the HMM input file! Must be DNA");
				assert(hmm.nuclAbc->getSymbol() == "ACGT");
				hmm.setHeaderOpt(tag, "DNA");
			} else if (tag == "STATS") {
				string mode;
				string distrib;
				iss >> mode >> distrib;
				tag += " " + mode + " " + distrib; // use STATS + mode + distribution as the new tag name
				string val;
				getline(iss, val);
				hmm.setHeaderOpt(tag, BandedHMMP7::trim(val));
			} else if(tag == "HMM") { /* HMM TAG */
				string tmp;
				getline(is, tmp); /* ignore the next line too */
			}
			else { /* non-mandotory lines */
				string val;
				getline(iss, val); // get the entire remaining part of this line
				if(!tag.empty())
					hmm.setHeaderOpt(tag, BandedHMMP7::trim(val)); // record this tag-value pair
			}
		} /* end of header section */
		else { /* Main body, starts with space */
			iss >> tag;
			if (tag == "COMPO" || BandedHMMP7::isInteger(tag)) { // A compo line can be treated as position 0
				assert((tag == "COMPO" && k == 0) || atoi(tag.c_str()) == k);
				/* process current emission line */
				VectorXf emitFreq = VectorXf(hmm.nuclAbc->getSize());
				for (VectorXf::Index i = 0; i < emitFreq.size(); ++i)
					iss >> emitFreq(i);
				if (tag == "COMPO") { // COMPO line
					emitFreq = (-emitFreq).array().exp();
					hmm.setSpEmissionFreq(emitFreq);
					hmm.hmmBg.setBgFreq(emitFreq);
				} else {
					/* Mk emission line */
					hmm.E_M_log.col(k) = -emitFreq;
					/* Make sure the MAP tag is set */
					string val;
					if(hmm.getHeaderOpt("MAP") != "yes") {
						cerr << "Error: HMM file must has the MAP flag set to 'yes'" << endl;
						abort();
					}
					iss >> tmp;
					hmm.cs2ProfileIdx[atoi(tmp.c_str())] = k;
					hmm.profile2CSIdx[k] = atoi(tmp.c_str());
					hmm.locOptTags["MAP"][k] = tmp;
					/* read other optional tags */
					if(!hmm.getHeaderOpt("CONS").empty()) { /* this tag is present, regarding yes or no */
						iss >> tmp;
						hmm.locOptTags["CONS"][k] = tmp;
					}
					if(!hmm.getHeaderOpt("RF").empty()) { /* this tag is present, regarding yes or no */
						iss >> tmp;
						hmm.locOptTags["RF"][k] = tmp;
					}
					if(!hmm.getHeaderOpt("MM").empty()) { /* this tag is present, regarding yes or no */
						iss >> tmp;
						hmm.locOptTags["MM"][k] = tmp;
					}
					if(!hmm.getHeaderOpt("CS").empty()) { /* this tag is present, regarding yes or no */
						iss >> tmp;
						hmm.locOptTags["CS"][k] = tmp;
					}
				}
				/* process the following Ik emission line */
				for (MatrixXf::Index i = 0; i < hmm.E_I_log.rows(); ++i)
					is >> hmm.E_I_log(i, k);
				/* process the following state K transition line */
					is >> tmp; hmm.T_MM_log(k, k+ 1) = tmp != "*" ? atof(tmp.c_str()) : BandedHMMP7::inf;  // Mk -> Mk+1
					is >> tmp; hmm.T_MI_log(k, k) = tmp != "*" ? atof(tmp.c_str()) : BandedHMMP7::inf;     // Mk -> Ik
					is >> tmp; hmm.T_MD_log(k, k + 1) = tmp != "*" ? atof(tmp.c_str()) : BandedHMMP7::inf; // Mk -> Dk+1
					is >> tmp; hmm.T_IM_log(k, k + 1) = tmp != "*" ? atof(tmp.c_str()): BandedHMMP7::inf;  // Ik -> Mk+1
					is >> tmp; hmm.T_II_log(k, k) = tmp != "*" ? atof(tmp.c_str()) : BandedHMMP7::inf;     // Ik -> Ik
					is >> tmp; hmm.T_DM_log(k, k + 1) = tmp != "*" ? atof(tmp.c_str()) : BandedHMMP7::inf; // Dk -> Mk+1
					is >> tmp; hmm.T_DD_log(k, k + 1) = tmp != "*" ? atof(tmp.c_str()) : BandedHMMP7::inf; // Dk -> Dk+1
			} /* combo line section or match state line section */
			else { // non-COMPO begin state line (M0)
				assert(k == 0);
				string tmp;
				/* process the BEGIN insert emission line */
				for (MatrixXf::Index i = 0; i < hmm.E_I_log.rows(); ++i)
					is >> hmm.E_I_log(i, k);
				/* process the B state K transition line */
				is >> tmp; hmm.T_MM_log(k, k + 1) = tmp != "*" ? atof(tmp.c_str()) : BandedHMMP7::inf; // Mk -> Mk+1
				is >> tmp; hmm.T_MI_log(k, k) = tmp != "*" ? atof(tmp.c_str()) : BandedHMMP7::inf;     // Mk -> Ik
				is >> tmp; hmm.T_MD_log(k, k + 1) = tmp != "*" ? atof(tmp.c_str()) : BandedHMMP7::inf; // Mk -> Dk+1
				is >> tmp; hmm.T_IM_log(k, k + 1) = tmp != "*" ? atof(tmp.c_str()): BandedHMMP7::inf;  // Ik -> Mk+1
				is >> tmp; hmm.T_II_log(k, k) = tmp != "*" ? atof(tmp.c_str()) : BandedHMMP7::inf;     // Ik -> Ik
				is >> tmp; hmm.T_DM_log(k, k + 1) = tmp != "*" ? atof(tmp.c_str()) : BandedHMMP7::inf; // Dk -> Mk+1
				is >> tmp; hmm.T_DD_log(k, k + 1) = tmp != "*" ? atof(tmp.c_str()) : BandedHMMP7::inf; // Dk -> Dk+1
			}
			k++;
		} /* end of main section */
	} /* end of each line */
	// reverse the signs
	hmm.E_I_log = -hmm.E_I_log;

	hmm.T_MM_log = -hmm.T_MM_log;
	hmm.T_MI_log = -hmm.T_MI_log;
	hmm.T_MD_log = -hmm.T_MD_log;
	hmm.T_IM_log = -hmm.T_IM_log;
	hmm.T_II_log = -hmm.T_II_log;
	hmm.T_DM_log = -hmm.T_DM_log;
	hmm.T_DD_log = -hmm.T_DD_log;

	// set the non-log matrices
	hmm.E_M = hmm.E_M_log.array().exp();
	hmm.E_I = hmm.E_I_log.array().exp();
	hmm.T_MM = hmm.T_MM_log.array().exp();
	hmm.T_MI = hmm.T_MI_log.array().exp();
	hmm.T_MD = hmm.T_MD_log.array().exp();
	hmm.T_IM = hmm.T_IM_log.array().exp();
	hmm.T_II = hmm.T_II_log.array().exp();
	hmm.T_DM = hmm.T_DM_log.array().exp();
	hmm.T_DD = hmm.T_DD_log.array().exp();
	hmm.wingRetract();
	return is;
}

ostream& operator<<(ostream& os, const BandedHMMP7& hmm) {
	/* write header tags */
	for(vector<string>::const_iterator it = hmm.headNames.begin(); it != hmm.headNames.end(); ++it)
		os << *it << "  " << hmm.getHeaderOpt(*it) << endl;
	/* write optional HMM tags */
	os << HMM_TAG << endl;
	for(int k = 0; k <= hmm.K; ++k) {
		/* write M or background emission line */
		if(k == 0)
			os << "\tCOMPO\t" << (-hmm.hmmBg.getBgEmitLogPr().transpose())  << endl;
		else {
			os << "\t" << k << "\t" << -hmm.E_M_log.col(k).transpose();
			/* write other optional tags, if present */
			if(!hmm.getHeaderOpt("MAP").empty())
				os << "\t" << hmm.getLocOptTag("MAP", k);
			if(!hmm.getHeaderOpt("CONS").empty())
				os << "\t" << hmm.getLocOptTag("CONS", k);
			if(!hmm.getHeaderOpt("RF").empty())
				os << "\t" << hmm.getLocOptTag("RF", k);
			if(!hmm.getHeaderOpt("MM").empty())
				os << "\t" << hmm.getLocOptTag("MM", k);
			if(!hmm.getHeaderOpt("CS").empty())
				os << "\t" << hmm.getLocOptTag("CS", k);
			os << endl;
		}
		/* write insert emission line */
		for(MatrixXf::Index i = 0; i != hmm.E_I_log.rows(); ++i)
			hmm.E_I_log(i, k) != BandedHMMP7::infV ? os << "\t\t" << hmm.E_I_log(i, k) : os << "\t\t*";
		os << endl;
		/* write state transition line */
		hmm.T_MM_log(k, k + 1) != BandedHMMP7::infV ? os << "\t\t" << -hmm.T_MM_log(k, k + 1) : os << "\t\t*";
		hmm.T_MI_log(k, k) != BandedHMMP7::infV ? os << "\t\t" << -hmm.T_MI_log(k, k) : os << "\t*";
		hmm.T_MD_log(k, k + 1) != BandedHMMP7::infV ? os << "\t\t" << -hmm.T_MD_log(k, k + 1) : os << "\t*";
		hmm.T_IM_log(k, k + 1) != BandedHMMP7::infV ? os << "\t\t" << -hmm.T_IM_log(k, k + 1) : os << "\t*";
		hmm.T_II_log(k, k) != BandedHMMP7::infV ? os << "\t\t" << -hmm.T_II_log(k, k) : os << "\t*";
		hmm.T_DM_log(k, k + 1) != BandedHMMP7::infV ? os << "\t\t" << -hmm.T_DM_log(k, k + 1) : os << "\t*";
		hmm.T_DD_log(k, k + 1) != BandedHMMP7::infV ? os << "\t\t" << -hmm.T_DD_log(k, k + 1) : os << "\t*";
		os << endl;
	}
	os << "//" << endl;
	return os;
}

string BandedHMMP7::trim(const string& str, const string& whitespace) {
	const string::size_type strBegin = str.find_first_not_of(whitespace);
	if(strBegin == string::npos) // no content
		return "";
	string::size_type strRange = str.find_last_not_of(whitespace) - strBegin + 1;
	return str.substr(strBegin, strRange);
}

} /* namespace EGriceLab */

/*void EGriceLab::BandedHMMP7::setSequenceMode(enum align_mode mode) {
 switch(mode) {
 case GLOBAL:
 }
 }*/

EGriceLab::BandedHMMP7::BandedHMMP7() :
		version(progName + "-" + progVersion), name("NULL"), K(0), nuclAbc(
				SeqCommons::nuclAbc), hmmBg(0), transInit(false), emisInit(false),
				spInit(false), limitInit(false), indexInit(false), wingRetracted(false) {
	/* Assert IEE559 at construction time */
	assert(numeric_limits<float>::is_iec559);
	/* set mandatory tags */
	setHeaderOpt("HMMER3/f", progName + "-" + progVersion);
	setHeaderOpt("NAME", "NULL");
	setHeaderOpt("LENG", "0");
	setHeaderOpt("ALPH", "DNA");
}

EGriceLab::BandedHMMP7::BandedHMMP7(const string& name, int K, const IUPACNucl* abc) :
		version(progName + "-" + progVersion), name(name), K(K), nuclAbc(abc),
		hmmBg(K), transInit(false), emisInit(false),
				spInit(false), limitInit(false), indexInit(false), wingRetracted(false) {
	/* Assert IEE559 at construction time */
	assert(numeric_limits<float>::is_iec559);
	init_transition_params();
	init_emission_params();
	init_special_params();
	init_index();
	init_limits();
	enableProfileLocalMode(); // always in profile local alignment mode
	setSpEmissionFreq(hmmBg.getBgEmitPr()); // Use bg emission freq by default
	/* set mandatory tags */
	setHeaderOpt("HMMER3/f", progName + "-" + progVersion);
	setHeaderOpt("NAME", "NULL");
	string size_s;
	ostringstream os(size_s);
	os << K;
	setHeaderOpt("LENG", size_s);
	setHeaderOpt("ALPH", "DNA");
}

void EGriceLab::BandedHMMP7::setProfileSize(int size) {
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

void EGriceLab::BandedHMMP7::setSequenceMode(enum align_mode mode) {
	switch (mode) {
	case GLOBAL:
		T_SP(N, N) = T_SP(C, C) = 0;
		break;
	case LOCAL:
		T_SP(N, N) = T_SP(C, C) = 1 - hmmBg.getBgTransPr();
		break;
	case NGCL:
		T_SP(N, N) = 0;
		T_SP(C, C) = 1 - hmmBg.getBgTransPr();
		break;
	case CGNL:
		T_SP(N, N) = 1 - hmmBg.getBgTransPr();
		T_SP(C, C) = 0;
		break;
	default:
		cerr
				<< "Unknown seqAlign type, must be one of 'GLOBAL', 'LOCAL', 'NGCL' or 'CGNL'"
				<< endl;
		abort();
	}
	T_SP(N, B) = 1 - T_SP(N, N);
	T_SP(E, C) = 1; // always exit from E->C
	T_SP_log = T_SP.array().log(); // Eigen3 handle array to matrix assignment automatically
}

void EGriceLab::BandedHMMP7::setSpEmissionFreq(const VectorXf& freq) {
	E_SP.col(N) = E_SP.col(C) = freq / freq.sum(); // re-do normalization, even if already done
	E_SP.col(B) = E_SP.col(E) = VectorXf::Zero(nuclAbc->getSize()); // no emission for state B and E
	E_SP_log = E_SP.array().log();
}

void EGriceLab::BandedHMMP7::init_transition_params() {
	if (transInit) // already initialized
		return;
	/*
	 * state 0 serves as the B state
	 * state K + 1 serve as E state
	 */
	T_MM = T_MM_log = MatrixXf(K + 2, K + 2);
	T_II = T_II_log = MatrixXf(K + 2, K + 2);
	T_DD = T_DD_log = MatrixXf(K + 2, K + 2);

	T_IM = T_IM_log = MatrixXf(K + 2, K + 2);
	T_MI = T_MI_log = MatrixXf(K + 2, K + 2);

	T_DM = T_DM_log = MatrixXf(K + 2, K + 2);
	T_MD = T_MD_log = MatrixXf(K + 2, K + 2);

	T_MM_retract = T_MM_retract_log = MatrixXf(K + 2, K + 2);
	transInit = true;
}

void EGriceLab::BandedHMMP7::init_emission_params() {
	if (emisInit) // already initialized
		return;
	/* state 0 serves as B state and will not be used */
	E_M = E_M_log = MatrixXf(nuclAbc->getSize(), K + 1);
	E_I = E_I_log = MatrixXf(nuclAbc->getSize(), K + 1);
	emisInit = true;
}

void EGriceLab::BandedHMMP7::init_special_params() {
	if (spInit) // already initialized
		return;
	E_SP = E_SP_log = MatrixXf(nuclAbc->getSize(), kNSP);
	T_SP = T_SP_log = MatrixXf(kNSP, kNSP);
	spInit = true;
}

void EGriceLab::BandedHMMP7::init_limits() {
	if (limitInit) // already initialized
		return;
/*	insBeforeLimit = insAfterLimit = VectorXi::Constant(K + 1, static_cast<int> (kMaxProfile));
	delBeforeLimit = delAfterLimit = VectorXi::Constant(K + 1, static_cast<int> (kMaxProfile));*/

	gapBeforeLimit = gapAfterLimit = VectorXi(K + 1);
	//delBeforeLimit = delAfterLimit = VectorXi(K + 1);
	for(VectorXi::Index j = 1; j <= K; ++j) {
		gapBeforeLimit(j) = j * kMaxGapFrac;
		gapAfterLimit(j) = (K - j) * kMaxGapFrac;
	}
	limitInit = true;
}

void EGriceLab::BandedHMMP7::init_index() {
	if(indexInit) // already initialized
		return;
	for(int i = 0; i <= kMaxProfile; ++i)
		cs2ProfileIdx[i] = 0;
	for(int i = 1; i <= kMaxCS; ++i)
		profile2CSIdx[i] = 0;
}

void EGriceLab::BandedHMMP7::enableProfileLocalMode() {
	/* set entering probabilities */
	T_MM(0, 0) = T_MM(0, K + 1) = 0; // B->B, B->E not allowed
	T_MM.row(0).segment(1, K).setConstant((1.0 - hmmBg.getBgTransPr()) / K); /* B->M1..MK equal probability */


	T_MI.block(0, 1, 1, K + 1).setConstant(0); // B -> I1..IK+1 not allowed; B->I0 leave as is

	T_MD(0, 0) = T_MD(0, K + 1) = 0; // B->B, B->E not allowed
	T_MD(0, 1) = hmmBg.getBgTransPr() / K; // B->D1 as background probability
	T_MD.row(0).segment(2, K - 1).setConstant(0); // B -> D2...DK not allowed

	T_IM.row(0).setConstant(0); // no I0 state
	T_II.row(0).setConstant(0); // no I0 state
	T_DM.row(0).setConstant(0); // no D0 state
	T_DD.row(0).setConstant(0); // no D0 state

	/* set exiting probabilities */
	T_MM(K + 1, K + 1) = 0; // E->E not allowed
	T_MM.col(K + 1).segment(1, K).setConstant((1.0 - hmmBg.getBgTransPr()) / K); /* M1..MK ->E equal probability */

	T_MI.col(K + 1).setConstant(0); // no IK+1 state
	T_MD.col(K + 1).setConstant(0); // no DK+1 state

	T_IM.col(K + 1).setConstant(0); // I0..IK+1 -> E not allowed
	T_II.col(K + 1).setConstant(0); // no IK+1 state

	T_DM(0, K + 1) = T_DM(K + 1, K + 1) = 0; // D0 -> E and DK+1 -> E not allowed
	T_DM(K, K + 1) = hmmBg.getBgTransPr() / K; // DK -> E as background probability
	T_DM.block(1, K + 1, K - 1, 1).setConstant(0); // D1..DK-1 -> E not allowed

	T_DD.col(K + 1).setConstant(0); // no DK+1 state

	/* set log version */
	T_MM_log.row(0) = T_MM.row(0).array().log();
	T_MI_log.row(0) = T_MI.row(0).array().log();
	T_MD_log.row(0) = T_MD.row(0).array().log();
	T_IM_log.row(0) = T_IM.row(0).array().log();
	T_II_log.row(0) = T_II.row(0).array().log();
	T_DM_log.row(0) = T_DM.row(0).array().log();
	T_DD_log.row(0) = T_DD.row(0).array().log();
}
EGriceLab::BandedHMMP7::ViterbiScores EGriceLab::BandedHMMP7::initViterbiScores(
		const PrimarySeq& seq) const {
	assert(nuclAbc == seq.getDegenAlphabet());
	ViterbiScores vs(seq);

	return resetViterbiScores(vs, seq);
/*	resetForwardBackwardScores(bHmm);
	resetTraceScores(bHmm);*/
}

EGriceLab::BandedHMMP7::ViterbiScores& EGriceLab::BandedHMMP7::resetViterbiScores(ViterbiScores& vs,
		const PrimarySeq& seq) const {
	vs.seq = &seq; // swap the pointer value but leave the pointing object untouched
	const int L = vs.seq->length();
	/* Use -Inf as initial values for Viterbi matrices */
/*	vs.DP_M = vs.DP_I = vs.DP_D = vs.S = MatrixXf::Constant(L + 1, K + 1,
			-numeric_limits<float>::infinity());*/
	/* set dimension */
	vs.DP_M.resize(L + 1, K + 1);
	vs.DP_I.resize(L + 1, K + 1);
	vs.DP_D.resize(L + 1, K + 1);

	/* set elements */
	vs.DP_M.setConstant(infV);
	vs.DP_I.setConstant(infV);
	vs.DP_D.setConstant(infV);
	return vs;
}

/*void EGriceLab::BandedHMMP7::resetForwardBackwardScores(
		ViterbiScores& vs) const {
	int L = vs.ds.length();
	 Use 0 as initial values for forward-backward matricies
	vs.F_M = vs.F_I = vs.F_D = MatrixXf::Zero(L + 1, K + 1);
	vs.B_M = vs.B_I = vs.B_D = MatrixXf::Zero(L + 1, K + 1);
	 Use 1 as default scaling factor value
	vs.SV_M = vs.SV_I = vs.SV_D = VectorXf::Constant(L + 1, 1);
}*/

/*void EGriceLab::BandedHMMP7::resetTraceScores(ViterbiScores& vs) const {
	int L = vs.ds.length();
	 Use -1 as initial values for trace scores
	vs.TRACE = MatrixXi::Constant(L + 1, K + 1, -1);
}*/


EGriceLab::BandedHMMP7::ViterbiAlignPath EGriceLab::BandedHMMP7::initViterbiAlignPath(
		int L) const {
	ViterbiAlignPath vpath;
	vpath.K = K; // K is fixed between different seq
	vpath.N = 0;
	return resetViterbiAlignPath(vpath, L);
}

EGriceLab::BandedHMMP7::ViterbiAlignPath& EGriceLab::BandedHMMP7::resetViterbiAlignPath(ViterbiAlignPath& vpath,
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

void EGriceLab::BandedHMMP7::calcViterbiScores(
		ViterbiScores& vs) const {
	assert(vs.seq != NULL);
	assert(wingRetracted); // make sure wing is retracted
	const int L = vs.seq->length();
/*	vs.DP_M = vs.DP_I = vs.DP_D = MatrixXf::Constant(L + 1, K + 1, -numeric_limits<float>::infinity());*/
	/* vs.DP_M(0, 0) = 0; */
	/* Initialize the M(,0), the B state */
	for (int i = 1; i <= L; ++i)
		if(i == 1) /* No N->N loop */
			vs.DP_M(i, 0) = 0;
		else
			vs.DP_M(i, 0) = T_SP_log(N, N) * (i - 1); /* (i-1) N state circles */
	vs.DP_M.col(0).array() += T_SP_log(N, B); /* N->B transition */
	/* Full Dynamic-Programming at row-first order */
	for (int j = 1; j <= K; ++j) {
		for (int i = 1; i <= L; ++i) {
			vs.DP_M(i, j) = E_M_log(vs.seq->encodeAt(i-1), j) + EGriceLab::BandedHMMP7::max(
							static_cast<float> (vs.DP_M(i, 0) + T_MM_retract_log(0, j)), // from M(i,0), the B state
							static_cast<float> (vs.DP_M(i - 1, j - 1) + T_MM_retract_log(j - 1, j)), // from Mi-1,j-1
							static_cast<float> (vs.DP_I(i - 1, j - 1) + T_IM_log(j - 1, j)), // from Ii-1,j-1
							static_cast<float> (vs.DP_D(i - 1, j - 1) + T_DM_log(j - 1, j))); // from Di-1,j-1
			vs.DP_I(i, j) = E_I_log(vs.seq->encodeAt(i - 1), j) + std::max(
							static_cast<float> (vs.DP_M(i - 1, j) + T_MI_log(j, j)), // from Mi-1,j
							static_cast<float> (vs.DP_I(i - 1, j) + T_II_log(j, j))); // from Ii-1,j
			vs.DP_D(i, j) = std::max(
							static_cast<float> (vs.DP_M(i, j - 1) + T_MD_log(j - 1, j)), // from Mi,j-1
							static_cast<float> (vs.DP_D(i, j - 1) + T_DD_log(j - 1, j))); // from Di,j-1
		}
	}
	vs.S = vs.DP_M; // final Viterbi scores, copied from the calculated DP_M;
	for (int j = 1; j <= K; j++)
		vs.S.col(j).array() += T_MM_retract_log(j, K + 1); // add M->E transition
	vs.S.array() += T_SP_log(E, C); // add E->C transition
	for (int i = 1; i < L; ++i)
		vs.S.row(i).array() += T_SP_log(C, C) * (L - i); // add L-i C->C circles
	vs.S.col(0).setConstant(infV); // set the B state scores to infV

}

void EGriceLab::BandedHMMP7::calcViterbiScores(
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
	for (int i = 1; i <= L; i++)
		if(i == 1) // no N->N loop if entering from 1
			vs.DP_M(i, 0) = 0;
		else
			vs.DP_M(i, 0) = T_SP_log(N, N) * (i - 1); /* (i-1) N states */

	vs.DP_M.col(0).array() += T_SP_log(N, B); /* N->B transition */
	//cerr << "M(,0)" << endl << vs.DP_M.col(0) << endl;
/*	cerr << "M(,0) initialized" << endl;*/
	/* process each known path upstream and themselves */
	for(vector<int>::size_type n = 0; n < vpath.N; ++n) {
		/* Determine banded boundaries */
		int upQLen = n == 0 ? vpath.from[n] - 1 : vpath.from[n] - vpath.to[n - 1];
		int up_start = n == 0 ? vpath.start[n] - upQLen * (1 + kMaxGapFrac) : vpath.end[n - 1];
		if (up_start < 1)
			up_start = 1;
		int up_from = n == 0 ? vpath.from[n] - upQLen * (1 + kMaxGapFrac) : vpath.to[n - 1];
		if (up_from < 1)
			up_from = 1;
		//cerr << "upQLen:" << upQLen << endl;
//		cerr << "up_start:" << up_start << " up_end:" << vpath.start[n] << endl;
//		cerr << "up_from:" << up_from << " up_to:" << vpath.from[n] << endl;

		/* Dynamic programming of upstream of this known path at row-first order */
		for (int j = up_start; j <= vpath.start[n]; ++j) {
			for (int i = up_from; i <= vpath.from[n]; ++i) {
				vs.DP_M(i, j) = E_M_log(vs.seq->encodeAt(i - 1), j)
							+ EGriceLab::BandedHMMP7::max(
									static_cast<float>(vs.DP_M(i, 0) + T_MM_retract_log(0, j)), // from Mi,0, the B state
									static_cast<float>(vs.DP_M(i - 1, j - 1)
											+ T_MM_retract_log(j - 1, j)), // from Mi-1,j-1
											static_cast<float>(vs.DP_I(i - 1, j - 1)
													+ T_IM_log(j - 1, j)), // from Ii-1,j-1
													static_cast<float>(vs.DP_D(i - 1, j - 1)
															+ T_DM_log(j - 1, j))); // from Di-1,j-1
				vs.DP_I(i, j) = E_I_log(vs.seq->encodeAt(i - 1), j)
							+ std::max(
									static_cast<float>(vs.DP_M(i - 1, j)
											+ T_MI_log(j, j)), // from Mi-1,j
											static_cast<float>(vs.DP_I(i - 1, j)
													+ T_II_log(j, j))); // from Ii-1,j
				vs.DP_D(i, j) =
						std::max(
								static_cast<float>(vs.DP_M(i, j - 1)
										+ T_MD_log(j - 1, j)), // from Mi,j-1
										static_cast<float>(vs.DP_D(i, j - 1)
												+ T_DD_log(j - 1, j))); // from Di,j-1
			}
		}
		/* Fill the score of the known alignment path */
		for (int j = vpath.start[n]; j <= vpath.end[n]; ++j) {
			for(int i = vpath.from[n]; i <= vpath.to[n]; ++i) {
				int dist = diagnalDist(i, j, vpath.from[n], vpath.start[n]);
				if(!(dist <= vpath.numIns[n] && dist >= -vpath.numDel[n]))
					continue;
				vs.DP_M(i, j) = E_M_log(vs.seq->encodeAt(i - 1), j)
							+ EGriceLab::BandedHMMP7::max(
									static_cast<float>(vs.DP_M(i, 0) + T_MM_retract_log(0, j)), // from M(,0), the B state
									static_cast<float>(vs.DP_M(i - 1, j - 1)
											+ T_MM_retract_log(j - 1, j)), // from Mi-1,j-1
											static_cast<float>(vs.DP_I(i - 1, j - 1)
													+ T_IM_log(j - 1, j)), // from Ii-1,j-1
													static_cast<float>(vs.DP_D(i - 1, j - 1)
															+ T_DM_log(j - 1, j))); // from Di-1,j-1
				vs.DP_I(i, j) = E_I_log(vs.seq->encodeAt(i - 1), j)
							+ std::max(
									static_cast<float>(vs.DP_M(i - 1, j)
											+ T_MI_log(j, j)), // from Mi-1,j
											static_cast<float>(vs.DP_I(i - 1, j)
													+ T_II_log(j, j))); // from Ii-1,j
				vs.DP_D(i, j) =
						std::max(
								static_cast<float>(vs.DP_M(i, j - 1)
										+ T_MD_log(j - 1, j)), // from Mi,j-1
										static_cast<float>(vs.DP_D(i, j - 1)
												+ T_DD_log(j - 1, j))); // from Di,j-1
			}
		}
		// assert(i == vpath.to + 1 && j == vpath.end + 1);
	} /* end of each known path segment
	/* Dynamic programming of the remaining downstream of the known paths, if any */
	int last_end = vpath.end[vpath.N - 1];
	int last_to = vpath.to[vpath.N - 1];
	int downQLen = L - last_to;
	int down_end = last_end + downQLen * (1 + kMaxGapFrac);
	int down_to = last_to + downQLen * (1 + kMaxGapFrac);

	for (int j = last_end; j <= down_end && j <= K; ++j) {
		for (int i = last_to; i <= down_to && i <= L; ++i) {
			vs.DP_M(i, j) = E_M_log(vs.seq->encodeAt(i - 1), j) + EGriceLab::BandedHMMP7::max(
							// from Mi,0, the B state is not possible
							static_cast<float>(vs.DP_M(i - 1, j - 1)
									+ T_MM_retract_log(j - 1, j)), // from Mi-1,j-1
							static_cast<float>(vs.DP_I(i - 1, j - 1)
									+ T_IM_log(j - 1, j)), // from Ii-1,j-1
							static_cast<float>(vs.DP_D(i - 1, j - 1)
									+ T_DM_log(j - 1, j))); // from Di-1,j-1
			vs.DP_I(i, j) = E_I_log(vs.seq->encodeAt(i - 1), j) + std::max(
							static_cast<float>(vs.DP_M(i - 1, j)
									+ T_MI_log(j, j)), // from Mi-1,j
							static_cast<float>(vs.DP_I(i - 1, j)
									+ T_II_log(j, j))); // from Ii-1,j
			vs.DP_D(i, j) =
					std::max(
							static_cast<float>(vs.DP_M(i, j - 1)
									+ T_MD_log(j - 1, j)), // from Mi,j-1
							static_cast<float>(vs.DP_D(i, j - 1)
									+ T_DD_log(j - 1, j))); // from Di,j-1
		}
	}
	//cerr << "downstream done " << vpath.end + 1 << "-" << down_end << " " << vpath.to + 1 << "-" << down_to << endl;
	vs.S = vs.DP_M; // final Viterbi scores, copied from the calculated DP_M;*/
/*	cerr << "vs assigned" << endl;*/

	for (int j = 1; j <= K; j++)
		vs.S.col(j).array() += T_MM_retract_log(j, K + 1); /* add M->E transition */
/*	if(vs.ds.getName() == "r131")
		cerr << "S(L,):" << endl << vs.S.row(L) << endl;
	cerr << "M->E added" << endl;*/

	vs.S.array() += T_SP_log(E, C); // add E->C transition
/*	if(vs.ds.getName() == "r131")
		cerr << "S(L,):" << endl << vs.S.row(L) << endl;
	cerr << "E->C added" << endl;*/

	for (int i = 1; i < L; ++i)
		vs.S.row(i).array() += T_SP_log(C, C) * (L - i); // L-i C->C circles
/*	if(vs.ds.getName() == "r131")
		cerr << "S(L,):" << endl << vs.S.row(L) << endl;
	cerr << "C-C circles added" << endl;*/

	vs.S.col(0).setConstant(infV); // set the B state scores to infV
	/*cerr << "DP_M:" << endl << vs.DP_M << endl;
	cerr << "DP_I:" << endl << vs.DP_I << endl;
	cerr << "DP_D:" << endl << vs.DP_D << endl;
	cerr << "S:" << endl << vs.S << endl;*/
	return;
}

void EGriceLab::BandedHMMP7::addKnownAlignPath(ViterbiAlignPath& vpath,
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

		bool nonGap = nuclAbc->isSymbol(*it) && !nuclAbc->isGap(*it);

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

float EGriceLab::BandedHMMP7::buildViterbiTrace(const ViterbiScores& vs, ViterbiAlignPath& vpath) const {
	assert(vs.seq != NULL);
	assert(vs.S.rows() == vpath.L + 1 && vs.S.cols() == vpath.K + 1);
	MatrixXf::Index maxRow, maxCol;
	float maxScore = vs.S.maxCoeff(&maxRow, &maxCol);
	if(maxScore == infV)
		return maxScore; // max score not found, do nothing
	/* do trace back in the vScore matrix */
	int i = maxRow, j = maxCol;
	vpath.alnStart = vpath.alnEnd = maxCol;
	vpath.alnFrom = vpath.alnTo = maxRow;
	//cerr << "maxRow:" << maxRow << " maxCol:" << maxCol << " maxScore:" << maxScore << " maxEnd:" << getCSLoc(maxCol) << endl;

	char s = 'M'; // current state (p7 always found maximum at M state)
	vpath.alnPath.push_back('E'); // ends with E
	while(i >= 0 && j >= 0) {
		//vpath.TRACE(i, j) = s;
		vpath.alnPath.push_back(s);
		if(j == 0) {
			//cerr << "i:" << i << " j:" << j << " s:" << decode(B) << ":" << static_cast<float> (vs.DP_M(i, j)) << endl;
			break;
		}
		// update the status
		if(s == 'M') {
			//cerr << "i:" << i << " j:" << j << " s:" << decode(s) << ":" << static_cast<float> (vs.DP_M(i, j));
			vpath.alnStart--;
			vpath.alnFrom--;
			s = BandedHMMP7::whichMax(
					static_cast<float> (vs.DP_M(i, 0) + T_MM_retract_log(0, j)), /* from M(i, 0), the B-state */
					static_cast<float> (vs.DP_M(i - 1, j - 1) + T_MM_retract_log(j - 1, j)), /* from M(i-1,j-1) */
					static_cast<float> (vs.DP_I(i - 1, j - 1) + T_IM_log(j - 1, j)), /* from I(i-1,j-1) */
					static_cast<float> (vs.DP_D(i - 1, j - 1) + T_DM_log(j - 1, j))); /* from D(i-1,j-1) */
			//cerr << "->" << decode(s) << endl;
			if(s == 'B')
				j = 0; // jump to B state
			else { // M, I or D state
				i--;
				j--;
			}

		}
		else if(s == 'I') {
			//cerr << "i:" << i << " j:" << j << " s:" << decode(s) << ":" << static_cast<float> (vs.DP_I(i, j));
			vpath.alnFrom--;
			s = BandedHMMP7::whichMax(
					static_cast<float> (vs.DP_M(i - 1, j) + T_MI_log(j, j)), /* from M(i-1,j) */
					static_cast<float> (vs.DP_I(i - 1, j) + T_MI_log(j, j)), /* from I(i-1,j) */
					"MI");
			//cerr << "->" << decode(s) << endl;
			i--;
		}
		else if(s == 'D') {
			vpath.alnStart--;
			//cerr << "i:" << i << " j:" << j << " s:" << decode(s) << ":" << static_cast<float> (vs.DP_D(i, j));
			s = BandedHMMP7::whichMax(
					static_cast<float> (vs.DP_M(i, j - 1) + T_MD_log(j - 1, j)), /* from M(i,j-1) */
					static_cast<float> (vs.DP_D(i, j - 1) + T_DD_log(j - 1, j)), /* from D(i,j-1) */
					"MD");
			//cerr << "->" << decode(s) << endl;
			j--;
		}
		else { // B state found
			//cerr << "i:" << i << " j:" << j << " s:" << decode(s) << ":" << static_cast<float> (vs.DP_M(i, j)) << endl;
			break;
		}
		//cerr << "vpath.alnStart:" << vpath.alnStart + 1 << " vpath.alnEnd:" << vpath.alnEnd << endl;
	} /* end of while */
	vpath.alnStart++; // 1-based
	vpath.alnFrom++; // 1-based
	if(vpath.alnPath[vpath.alnPath.length() - 1] != 'B')
		vpath.alnPath.push_back('B');
	reverse(vpath.alnPath.begin(), vpath.alnPath.end()); // reverse the alnPath string
	return maxScore;
}

std::string EGriceLab::BandedHMMP7::buildGlobalAlignSeq(const ViterbiScores& vs,
		ViterbiAlignPath& vpath) const {
	assert(vs.seq != NULL);
	string aSeq;
	int profileNLen = vpath.alnStart - 1;
	int profileCLen = vpath.K - vpath.alnEnd;
	int seqNLen = vpath.alnFrom - 1;
	int seqCLen = vpath.L - vpath.alnTo;

	/* place N state residues to the beginning and ending of the unmatched up-stream of the profile */
	int NHalf = std::min(profileNLen / 2, seqNLen / 2);
	int CHalf = std::min(profileCLen / 2, seqCLen / 2);
	for(int i = 0; i < NHalf; ++i) /* put half the N residues at the beginning of the N' */
		aSeq.push_back(vs.seq->charAt(i - 1));
	aSeq.append(profileNLen - 2 * NHalf, '-'); /* add the middle gaps, if any */
	for(int i = profileNLen - NHalf; i < profileNLen; ++i) /* put half the N residues at the beginning */
		aSeq.push_back(vs.seq->charAt(i - 1)); /* put half the N residues at the end of the N' */

	/* place aligned residules */
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
	for(int i = vpath.alnTo; i < vpath.alnTo + CHalf; ++i) /* put half the N residues at the begining of the N' */
		aSeq.push_back(vs.seq->charAt(i - 1));
	aSeq.append(profileCLen - 2 * CHalf, '-'); /* add the middle gaps, if any */
	for(int i = vpath.K - NHalf; i < vpath.K; ++i) /* put half the N residues at the begining */
		aSeq.push_back(vs.seq->charAt(i - 1)); /* put half the N residues at the end of the N' */
	return aSeq;
}

void EGriceLab::BandedHMMP7::wingRetract() {
	if(wingRetracted) // already wing-retracted
		return;
	T_MM_retract = T_MM; // copy the transition
	/* retract entering probabilities B->Mj */
	for(MatrixXf::Index j = 2; j <= K; ++j) {
		float logP = 0; // additional retract probability in log-scale
		logP += T_MD_log(0, 1); // B->D1
		for(MatrixXf::Index i = 1; i < j - 1; ++i)
			logP += T_DD_log(i, i + 1); // Di->Di+1
		logP += T_DM_log(j - 1, j); // Dj-1->Mj
		T_MM_retract(0, j) += exp(logP); // retract B->D1->D2...Dj-1->Mj to B->Mj
		if(T_MM_retract(0, j) > 1)
			T_MM_retract(0, j) = 1;
	}
	/* retract exiting probabilities Mi -> E */
	for(MatrixXf::Index i = 1; i <= K - 1; ++i) {
		float logP = 0; // additional retract probability in log-scale
		logP += T_MD_log(i, i + 1); // Mi -> Di+1
		for(MatrixXf::Index j = i + 1; j < K; ++j)
			logP += T_DD_log(j, j + 1); // Dj->Dj+1
		logP += T_DM_log(K, K + 1); // DK -> E
		T_MM_retract(i, K + 1) += exp(logP); // retract Mi->Di+1->Di+2...->DK->E to Mi->E
		if(T_MM_retract(i, K + 1) > 1)
			T_MM_retract(i, K + 1) = 1;
	}
	T_MM_retract_log = T_MM_retract.array().log();
	wingRetracted = true;
}

bool EGriceLab::BandedHMMP7::isValidAlignPath(const ViterbiAlignPath& vpath) const {
	bool flag = true;
	for(int n = 0; n < vpath.N; ++n)
		if(!(vpath.start[n] >= 1 && vpath.start[n] <= vpath.end[n] && vpath.end[n] <= vpath.K
			&& vpath.from[n] >= 1 && vpath.from[n] <= vpath.to[n] && vpath.to[n] <= vpath.L) && // this path is valid
				(n == 0 || (vpath.start[n] > vpath.end[n - 1] && vpath.from[n] > vpath.to[n - 1]))) // path is properly ordered
			return false;
	return true;
}
