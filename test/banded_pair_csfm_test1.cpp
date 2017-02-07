/*
 * banded_align_test.cpp
 *
 *  Created on: Jun 18, 2015
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include "HmmUFOtu_common.h"
#include "HmmUFOtu_hmm.h"

using namespace std;
using namespace EGriceLab;

int main(int argc, char *argv[]) {
	if(argc != 5) {
		cerr << "Usage:  " << argv[0] << " HMM-FILE CSFM-DB TRAINING-FILE OUTFILE" << endl;
		return 0;
	}

	ifstream hmm_in(argv[1]);
	ifstream idx_in(argv[2]);
	ifstream in(argv[3]);
	ofstream out(argv[4]);

	// read in HMM profile
	BandedHMMP7 hmm;
	//cerr << "hmm initiated" << endl;
	hmm_in >> hmm;
	//cerr << "hmm profile read" << endl;

	hmm.setSequenceMode(BandedHMMP7::GLOBAL);

	//cerr << "hmm sequence mode set to global" << endl;

	hmm.wingRetract();
	//cerr << "wing retracted" << endl;

	// load CSFM index
	CSFMIndex idx;
	idx.load(idx_in);

	// process each training read
	string line;
	char id[64];
	int seq_start, seq_end, fwd_seed_from, fwd_seed_to, rev_seed_from, rev_seed_to;
	char fwd_seed_CS[4096], rev_seed_CS[4096], seq_CS[4096];
	getline(in, line); // ignore header line
	out << line << "\tknown_CS_start\tknown_CS_end\tknown_hmm_start\tknown_hmm_end\tbHmm_start\tbHmm_end\tCS_start\tCS_end" << endl;
	BandedHMMP7::ViterbiScores vscore = hmm.initViterbiScores(); // construct an empty score object
	BandedHMMP7::ViterbiAlignPath vpath = hmm.initViterbiAlignPath(0); // a reusable object

	while(getline(in, line)) {
		sscanf(line.c_str(), "%*s\t%*s\t%*d\t%*d\t%*s\t%*s\t%*d\t%*d\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s",
				&seq_start, &seq_end, id,
				&fwd_seed_from, &fwd_seed_to, fwd_seed_CS,
				&rev_seed_from, &rev_seed_to, rev_seed_CS,
				seq_CS);
		// Determine seed_start and seed_end
		PrimarySeq fwd_seed(hmm.getNuclAbc(), id, fwd_seed_CS);
		PrimarySeq rev_seed(hmm.getNuclAbc(), id, rev_seed_CS);
		fwd_seed.removeGaps();
		rev_seed.removeGaps();

		const CSLoc& fwd_loc = idx.locateOne(fwd_seed.getSeq());
		const CSLoc& rev_loc = idx.locateOne(rev_seed.getSeq());
		//cerr << "CS: start:" << loc.start << " end:" << loc.end << " CS:" << loc.CS << endl;

		//cerr << "known_start:" << known_start << " known_end:" << known_end << endl;
		//cerr << "seed_start:" << seed_start << " seed_end:" << seed_end << endl;

		/* calculate known seq start and end */
		while(hmm.getProfileLoc(seq_start) == 0) // D loc
			seq_start++;
		while(hmm.getProfileLoc(seq_end) == 0) // D loc
			seq_end--;

		/* calculate known hmm start and end */
		int hmm_start = hmm.getProfileLoc(seq_start);
		int hmm_end = hmm.getProfileLoc(seq_end);

//		cerr << "id:" << id << " length:" << seq.length() << endl;
//		cerr << "constructed the seq" << endl;
		PrimarySeq seq(hmm.getNuclAbc(), id, seq_CS);
		seq.removeGaps();
		hmm.resetViterbiScores(vscore, seq);
		/*cerr << "constructed the ViterbiScore" << endl;*/

		hmm.resetViterbiAlignPath(vpath, seq.length());
		//cerr << "vpath.start:" << vpath.start << " vpath.end:" << vpath.end << endl;
	/*cerr << "constructed the ViterbiAlignPath" << endl;*/

		if(fwd_loc.start > 0 && fwd_loc.end > 0)
			hmm.addKnownAlignPath(vpath, fwd_loc, fwd_seed_from, fwd_seed_to);
		if(rev_loc.start > 0 && rev_loc.end > 0)
			hmm.addKnownAlignPath(vpath, rev_loc, rev_seed_from, rev_seed_to);
		//cerr << "set the known path" << endl;

		hmm.calcViterbiScores(vscore, vpath);
		//hmm.calcViterbiScores(vscore);
		//cerr << "Bhmm aligned" << endl;

		float maxScore = hmm.buildViterbiTrace(vscore, vpath);
		if(maxScore == infV) {
			cerr << id << " no maxScore" << endl;
			continue;
		}
		//cerr << "back traced" << endl;

/*		string align = hmm.buildGlobalAlignSeq(vscore, vpath);*/
		// output
		out << line << "\t" << seq_start << "\t" << seq_end << "\t" << hmm_start << "\t" << hmm_end << "\t" <<
				vpath.alnStart << "\t" << vpath.alnEnd << "\t" <<
				hmm.getCSLoc(vpath.alnStart) << "\t" << hmm.getCSLoc(vpath.alnEnd) << endl;
				// "\t" << vpath.alnPath << endl;
	}

	in.close();
	out.close();
	return 0;
}

