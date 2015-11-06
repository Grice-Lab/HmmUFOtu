/*
 * CSFMIndex.h
 *
 *  Created on: Nov 5, 2015
 *      Author: zhengqi
 */

#ifndef CSFMINDEX_H_
#define CSFMINDEX_H_
#include "divsufsort.h"
#include "MSA.h"

namespace EGriceLab {

/**
 * A Consensus-Sequence FM-index for ultra-fast indexing the consensus positions of a multiple-sequence alignment
 */
class CSFMIndex {
public:
	CSFMIndex();
	virtual ~CSFMIndex();

private:
	/* private default constructor */
	CSFMIndex();
	/* disable the copy and assignment constructor */
	CSFMIndex(const CSFMIndex& other);
	CSFMIndex& operator=(const CSFMIndex& other);


};

} /* namespace EGriceLab */

#endif /* CSFMINDEX_H_ */
