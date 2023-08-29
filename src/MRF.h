/*
 * MRF.h
 *
 *  Created on: 4 Feb 2018
 *      Author: yanlin
 */

#ifndef MRF_H_
#define MRF_H_
#include <queue>
#include "contactMap.h"

class MRF
{
public:
	int initialize();
	int runMRF();

private:
	UpperDiag<char> *hasChanged;
	unsigned long nChanges;
	int rep;
	int rowStep, colStep;
	std::queue<std::pair<int, int>> job;
	static void *blockMRF(void *object)
	{
		reinterpret_cast<MRF *>(object)->blockMRF();
		return 0;
	}
	int blockMRF();
};

#endif /* MRF_H_ */
