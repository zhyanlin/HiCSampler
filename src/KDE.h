/*
 * KDE.h
 *
 *  Created on: 4 Feb 2018
 *      Author: yanlin
 */

#ifndef KDE_H_
#define KDE_H_
#include "contactMap.h"
#include <utility>
#include <queue>

class KDE
{
public:
	KDE();
	virtual ~KDE();
	int initialize();
	int runKDE();

private:
	int nbStdDevTh;
	int nSet;
	double gaussianSum;
	float maxRadius2;
	float maxRadius;
	float **gaussian;
	int *minb;
	//static void* applyGaussianFilter(void* object);

	static void *applyGaussianFilter(void *object)
	{
		reinterpret_cast<KDE *>(object)->applyGaussianFilter();
		return 0;
	}
	int applyGaussianFilter();
	UpperDiag<unsigned long> *cumO;
	UpperDiag<int> *nextNonEmpty;
	float mini_h;
	float maxi_h;
	int rowStep, colStep;
	std::queue<std::pair<int, int>> job;
	double multFactor;
	bool lastIter;
};

#endif /* KDE_H_ */
