//This file is based on the paper "Fast wavefront propagation (FWP) for computing exact geodesic distances on meshes." in 2015. (author: Xu, Chunxu et al.) 


#ifndef FWP_QUEUE_H
#define FWP_QUEUE_H

//#define OUTPUTINFO

#include <vector>
#include "MyList.h"
#include "ListGroup.h"
#include "geodesic_memory.h"
#include <cstdio>

template<typename T>
class FWP_Queue
{
public:

	struct BinItem
	{
		T element;
		unsigned index;
		unsigned birthTime;
	};

	FWP_Queue();
	FWP_Queue(int w, int h);
	~FWP_Queue();

public:
	void setParametersToInitial();
	void setParameters1(double _binWidth, unsigned _binNum);
	void setParameters2(unsigned _Kmin, unsigned _Kmax, double _step);

	void setBinRange();

	void setPhase(unsigned _phase);
	unsigned getPhase() {return phase;}

	void* push(const T &elem, double key);
	T top();
	T top_for_parallel();

	void updatePhase();
	T pop();
	T pop_for_parallel();

	void remove(void *ptr);
	void remove_without_record(void *ptr);

	int size();
	bool empty();

//private:
	void updateParamenters();

private:
// 	std::vector< mylist<BinItem*> > bins;
// 	std::vector<unsigned> binSizes;
	ListGroup<BinItem> bins;

	//geodesic::MemoryAllocator<BinItem> m_memory_bin;

	unsigned totalSize;

	unsigned firstBin;				//first bin that is not empty;
	unsigned lastBin;				//last bin that is not empty;
	unsigned firstNonEmptyBin;
	unsigned midNonEmptyBin;
	double binWidth;

	int K;
	double threshold;
	int Kmin, Kmax;
	int smallCount, largeCount;
	double step;
	unsigned phase;

	int rmCounter1, rmCounter2;
	int popCounter;
	int estimateNum;

	BinItem binItem;
	T tmpElem;
};

template <typename T>
FWP_Queue<T>::FWP_Queue()
{
	phase = 0;
	largeCount = 0; smallCount = 0;
	totalSize = 0;
	rmCounter1 = 0, rmCounter2 = 0;
	popCounter = 0;
	estimateNum = 0;
}

template <typename T>
FWP_Queue<T>::FWP_Queue(int w, int h) : 
bins(w, h)
{
	phase = 0;
	largeCount = 0; smallCount = 0;
	totalSize = 0;
	rmCounter1 = 0, rmCounter2 = 0;
	popCounter = 0;
	estimateNum = 0;
}

template <typename T>
FWP_Queue<T>::~FWP_Queue()
{
	
}

template<typename T>
void FWP_Queue<T>::setParametersToInitial()
{
	//make the firstNonEmptyBin = firstBin
	//phase++;

	while (bins.empty(firstBin) && firstBin < lastBin) ++firstBin;
	while (bins.empty(lastBin - 1) && lastBin > firstBin) --lastBin;

	//update firstNonEmptyBin, midNonEmptyBin
	firstNonEmptyBin = firstBin;
	midNonEmptyBin = lastBin;

}

template <typename T>
void FWP_Queue<T>::setParameters1(double _binWidth, unsigned _binNum)
{
	binWidth = _binWidth;

	bins.resize(_binNum);

	firstBin = -1;
	lastBin = -1;
	firstNonEmptyBin = -1;
	midNonEmptyBin = -1;
}

template <typename T>
void FWP_Queue<T>::setParameters2(unsigned _Kmin, unsigned _Kmax, double _step)
{
	K = totalSize;
	Kmin = _Kmin; Kmax = _Kmax;
	step = _step;
}

template <typename T>
void FWP_Queue<T>::setBinRange()
{
	unsigned skipIntervals = 0;
	firstNonEmptyBin = firstBin;
	for (midNonEmptyBin = firstBin; midNonEmptyBin < lastBin; ++midNonEmptyBin)
	{
		skipIntervals += bins.sizesOf(midNonEmptyBin);
		if (skipIntervals >= K) 
		{
			++midNonEmptyBin;
			break;
		}
	}
	//update threshold
	threshold = midNonEmptyBin * binWidth;
}

template <typename T>
void FWP_Queue<T>::setPhase(unsigned _phase)
{
	phase = _phase;
}

template <typename T>
void* FWP_Queue<T>::push(const T &elem, double key)
{
	
	binItem.birthTime = phase;
	binItem.element = elem;
	binItem.index = (unsigned)(key / binWidth);

	void *rtPtr = bins.push_back(binItem, binItem.index);

	++totalSize;

	//update largeCount or smallCount
	if (key > threshold) ++largeCount;
	else ++smallCount;

	//update first & lastBin (if it need to be updated)
	if (lastBin == -1 || binItem.index >= lastBin) lastBin = binItem.index + 1;
	if (firstBin == -1 || binItem.index < firstBin) firstBin = binItem.index;

	if (firstNonEmptyBin == -1) firstNonEmptyBin = binItem.index;
	if (midNonEmptyBin == -1) midNonEmptyBin = binItem.index + 1;

	return rtPtr;
}

template <typename T>
T FWP_Queue<T>::top()
{
	while (firstNonEmptyBin < midNonEmptyBin && 
		(bins.empty(firstNonEmptyBin) || ((ListGroup<BinItem>::ListItem*)bins.begin(firstNonEmptyBin))->element.birthTime == phase))
	{
		++firstNonEmptyBin;
	}
	if (firstNonEmptyBin >= midNonEmptyBin) updateParamenters();

	void *headOfBinDummy = bins.begin(firstNonEmptyBin);
	return ((ListGroup<BinItem>::ListItem*)headOfBinDummy)->element.element;
}

template <typename T>
T FWP_Queue<T>::top_for_parallel()
{
	void *headOfBinDummy = bins.begin(firstNonEmptyBin);
	return ((ListGroup<BinItem>::ListItem*)headOfBinDummy)->element.element;
}

template <typename T>
void FWP_Queue<T>::updatePhase()
{
	//bool updated = false;
	//int update = 0;
	while (firstNonEmptyBin < midNonEmptyBin &&
		(bins.empty(firstNonEmptyBin) || ((ListGroup<BinItem>::ListItem*)bins.begin(firstNonEmptyBin))->element.birthTime == phase))
	{
		++firstNonEmptyBin;
		//updated = true;
	}
	//if (updated)phase++;
	if (firstNonEmptyBin >= midNonEmptyBin) phase++;
	//if (firstNonEmptyBin >= midNonEmptyBin) updateParamenters();
}

template <typename T>
T FWP_Queue<T>::pop()
{
#ifdef OUTPUTINFO
	++popCounter;
#endif
	while (firstNonEmptyBin < midNonEmptyBin && 
		(bins.empty(firstNonEmptyBin) || ((ListGroup<BinItem>::ListItem*)bins.begin(firstNonEmptyBin))->element.birthTime == phase))
	{
		++firstNonEmptyBin;
	}
	if (firstNonEmptyBin >= midNonEmptyBin) updateParamenters();

	void *headOfBinDummy = NULL;

	headOfBinDummy = bins.begin(firstNonEmptyBin);
	binItem = ((ListGroup<BinItem>::ListItem*)headOfBinDummy)->element;

	--totalSize;

	tmpElem = binItem.element;
	bins.erase(headOfBinDummy, binItem.index);

	return tmpElem;
}

template <typename T>
T FWP_Queue<T>::pop_for_parallel()
{
	void *headOfBinDummy = NULL;

	headOfBinDummy = bins.begin(firstNonEmptyBin);
	binItem = ((ListGroup<BinItem>::ListItem*)headOfBinDummy)->element;

	--totalSize;

	tmpElem = binItem.element;
	bins.erase(headOfBinDummy, binItem.index);

	return tmpElem;
}


template <typename T>
void FWP_Queue<T>::remove(void *ptr)
{
	binItem = ((ListGroup<BinItem>::ListItem*)ptr)->element;
	unsigned index = binItem.index;
	
	if (index >= midNonEmptyBin)
		largeCount--;
	else
		smallCount--;

#ifdef OUTPUTINFO
	if (binItem.index >= firstNonEmptyBin 
		&& binItem.index < midNonEmptyBin
		&& binItem.birthTime != phase) 
		++rmCounter2;
	else ++rmCounter1;
#endif

	bins.erase(ptr, index);
	--totalSize;
}

template <typename T>
void FWP_Queue<T>::remove_without_record(void *ptr)
{
	binItem = ((ListGroup<BinItem>::ListItem*)ptr)->element;
	unsigned index = binItem.index;
	bins.erase(ptr, index);
	--totalSize;
}

template <typename T>
int FWP_Queue<T>::size()
{
	return totalSize;
}

template <typename T>
bool FWP_Queue<T>::empty()
{
	return totalSize == 0;
}

template <typename T>
inline void FWP_Queue<T>::updateParamenters()
{
	//if it runs here, it means that the K need to be updated.
	//update K

#ifdef OUTPUTINFO
	printf("K: %d, estimateNum: %d, actualTouched: %d, actual poped: %d, removePhase: %d, removeNotPhase: %d, size: %d\n", 
		K, estimateNum, popCounter + rmCounter2, popCounter, rmCounter1, rmCounter2, totalSize);
	popCounter = 0; rmCounter1 = 0; rmCounter2 = 0;
#endif

	//update K
	K += step  * (largeCount - smallCount);
	if (K < Kmin) K = Kmin;
	if (K > Kmax) K = Kmax;

	//update large&smallCount
	largeCount = 0; smallCount = 0;

	//update phase
	++phase;

	while (bins.empty(firstBin) && firstBin < lastBin) ++firstBin;
	while (bins.empty(lastBin - 1) && lastBin > firstBin) --lastBin;


	//update firstNonEmptyBin, midNonEmptyBin
	firstNonEmptyBin = firstBin;

	unsigned skipIntervals = bins.sizesOf(firstNonEmptyBin);
	for (midNonEmptyBin = firstNonEmptyBin + 1; midNonEmptyBin < lastBin; ++midNonEmptyBin)
	{
		skipIntervals += bins.sizesOf(midNonEmptyBin);
		if (skipIntervals >= K) 
		{
			++midNonEmptyBin;
			break;
		}
	}
	estimateNum = skipIntervals;

	//update threshold
	threshold = midNonEmptyBin * binWidth;

	
}

#endif
