#ifndef LIST_GROUP_H
#define LIST_GROUP_H

#include <vector>
#include "MyList.h"
#include "geodesic_memory.h"

template<typename T>
class ListGroup
{
public:
	struct ListItem
	{
		T element;
		ListItem *next;
		ListItem *prev;
	};

	ListGroup();
	ListGroup(int w, int h);
	~ListGroup();

	void resize(int _listNum);

	void* push_back(const T &val, unsigned index);

	void* insert(void *pos, const T& val, unsigned index);
	void erase(void *ptr, unsigned index);

	void* begin(unsigned index);
	void* end(unsigned index);

	T pop();
	T top();

	int layers();

	int size();
	int sizesOf(unsigned index);
	bool empty();
	bool empty(unsigned index);

private:
// 	ListItem **heads;
// 	ListItem **tails;
// 	unsigned *_sizes;
	std::vector<ListItem*> heads;
	std::vector<ListItem*> tails;
	std::vector<unsigned> _sizes;

	unsigned totalSize;
	int listNum;

	geodesic::MemoryAllocator<ListItem> m_memory;
};

template<typename T>
ListGroup<T>::ListGroup()
{
	listNum = 0;
	totalSize = 0;
}

template<typename T>
ListGroup<T>::ListGroup(int w, int h) : 
m_memory(w, h)
{
	listNum = 0;
	totalSize = 0;
}

template<typename T>
ListGroup<T>::~ListGroup()
{
}

template<typename T>
void ListGroup<T>::resize(int _listNum)
{
	heads.resize(_listNum); tails.resize(_listNum); _sizes.resize(_listNum, 0);
	unsigned originalListNum = listNum;
	listNum = _listNum;

	for (unsigned i = originalListNum; i < listNum; ++i)
	{
		tails[i] = m_memory.allocate();
		heads[i] = tails[i];

		tails[i]->next = 0; tails[i]->prev = 0; _sizes[i] = 0;
		
	}
}

template<typename T>
inline void* ListGroup<T>::insert(void *pos, const T& val, unsigned index)
{
	if(!pos) return NULL;

	ListItem *position = (ListItem*)pos;

	ListItem *curItem = m_memory.allocate();
	curItem->element = val;

	curItem->next = position;
	curItem->prev = position->prev;

	if(position->prev) position->prev->next = curItem;
	else heads[index] = curItem;

	position->prev = curItem;

	++_sizes[index];
	++totalSize;

	return curItem;
}

template<typename T>
inline void ListGroup<T>::erase(void *pos, unsigned index)
{
	if(!pos) return;
	ListItem *position = (ListItem*)pos;
	if(position == tails[index]) return;

	if(position->prev)
		position->prev->next = position->next;
	else
		heads[index] = position->next;
	position->next->prev = position->prev;

	m_memory.deallocate(position);
	--_sizes[index];
	-- totalSize;
}

template<typename T>
inline void* ListGroup<T>::push_back(const T &val, unsigned index)
{
	return insert(end(index), val, index);
}

template<typename T>
inline void* ListGroup<T>::begin(unsigned index)
{
	if (index + 1 > listNum) resize(index+1);
	return (void*)heads[index];
}

template<typename T>
inline void* ListGroup<T>::end(unsigned index)
{
	if (index + 1 > listNum) resize(index+1);
	return (void*)tails[index];
}

template<typename T>
inline int ListGroup<T>::layers()
{
	return listNum;
}

template<typename T>
inline int ListGroup<T>::size()
{
	return totalSize;
}

template<typename T>
inline int ListGroup<T>::sizesOf(unsigned index)
{
	return _sizes[index];
}

template<typename T>
inline bool ListGroup<T>::empty()
{
	return totalSize == 0;
}

template<typename T>
inline bool ListGroup<T>::empty(unsigned index)
{
	return _sizes[index] == 0;
}

#endif
