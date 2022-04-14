#ifndef MY_LIST_H
#define MY_LIST_H

template<typename T>
class mylist
{
public:
	struct ListItem
	{
		T element;
		ListItem *next;
		ListItem *prev;
	};

	mylist();
	mylist(mylist<T> &_list);
	~mylist();

	void* insert(void *pos, const T& val);
	void erase(void *pos);

	void* push_back(const T &val);
	void* push_front(const T &val);

	T top();
	T pop_front();

	void* begin();
	void* end();

	void* next(void* iter);
	void* prev(void* iter);

	T getElement(void* iter);

	int size();
	bool empty();

private:
	ListItem *head;
	ListItem *tail;
	int _size;
};

template<typename T>
mylist<T>::mylist()
{
	tail = new ListItem;
	//tail = (ListItem*)malloc(sizeof(ListItem));
	head = tail;

	tail->next = 0;
	tail->prev = 0;
	_size = 0;
}

template<typename T>
mylist<T>::mylist(mylist<T> &_list)
{
	ListItem *_head = (ListItem*)(_list.begin());
	ListItem *_tail = (ListItem*)(_list.end());

	tail = new ListItem;
	//tail = (ListItem*)malloc(sizeof(ListItem));
	head = tail;
	tail->next = 0;
	tail->prev = 0;
	_size = 0;

	for (; _head != _tail; _head = _head->next)
	{
		push_back(_head->element);
	}
}

template<typename T>
mylist<T>::~mylist()
{
	while(head != tail)
	{
		ListItem *next = head->next;
		//free(head);
		delete next;
		head = next;
	}
	//free(head);
	delete head;
}

template<typename T>
inline void* mylist<T>::insert(void *pos, const T& val)
{
	if(!pos) return NULL;
	ListItem *position = (ListItem*)pos;

	ListItem *curItem = new ListItem;
	//ListItem *curItem = (ListItem*)malloc(sizeof(ListItem));
	curItem->element = val;

	curItem->next = position;
	curItem->prev = position->prev;

	if(position->prev) position->prev->next = curItem;
	else head = curItem;

	position->prev = curItem;

	++_size;

	return curItem;
}

template<typename T>
inline void mylist<T>::erase(void *pos)
{
	if(!pos) return;
	ListItem *position = (ListItem*)pos;
	if(position == tail) return;

	if(position->prev)
		position->prev->next = position->next;
	else
		head = position->next;
	position->next->prev = position->prev;

	//free(position);
	delete position;

	--_size;
}

template<typename T>
inline void* mylist<T>::push_back(const T &val)
{
	return insert(end(), val);
}

template<typename T>
inline void* mylist<T>::push_front(const T &val)
{
	return insert(begin(), val);
}

template<typename T>
inline T mylist<T>::top()
{
	return head->element;
}

template<typename T>
inline T mylist<T>::pop_front()
{
	T rt = head->element;
	ListItem *rmPtr = head;
	head = head->next;
	//no need to judge whether head is NULL
	head->prev = NULL;
	--_size;
	delete rmPtr;
	return rt;
}

template<typename T>
inline void* mylist<T>::begin()
{
	return (void*)head;
}

template<typename T>
inline void* mylist<T>::end()
{
	return (void*)tail;
}

template<typename T>
inline void* mylist<T>::next(void* iter)
{
	if(!iter) return 0;
	ListItem *i = (ListItem*)iter;
	i = i->next;
	return (void*)i;
}

template<typename T>
inline void* mylist<T>::prev(void* iter)
{
	if(!iter) return 0;
	ListItem *i = (ListItem*)iter;
	i = i->prev;
	return (void*)i;
}

template<typename T>
inline T mylist<T>::getElement(void* iter)
{
	if(!iter) return 0;
	ListItem *i = (ListItem*)iter;
	return i->element;
}

template<typename T>
inline int mylist<T>::size()
{
	return _size;
}

template<typename T>
inline bool mylist<T>::empty()
{
	return _size == 0;
}

#endif