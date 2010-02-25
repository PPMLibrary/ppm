// 
#ifndef LIST_H
#define LIST_H
#include "Boolean.h"
#include <assert.h>
#include <iostream>

class ListNode {
    friend class List;
    friend class ListIterator;
private:
    ListNode() {next = 0;}	// Private, since only needed by List
    int item;
    ListNode  *next;
};

class List {
    friend class ListIterator;
public:
    List();
    ~List();
    List(const List &);
    const List & operator=(const List &);

    Boolean isEmpty() const;
    int length() const;
    void insertNth(int, int);
    void Append(int);
    void Prepend(int item){insertNth(1, item);}
    void deleteNth(int);
    int nth (int) const;

private:
    int size;
    ListNode  *head, *tail;
    inline ListNode *ptrToNth (int) const;
};

class ListIterator {
public:
    ListIterator(){current = NULL;}
    int operator()(){return current->item;}
    void start(const List & l) {current = l.head;}
    void operator++(){current = current->next;}
    Boolean done() {return current == 0 ? TRUE : FALSE;}
private:
    ListNode *current;
    Boolean JustStarted;
};

/*---------------------implementation-------------------------*/

List::List() {head = 0; tail = 0; size=0;}

List::List(const List &old) {
    head = 0;
    tail = 0;
    size = 0;
    *this = old;
}
  
const List & List::operator=(const List & old)
{
    if (&old == this) {}
	// Self assignment; do nothing
    else {
	    while (!isEmpty()) {
	        deleteNth (1);
        }

        size = old.size;
        head = 0;
        tail = 0;

        if (old.head != 0) {
	        head = new ListNode;
	        assert (head != 0);
			ListNode *newPtr, *oldPtr;
	        newPtr = head;
	        oldPtr = old.head;
	        newPtr->item = oldPtr->item;

	        while (oldPtr->next != 0) {
		        newPtr->next = new ListNode;
		        newPtr = newPtr->next;
		        oldPtr = oldPtr->next;
		        newPtr->item = oldPtr->item;
	        }

	        newPtr->next = 0;
	        tail = newPtr;
	        newPtr = 0;
	        oldPtr = 0;
	    }
    }
    return *this;
}

List::~List()
{
    while (!isEmpty()) {
		deleteNth (1);
    }
}

Boolean List::isEmpty() const {return head == 0 ? TRUE : FALSE;}

int List::length() const {return size;}

void List::insertNth (int position, int item)
{
    assert (position >= 1 && position <= size+1);
    size++;
    ListNode *newPtr = new ListNode;
    newPtr->item = item;
    if (position == 1) {
	    newPtr->next = head;
	    head = newPtr;
    } else {
	    ListNode *prev = ptrToNth (position-1);
	    newPtr->next = prev->next;
	    prev->next = newPtr;
    }
    if (tail -> next != 0)
        tail = tail -> next;
}

void List::Append (int item)
{
    size++;
    ListNode *newPtr = new ListNode;
    newPtr->item = item;
    
    if (isEmpty())
        tail = head = newPtr;
    else {
	    tail->next = newPtr;
	    tail = newPtr;
	    newPtr = 0;
    }
}

void List::deleteNth (int position)
{
    assert (1 <= position && position <= size);
    size--;
    ListNode *old;
    if (position == 1) {
	    old = head;
	    head = head->next;
    } else {
	    ListNode *prev = ptrToNth (position -1);
	    old = prev->next;

        if (old == tail)
	        tail = prev;

	    prev->next = old->next;
    }
    delete old;
}

int List::nth (int n) const {return ptrToNth(n)->item;}

inline ListNode *List::ptrToNth(int position) const
{
    assert (position >= 1 && position <= size);
    ListNode *p = head;
    for (int i=1; i<position; i++)
	    p = p->next;
    return p;
}

#endif
