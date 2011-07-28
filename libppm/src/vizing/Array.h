// Templated safe Array.  An Array are resizable by calling */
// resize(newsize), making the array smaller or larger.  If the array */
// is made smaller, only a partial copy is performed.  All arrays */
// initialize to 0, so the Templated class must have a way to */
// initialize to 0 by writing george = 0;
//
// Bounds checking is done only via assert() calls.
#ifndef ARRAY_H
#define ARRAY_H

#include "Boolean.h"
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include "edge.h"
using namespace std;

template<class T>
class Array {
public:
    Array(int size = 0);
    Array(const Array<T> & source);
    ~Array();
    void resize(const int newsize);       
    // OK to resize the array to make it smaller;
    int getSize() const;
    T &operator[](int);
    T &operator()(int) const;       
    // Just like operator[], but cannot be an lvalue;
    void Append(T);
    void insert(int, T);

private:
    T *ptr; // pointer to first element of array
    int size; // size of the array
    void getNewArray (int size, Boolean initializeToZero);
};

/*---------------------implementation-------------------------*/

template<class T>
Array<T>::Array(int arraySize)
{
    getNewArray (arraySize, TRUE);
}

template<class T>
Array<T>::Array(const Array<T> & source)
{
    getNewArray (source.size, FALSE);
    for (int i = 0; i < size; i++)
	ptr[i] = source.ptr[i];
}

template<class T>
void Array<T>::insert(int i, T item)
{
    assert( i < size );
    ptr[i] = item;
}

template<class T>
Array<T>::~Array()
{
    delete [] ptr;
    ptr = NULL;
}

template<class T>
void Array<T>::resize(const int newSize) {
    T *oldPtr = ptr;
    int oldSize = size;
    getNewArray (newSize, TRUE);
    for (int i = 0; i < newSize && i < oldSize; i++)
	    ptr[i] = oldPtr[i];
    delete [] oldPtr;
}

template<class T>
void Array<T>::Append(T item) {
    int newSize = size + 1;
    resize(newSize);
    ptr[newSize-1] = item;
}

template<class T>
int Array<T>::getSize() const {
    return size;
}

template<class T>
T &Array<T>::operator[](int subscript)
{
    // check for subscript out of range error
    assert(0 <= subscript && subscript < size);
    return ptr[subscript];
}

// operator() is just like operator[] 
// but is const'ed and cannot be an lvalue
template<class T>
T &Array<T>::operator()(int subscript) const
{
    // check for subscript out of range error
    assert(0 <= subscript && subscript < size);
    return ptr[subscript];
}

template<class T>
void Array<T>::getNewArray (int arraySize, Boolean initializeToZero)
{
    size = arraySize;
    ptr = new T[size];
    assert(ptr != 0);
}

#endif
