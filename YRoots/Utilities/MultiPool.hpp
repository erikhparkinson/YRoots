//
//  concurrentPool.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/29/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef MultiPool_h
#define MultiPool_h

#include <atomic>
#include <vector>
#include "../Utilities/macros.hpp"

//This pool is only ever accessed by a single thread, but it can
//have objects that were originaly in another pool pushed into it
//and objects from this pool can be pushed into other pools.
template<typename T>
class ObjectPool {
public:
    ObjectPool(T _defaultObject, size_t _initialSize) :
    m_defaultObject(_defaultObject),
    m_head(0),
    m_tail(0),
    m_numPointers(0),
    m_mod(0) //Should always be m_numPointers - 1
    {
        //_initialSize must be power of 2!
        
        //Create more pointers
        m_pointers.resize(_initialSize, nullptr);
        m_numPointers = _initialSize;
        m_mod = m_numPointers - 1;

        //Create more objects
        m_objects.emplace_back(m_numPointers, m_defaultObject);
        size_t temp = m_objects.size() - 1;
        for(size_t i = 0; i < m_numPointers; i++) {
            m_pointers[i] = &m_objects[temp][i];
        }
        
        //Reset head and tail
        m_head = 0;
        m_tail = m_numPointers;
    }

    T* pop() {
        if(unlikely(m_tail == m_head)){
            addObjects();
        }
        return m_pointers[(m_head++)&m_mod];
    }
    
    void push(T* ptr) {
        if(unlikely(m_tail == m_head + m_numPointers)){
            addPointers();
        }
        m_pointers[(m_tail++)&m_mod] = ptr;
    }
    
    size_t size() {
        return m_numPointers;
    }
    
private:
    void addObjects() {
        //There are no pointers left in the original vector.
        assert(m_head == m_tail);
        
        //Create more objects
        m_objects.emplace_back(m_numPointers, m_defaultObject);
        size_t temp = m_objects.size() - 1;
        for(size_t i = 0; i < m_numPointers; i++) {
            m_pointers[i] = &m_objects[temp][i];
        }
        
        //Reset head and tail
        m_head = 0;
        m_tail = m_numPointers;
    }
    
    void addPointers() {
        //There is no room to push a ptr into the vector.
        assert(m_head + m_numPointers == m_tail);
        
        //Create more pointers
        m_head = 0;
        m_tail = m_numPointers;
        m_numPointers *= 2;
        m_mod = m_numPointers - 1;
        m_pointers.resize(m_numPointers, nullptr);
    }
    
private:
    T                               m_defaultObject;
    std::vector<std::vector<T> >     m_objects;          //The pool of objects
    std::vector<T*>                 m_pointers;
    uint64_t                        m_head;             //The next place to pop from.
    uint64_t                        m_tail;             //Next place to push to
    uint64_t                        m_numPointers;
    uint64_t                        m_mod;
};

template<typename T>
using MultiPool = std::vector<ObjectPool<T> >;

#endif /* MultiPool_h */
