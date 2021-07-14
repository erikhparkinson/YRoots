//
//  ConcurrentStack.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 7/29/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef ConcurrentStack_h
#define ConcurrentStack_h

#include <atomic>
#include "Utilities/MultiPool.hpp"

template<typename T>
struct StackNode {
    StackNode<T>* prev;
    T* value;
    
    StackNode() : prev(nullptr), value(nullptr) {}
};

template<typename T>
class ConcurrentStack {
public:
    ConcurrentStack(size_t _numThreads)
    {
        StackNode<T> defaultNode;
        for(size_t i = 0; i < _numThreads; i++) {
            m_pools.emplace_back(defaultNode, 1024);
        }
        
        StackNode<T>* root = m_pools[0].pop();
        root->value = nullptr;
        root->prev = root;
        m_top.store(root);
    }
        
    T* pop(size_t threadNum) {
        //Returns a nullptr if the stack is empty.
        T* ptr;
        StackNode<T>* oldTop = m_top.load(std::memory_order_relaxed);
        while(!m_top.compare_exchange_weak(oldTop, oldTop->prev, std::memory_order_release, std::memory_order_relaxed));
        //Push back into the pool is it's not null
        if(likely(ptr = oldTop->value)) {
            m_pools[threadNum].push(oldTop);
        }
        return ptr;
    }
    
    void push(size_t threadNum, T* ptr) {
        StackNode<T>* newNode = m_pools[threadNum].pop();
        newNode->value = ptr;
        newNode->prev = m_top.load(std::memory_order_relaxed);
        while (!m_top.compare_exchange_weak(newNode->prev, newNode, std::memory_order_release, std::memory_order_relaxed));
    }
    
private:
    ConcurrentStack<T>(const ConcurrentStack<T>&) =delete;
    ConcurrentStack<T>& operator=(const ConcurrentStack<T>&) =delete;
    
private:
    std::atomic<StackNode<T>*>          m_top;
    MultiPool<StackNode<T> >             m_pools;
};

#endif /* ConcurrentStack_h */
