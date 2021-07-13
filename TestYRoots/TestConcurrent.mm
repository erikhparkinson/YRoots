//
//  TestConcurrent.m
//  TestYRoots
//
//  Created by Erik Hales Parkinson on 7/29/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#import <XCTest/XCTest.h>
#include <random>
#include "Utilities/MultiPool.hpp"
#include "Utilities/ConcurrentStack.hpp"

@interface TestConcurrent : XCTestCase

@end

@implementation TestConcurrent

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
    Function::clearSavedFunctions();
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

uint64_t m_timingTrials = 1000000;

struct tempStruct {
    size_t a;
    double b;
    double c;
    std::vector<uint64_t> d;
    std::string e;
    tempStruct() : a(0), b(0), c(0) {}
};

- (void)testMultiPool {
    tempStruct defualtStruct;
    MultiPool<tempStruct> myPools;
    
    int numPools = 10;
    size_t numPops = 1000000;
    
    for(int i = 0; i < numPools; i++) {
        myPools.emplace_back(defualtStruct, 1);
    }
    
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, numPools-1); // define the range

    for(size_t pop=0; pop<numPops; pop++) {
        int poolNum  = distr(gen);
        tempStruct* temp = myPools[poolNum].pop();
        (temp->a)++;
        for (size_t j = 0; j < 4; j++) {
            (temp->e) += '0';
        }
        (temp->e) += '1';
        
        poolNum  = distr(gen);
        myPools[poolNum].push(temp);
    }
    
    //Everything in the pools should have a string of 4 zeros then a 1 repeated a times.
    //Otherwise something was being modified by two objects at once.
    for(size_t poolNum = 0 ; poolNum < myPools.size(); poolNum++) {
        for(size_t popNum = 0; popNum < myPools[poolNum].size(); popNum++) {
            tempStruct* temp = myPools[poolNum].pop();
            XCTAssert((temp->a)*5 == (temp->e).size());
            size_t spot = 0;
            for (size_t i = 0; i < temp->a; i++) {
                for (size_t j = 0; j < 4; j++) {
                    XCTAssert((temp->e)[spot] == '0');
                    spot++;
                }
                XCTAssert((temp->e)[spot] == '1');
                spot++;
            }
            myPools[poolNum].push(temp);
        }
    }
        
    //Time it
    tempStruct* temp;
    ObjectPool<tempStruct>& timePool = myPools[0];
    
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < m_timingTrials; i++) {
        temp = timePool.pop();
        timePool.push(temp);
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

    double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    
    std::cout << "\nMulti Pool push pop takes " << formatTimePretty(nanos/(m_timingTrials)) << ".\n\n";
}

- (void)testStack {
    size_t numThreads = 1;

    tempStruct defualtStruct;
    MultiPool<tempStruct> myPools;
    ConcurrentStack<tempStruct> myStack(numThreads);
    
    for(int i = 0; i < numThreads; i++) {
        myPools.emplace_back(defualtStruct, 2);
    }
    
    for (size_t i = 0; i < 100; i++) {
        tempStruct* temp = myPools[0].pop();
        (temp->a)++;
        myStack.push(0, temp);
        tempStruct* temp2 = myStack.pop(0);
        tempStruct* temp3 = myStack.pop(0);
        if(temp2) {
            myStack.push(0, temp2);
        }
        if(temp3) {
            myStack.push(0, temp3);
        }
        tempStruct* temp4 = myStack.pop(0);
        if(temp4) {
            (temp4->a)--;
            myPools[0].push(temp4);
        }
    }

    auto threadFunc = [](size_t threadNum, size_t n, ObjectPool<tempStruct>& pool, ConcurrentStack<tempStruct>& stack) {
        for (size_t i = 0; i < n; i++) {
            tempStruct* temp = pool.pop();
            (temp->a)++;
            stack.push(threadNum, temp);
            tempStruct* temp2 = stack.pop(threadNum);
            tempStruct* temp3 = stack.pop(threadNum);
            if(temp2) {
                stack.push(threadNum, temp2);
            }
            if(temp3) {
                stack.push(threadNum, temp3);
            }
            tempStruct* temp4 = stack.pop(threadNum);
            if(temp4) {
                (temp4->a)--;
                pool.push(temp4);
            }
        }
    };
    
    std::vector<std::unique_ptr<std::thread> > threadPool;
    for(size_t i = 0 ; i< numThreads; i++) {
        threadPool.emplace_back(::make_unique<std::thread>(threadFunc, i, 1000, std::ref(myPools[i]), std::ref(myStack)));
    }
    for(size_t i = 0 ; i< numThreads; i++) {
        threadPool[i]->join();
    }
    
    //Pop everything else off the stack
    tempStruct* temp5;
    while((temp5 = myStack.pop(0))) {
        XCTAssert(temp5->a == 1);
        (temp5->a)--;
        myPools[0].push(temp5);
    }
    
    //Everything in the pool should have a = 0. Otherwise it was lost.
    for(size_t poolNum = 0; poolNum < myPools.size(); poolNum++) {
        for(size_t i = 0 ; i < myPools[poolNum].size(); i++) {
            tempStruct* temp = myPools[poolNum].pop();
            XCTAssert(temp->a == 0);
            myPools[poolNum].push(temp);
        }
    }
    
    //Time it
    tempStruct* temp = myPools[0].pop();
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < m_timingTrials; i++) {
        myStack.push(0, temp);
        temp = myStack.pop(0);
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

    myPools[0].push(temp);
    double nanos = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    
    std::cout << "\nConcurrent Stack push pop takes " << formatTimePretty(nanos/(m_timingTrials)) << ".\n\n";
}

@end
