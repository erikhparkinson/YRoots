//
//  InputFileParser.h
//  YRoots
//
//  Created by Erik Hales Parkinson on 5/23/20.
//  Copyright © 2020 Erik Hales Parkinson. All rights reserved.
//

#ifndef InputFileParser_h
#define InputFileParser_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <thread>
#include <algorithm>
#include "Function.h"
#include "Timer.h"

class InputFileParser {
public:
    InputFileParser(std::string _filename):
    m_filename(_filename)
    {
        m_timer.registerTimer(m_timerInputParserIndex, "Input Parser");
    }
    
    void parse() {
        m_timer.startTimer(m_timerInputParserIndex);

        //Read the whole file into a string
        std::ifstream inputFile;
        std::string line;
        std::string inputString;
        inputFile.open(m_filename);
        while(std::getline(inputFile, line)) {
            inputString += line;
        }
        inputFile.close();
        
        //Remove whitespace
        inputString.erase(remove_if(inputString.begin(), inputString.end(), isspace), inputString.end());
        
        //Split by colons, remove everything after the last colon
        std::vector<std::string> lines = split(inputString, ";");
        lines.pop_back();
        
        //Make sure the end line is there and remove it
        if(lines.size()  == 0 || lines[lines.size() - 1] != "END") {
            printAndThrowRuntimeError("Incorrect input format! END not found!");
        }
        lines.pop_back();
        size_t parseSpot = 0;
        
        //Parse the PARAMETERS
        if(lines.size() > parseSpot && lines[parseSpot] == "PARAMETERS") {
            parseSpot++;
            parseParameters(lines, parseSpot);
        }
        
        //Parse the INTERVAL
        if(lines.size() > parseSpot && lines[parseSpot] == "INTERVAL") {
            parseSpot++;
            parseInterval(lines, parseSpot);
        }
        else {
            printAndThrowRuntimeError("Parser Error! No INTERVAL Found!");
        }

        //Parse the FUNCTIONS
        if(lines.size() > parseSpot && lines[parseSpot] == "FUNCTIONS") {
            parseSpot++;
            parseFunctions(lines, parseSpot);
        }
        else {
            printAndThrowRuntimeError("Parser Error! No FUNCTIONS Found!");
        }
        
        //Make sure the functions size and interval size matches!
        if(m_functions[0].size() != m_interval.lowerBounds.size()) {
            printAndThrowRuntimeError("Parser Error! Number of Functions must match interval dimension!");
        }
        
        m_timer.stopTimer(m_timerInputParserIndex);
    }
    
    //Getters
    std::vector<std::vector<Function::SharedFunctionPtr>>& getFunctions() {
        return m_functions;
    }
    
    Interval& getInterval() {
        return m_interval;
    }
        
    const SubdivisionParameters& getSubdivisionParameters() const {
        return m_subdivisionParameters;
    }
    
    const GeneralParameters& getGeneralParameters() const {
        return m_generalParameters;
    }
    
private:
    bool parseBool(std::string inputString) {
        std::transform(inputString.begin(), inputString.end(), inputString.begin(), ::tolower);
        if(inputString == "true" || inputString == "t" || inputString == "y" || inputString == "yes") {
            return true;
        }
        else if (inputString == "false" || inputString == "f" || inputString == "n" || inputString == "no") {
            return false;
        }
        printAndThrowRuntimeError("Parser Error! Invalid bool format!");
        return false;
    }
    
    int parseInteger(const std::string& inputString) {
        try{
            return std::stoi(inputString);
        }
        catch(...) {
            printAndThrowRuntimeError("Parser Error! Invalid integer format!");
        }
        return 0;
    }
    
    double parseDouble(const std::string& inputString) {
        try{
            return std::stod(inputString);
        }
        catch(...) {
            printAndThrowRuntimeError("Parser Error! Invalid integer format!");
        }
        return 0;
    }
    
    void parseParameters(const std::vector<std::string>& lines, size_t& parseSpot){
        while(parseSpot < lines.size()) {
            //Check for the END
            if(lines[parseSpot] == "PARAMETERS_END") {
                parseSpot++;
                return;
            }
            
            //Parse the Parameters
            std::vector<std::string> parameterString = split(lines[parseSpot], "=");
            if(parameterString.size() != 2) {
                printAndThrowRuntimeError("Parser Error! Invalid Parameter Format!");
            }
            else if(parameterString[0] == "numThreads") {
                int numThreads = parseInteger(parameterString[1]);
                if(numThreads == -1) {
                    m_generalParameters.numThreads = (size_t)std::thread::hardware_concurrency();
                }
                else if (numThreads <= 0) {
                    printAndThrowRuntimeError("Parser Error! Invalid numThreads " + std::to_string(numThreads) + "!");
                }
                else {
                    m_generalParameters.numThreads = numThreads;
                }
            }
            else if(parameterString[0] == "relApproxTol") {
                m_subdivisionParameters.relApproxTol = parseConstantNum(parameterString[1]);
                if(m_subdivisionParameters.relApproxTol < 0) {
                    printAndThrowRuntimeError("Parameter Error! relApproxTol must be >= 0!");
                }
            }
            else if(parameterString[0] == "absApproxTol") {
                m_subdivisionParameters.absApproxTol = parseConstantNum(parameterString[1]);
                if(m_subdivisionParameters.absApproxTol < 0) {
                    printAndThrowRuntimeError("Parameter Error! relApproxTol must be >= 0!");
                }
            }
            else if(parameterString[0] == "goodZerosFactor") {
                m_subdivisionParameters.goodZerosFactor = parseConstantNum(parameterString[1]);
                if(m_subdivisionParameters.goodZerosFactor < 0) {
                    printAndThrowRuntimeError("Parameter Error! goodZerosFactor must be >= 0!");
                }
            }
            else if(parameterString[0] == "minGoodZerosTol") {
                m_subdivisionParameters.minGoodZerosTol = parseConstantNum(parameterString[1]);
                if(m_subdivisionParameters.minGoodZerosTol < 0) {
                    printAndThrowRuntimeError("Parameter Error! minGoodZerosTol must be >= 0!");
                }
            }
            else if(parameterString[0] == "approximationDegree") {
                m_subdivisionParameters.approximationDegree = parseInteger(parameterString[1]);
                if(m_subdivisionParameters.approximationDegree < 1) {
                    printAndThrowRuntimeError("Parameter Error! minGoodZerosTol must be >= 1!");
                }
            }
            else if(parameterString[0] == "maxLevel") {
                m_subdivisionParameters.maxLevel = parseInteger(parameterString[1]);
                if(m_subdivisionParameters.maxLevel < 0) {
                    printAndThrowRuntimeError("Parameter Error! minGoodZerosTol must be >= 0!");
                }
            }
            else if(parameterString[0] == "trackIntervals") {
                m_generalParameters.trackIntervals = parseBool(parameterString[1]);
            }
            else if(parameterString[0] == "trackProgress") {
                m_generalParameters.trackProgress = parseBool(parameterString[1]);
            }
            else if(parameterString[0] == "useTimer") {
                m_generalParameters.useTimer = parseBool(parameterString[1]);
            }
            else {
                printAndThrowRuntimeError("Parser Error! Unrecognized Parameter " + parameterString[0] + "!");
            }
            parseSpot++;
        }
        printAndThrowRuntimeError("Parser Error! No PARAMETERS_END Found!");
    }
    
    double parseConstantNum(const std::string& inputString) {
        try{
            return std::stod(inputString);
        }
        catch(...) {}
        
        std::vector<std::string> emptyVector;
        Function constantFunction("", inputString, emptyVector);
        if(constantFunction.getFunctionType() != FunctionType::CONSTANT) {
            printAndThrowRuntimeError("Parser Error! Invalid Constant Number!");
        }
        else {
            return constantFunction.getValue();
        }
        return 0.0;
    }

    void parseInterval(const std::vector<std::string>& lines, size_t& parseSpot){
        while(parseSpot < lines.size()) {
            //Check for the END
            if(lines[parseSpot] == "INTERVAL_END") {
                parseSpot++;
                return;
            }
            //Parse the interval
            const std::string& currIntervalString = lines[parseSpot];
            if(currIntervalString.length() < 2 || currIntervalString.front() != '[' || currIntervalString.back() != ']') {
                printAndThrowRuntimeError("Parser Error! Incorrect Interval Format!");
            }
            else {
                std::vector<std::string> intervalNums = split(currIntervalString.substr(1, currIntervalString.length()-2), ",");
                if(intervalNums.size() != 2) {
                    printAndThrowRuntimeError("Parser Error! Incorrect Interval Format! Should have 2 Numbers!");
                }
                else {
                    m_interval.lowerBounds.push_back(parseConstantNum(intervalNums[0]));
                    m_interval.upperBounds.push_back(parseConstantNum(intervalNums[1]));
                }
            }
            parseSpot++;
        }
        printAndThrowRuntimeError("Parser Error! No INTERVAL_END Found!");
    }

    void parseFunctions(const std::vector<std::string>& lines, size_t& parseSpot) {
        //Get the function names
        if(lines.size() <= parseSpot || lines[parseSpot].substr(0,8) != "function") {
            printAndThrowRuntimeError("Function definitions not found!");
        }
        std::vector<std::string> functionNames = split(lines[parseSpot].substr(8), ",");
        parseSpot++;
        
        //Get the variable names
        if(lines.size() <= parseSpot || lines[parseSpot].substr(0,14) != "variable_group") {
            printAndThrowRuntimeError("Function definition not found!");
        }
        std::vector<std::string> variableNames = split(lines[parseSpot].substr(14), ",");
        parseSpot++;
        
        //Make sure there are the same number of functions are variables
        if(variableNames.size() != functionNames.size()) {
            printAndThrowRuntimeError("Parser Error! #Functions != #Variables!");
        }
        
        //Get the individual functions
        bool foundEnd = false;
        m_functions.resize(m_generalParameters.numThreads);
        while(parseSpot < lines.size()) {
            if(lines[parseSpot] == "FUNCTIONS_END") {
                parseSpot++;
                foundEnd = true;
                break;
            }
            std::vector<std::string> functionData = split(lines[parseSpot], "=");
            if(functionData.size() != 2) {
                printAndThrowRuntimeError("Incorrect function definition format!");
            }
            Function::addFunction(functionData[0], functionData[1], variableNames);
            parseSpot++;
        }
        
        if(!foundEnd) {
            printAndThrowRuntimeError("Parser Error! No FUNCTIONS_END Found!");
        }
        
        //Duplicate the functions for all the threads
        Function::addThreadFunctions(m_generalParameters.numThreads);
        
        //Make sure we found definitions for all the functions
        for(size_t threadNum = 0; threadNum < m_generalParameters.numThreads; threadNum++) {
            for(size_t funcNum = 0; funcNum < functionNames.size(); funcNum++) {
                m_functions[threadNum].push_back(Function::getThreadFunctionByName(threadNum, functionNames[funcNum]));
            }
        }
    }
    
private:
    std::string m_filename;
    
    //Parsed Variables
    std::vector<std::vector<Function::SharedFunctionPtr>> m_functions;
    Interval m_interval;
    SubdivisionParameters m_subdivisionParameters;
    GeneralParameters     m_generalParameters;
    
    static size_t           m_timerInputParserIndex;
    Timer&                  m_timer = Timer::getInstance();
};

size_t InputFileParser::m_timerInputParserIndex = -1;

#endif /* InputFileParser_h */
