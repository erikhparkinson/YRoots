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
#include "PowerBasisPolynomial.h"

class InputFileParser {
public:
    InputFileParser(std::string _filename):
    m_filename(_filename)
    {}
    
    std::vector<std::vector<std::unique_ptr<FunctionInterface>>> parseFunctions(size_t _numThreads) {
        std::vector<std::vector<std::unique_ptr<FunctionInterface>>> functions(_numThreads);
        
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
        
        //Make sure the end line is there
        if(lines.size()  == 0 || lines[lines.size() - 1] != "END") {
            throw std::runtime_error("Incorrect input format! END not found!");
        }
        
        //Get the function names
        if(lines.size() == 0 || lines[0].substr(0,8) != "function") {
            throw std::runtime_error("function definitions not found!");
        }
        std::set<std::string> functionNames;
        std::vector<std::string> functionNamesVector = split(lines[0].substr(8), ",");
        for(size_t i = 0; i < functionNamesVector.size(); i++) {
            functionNames.insert(functionNamesVector[i]);
        }
        
        //Get the variable names
        if(lines.size() == 1 || lines[1].substr(0,14) != "variable_group") {
            throw std::runtime_error("function definition not found!");
        }
        std::vector<std::string> variableNames = split(lines[1].substr(14), ",");
        
        //Get the individual functions
        for(size_t lineNum = 2; lineNum  + 1 < lines.size(); lineNum++) {
            std::vector<std::string> functionData = split(lines[lineNum], "=");
            if(functionData.size() != 2) {
                throw std::runtime_error("Incorrect function definition format!");
            }
            std::string functionName = functionData[0];
            std::set<std::string>::iterator nameFind = functionNames.find(functionName);
            if (nameFind == functionNames.end()) {
                throw std::runtime_error("Subfunctions not implemented!");
            }
            else {
                for(size_t i = 0 ; i < _numThreads; i++) {
                    functions[i].push_back(std::make_unique<PowerBasisPolynomial>(functionData[1], variableNames));
                }
                functionNames.erase(nameFind);
            }
        }
        if(functionNames.size() != 0) {
            for(std::set<std::string>::iterator it = functionNames.begin(); it != functionNames.end(); it++) {
                std::cout<<"No definition found for function " + *it + "\n";
            }
            throw std::runtime_error("Missing Function Definitions!");
        }
        
        return functions;
    }
    
    
    
private:
    std::string m_filename;
    
};

#endif /* InputFileParser_h */
