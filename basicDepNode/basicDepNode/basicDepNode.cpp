//
//  basicDepNode.cpp
//  basicDepNode
//
//  Created by Gary Tse on 3/7/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#include "basicDepNode.h"

MTypeId             Sine::id(0x00000001);
MObject             Sine::input;
MObject             Sine::output;

void* Sine::creator() {
    return  new     Sine;
}

MStatus Sine::initialize() {
    
    MFnNumericAttribute nAttr;
    
    output = nAttr.create("output", "out", MFnNumericData::kFloat, 0.0);
    nAttr.setWritable(false);
    nAttr.setStorable(false);
    
    input = nAttr.create("input", "in", MFnNumericData::kFloat, 0.0);
    nAttr.setStorable(true);
    
    addAttribute(input);
    attributeAffects(input, output);
    addAttribute(output);
    
    return MS::kSuccess;
}

Sine::Sine() {};
Sine::~Sine() {};

MStatus Sine::compute(const MPlug &plug, MDataBlock &data) {
    MStatus stat;
    
    if (plug == output) {
        MDataHandle     inputData       = data.inputValue(input, &stat);
    
        if (MS::kSuccess != stat) {
            MGlobal::displayInfo("Error getting data.");
        }
        else {
            float       result          = sin(inputData.asFloat());     //the as*() type must match what it was created as
            MDataHandle outputHandle    = data.outputValue(output);
            
            outputHandle.set(result);
            data.setClean(plug);
        }
    }
    return stat;
}