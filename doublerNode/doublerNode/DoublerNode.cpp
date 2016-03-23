//
//  DoublerNode.cpp
//  doublerNode
//
//  Created by Gary Tse on 3/7/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#include "DoublerNode.h"

MTypeId DoublerNode::id(0x00000001);
MObject DoublerNode::output;
MObject DoublerNode::input;

void* DoublerNode::creator() { return new DoublerNode;}

MStatus DoublerNode::compute(const MPlug &plug, MDataBlock &data) {
    if (plug != output) return MS::kUnknownParameter;
    
    //float inputValue = data.inputValue(input).asFloat();
    MDataHandle inputData = data.inputValue(input);
    
    float result = inputData.asFloat();
    
    result *= 2.0f;
    
    MDataHandle outputHandle = data.outputValue(output);
    outputHandle.setFloat(result);
    data.setClean(plug);
    
    return MS::kSuccess;
}

MStatus DoublerNode::initialize() {
    MFnNumericAttribute nAttr;
    
    output = nAttr.create("output", "out", MFnNumericData::kFloat, 0.0);
    nAttr.setWritable(false);
    nAttr.setStorable(false);
    addAttribute(output);
    
    input = nAttr.create("input", "in", MFnNumericData::kFloat, 0.0);
    nAttr.setKeyable(true);
    addAttribute(input);
    attributeAffects(input, output);
    
    return MS::kSuccess;
}

DoublerNode::DoublerNode() {};
DoublerNode::~DoublerNode() {};


