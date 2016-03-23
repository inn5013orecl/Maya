//
//  helloWorld.hpp
//  mayaHelloWorld
//
//  Created by Gary Tse on 2/22/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#ifndef HELLOWORLD_H
#define HELLOWORLD_H

#include <maya/MArgList.h>
#include <maya/MObject.h>
#include <maya/MGlobal.h>
#include <maya/MPxCommand.h>

class HelloWorld : public MPxCommand {
public:
    HelloWorld() {};
    virtual MStatus doIt(const MArgList& argList);
    static void* creator();
};
#endif
