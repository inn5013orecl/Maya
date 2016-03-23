// TM and (c) 2016 Gary Tse. All Rights Reserved.
// Reproduction in whole or in part without prior written permission of a
// duly authorized representative is prohibited.

#ifndef Plugin_cpp
#define Plugin_cpp

#include "doHelix.h"
#include <maya/MFnPlugin.h>

MStatus
initializePlugin(MObject plugin)
{
    //MStatus status;
    MFnPlugin pluginFn(plugin, "Gary Tse", "1.0.0", "Any");     // object, author, plugin version, Maya version required
    
    // Register elements here
    MStatus status = pluginFn.registerCommand("helix", DoHelix::creator);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    return status;
}

MStatus
uninitializePlugin(MObject plugin)
{
    //MStatus status;
    MFnPlugin pluginFn(plugin);
    
    // Deregister elements here
    MStatus status = pluginFn.deregisterCommand("helix");
    CHECK_MSTATUS_AND_RETURN_IT(status);
    return status;
}
// TM and (c) 2016 Gary Tse. All Rights Reserved.
// Reproduction in whole or in part without prior written permission of a
// duly authorized representative is prohibited.

#endif