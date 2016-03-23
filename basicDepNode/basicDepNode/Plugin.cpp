// TM and (c) 2016 Gary Tse. All Rights Reserved.
// Reproduction in whole or in part without prior written permission of a
// duly authorized representative is prohibited.
#include <maya/MFnPlugin.h>
#include "basicDepNode.h"

MStatus
initializePlugin(MObject plugin)
{
    MStatus status;
    MFnPlugin pluginFn(plugin, "Gary Tse", "1.0.0", "Any");
    
    // Register elements here
    
    status = pluginFn.registerNode("sine", Sine::id, Sine::creator, Sine::initialize);
    
    if (!status) {
        status.perror("registerNode");
        return status;
    }
    
    return status;
}

MStatus
uninitializePlugin(MObject plugin)
{
    MStatus status;
    MFnPlugin pluginFn(plugin);
    
    // Deregister elements here
    
    status = pluginFn.deregisterNode(Sine::id);
    
    if (!status) {
        status.perror("deregisterNode");
        return status;
    }
    
    return status;
}
// TM and (c) 2016 Gary Tse. All Rights Reserved.
// Reproduction in whole or in part without prior written permission of a
// duly authorized representative is prohibited.
