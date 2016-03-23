// TM and (c) 2016 Gary Tse. All Rights Reserved.
// Reproduction in whole or in part without prior written permission of a
// duly authorized representative is prohibited.
#include <maya/MFnPlugin.h>
#include "DoublerNode.h"

MStatus
initializePlugin(MObject plugin)
{
    MStatus status;
    MFnPlugin pluginFn(plugin, "Gary Tse", "1.0.0", "Any");
    
    // Register elements here
    
    status = pluginFn.registerNode("doublerNode", DoublerNode::id, DoublerNode::creator, DoublerNode::initialize);
    
    CHECK_MSTATUS_AND_RETURN_IT(status);
    
    return status;
}

MStatus
uninitializePlugin(MObject plugin)
{
    MStatus status;
    MFnPlugin pluginFn(plugin);
    
    // Deregister elements here
    
    status = pluginFn.deregisterNode(DoublerNode::id);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    
    return status;
}
// TM and (c) 2016 Gary Tse. All Rights Reserved.
// Reproduction in whole or in part without prior written permission of a
// duly authorized representative is prohibited.
