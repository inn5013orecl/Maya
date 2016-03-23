//
//  doHelix.cpp
//  curveTest
//
//  Created by Gary Tse on 2/22/16.
//  Copyright © 2016 Gary Tse. All rights reserved.
//

#ifndef doHelix_cpp
#define doHelix_cpp

#include "doHelix.h"

void* DoHelix::creator() { return new DoHelix; }

//DeclareSimpleCommand( helix, "Gary Tse", "3.0");  //the command, the author, the min version #
MStatus DoHelix::doIt( const MArgList& args )
{
    MStatus stat;
    const unsigned	deg 	= 3;                // Curve Degree
    const unsigned	ncvs 	= 20;               // Number of CVs
    const unsigned	spans 	= ncvs - deg;       // Number of spans
    const unsigned	nknots	= spans+2*deg-1;    // Number of knots
    double	radius			= 4.0;              // Helix radius
    double	pitch 			= 0.5;              // Helix pitch
    unsigned	i;
    // Parse the arguments.
    for ( i = 0; i < args.length(); i++ )
        if ( MString( "-p" ) == args.asString( i, &stat )
            && MS::kSuccess == stat)
        {
            double tmp = args.asDouble( ++i, &stat );
            if ( MS::kSuccess == stat )
                pitch = tmp;
        }
        else if ( MString( "-r" ) == args.asString( i, &stat )
                 && MS::kSuccess == stat)
        {
            double tmp = args.asDouble( ++i, &stat );
            if ( MS::kSuccess == stat )
                radius = tmp;
        }

    // Set up CVs and knots for the helix
    MPointArray	 controlVertices;
    MDoubleArray knotSequences;

    for (i = 0; i < ncvs; i++)
        controlVertices.append( MPoint(radius * cos( (double)i ),
                                       pitch * (double)i, radius * sin( (double)i ) ) );
    for (i = 0; i < nknots; i++)
        knotSequences.append( (double)i );
    
    // Create the curve
    MFnNurbsCurve curveFn;
    MObject curve = curveFn.create(controlVertices,
                                   knotSequences, deg,
                                   MFnNurbsCurve::kOpen,
                                   false, false,
                                   MObject::kNullObj,
                                   &stat );
    if ( MS::kSuccess != stat )
        cout << "Error creating curve.\n";
    return stat;
}

DoHelix::DoHelix() {}                           // Constructor does nothing
DoHelix::~DoHelix() {}                          // Destructor does nothing

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

#endif