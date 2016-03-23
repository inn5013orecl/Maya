#!/bin/sh
INSTALL_ROOT=$HOME/Library/Preferences/Autodesk/maya/$MAYA_VERSION
mkdir -p $INSTALL_ROOT/plug-ins
cp $BUILT_PRODUCTS_DIR/$PRODUCT_NAME.bundle $INSTALL_ROOT/plug-ins
										
