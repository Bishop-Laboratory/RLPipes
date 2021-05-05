#!/bin/bash

BINARY_HOME=$PREFIX/bin
RSEQ_HOME=$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM

# Copy source to the conda environment
mkdir -p $RSEQ_HOME
cp -R $SRC_DIR/* $RSEQ_HOME/

# Create symbolic links for RSeq launch script
mkdir -p $BINARY_HOME
chmod +x $RSEQ_HOME/bin/RSeq
ln -s $RSEQ_HOME/bin/RSeq $BINARY_HOME/
