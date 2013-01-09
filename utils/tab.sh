#!/bin/sh

# Copyright (C) 2012 Anthony Bradford.
# Check source code for TABs

grep -P '	' ../src/*.h ../src/*.cpp
