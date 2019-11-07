#!/bin/bash

make EXEC=swe_default.exe CONFIG_DIR=gnu CFDL=0 SFDL=0 SPLIT_DEV=0
make EXEC=swe_sfdl.exe CONFIG_DIR=gnu CFDL=0 SFDL=1 SPLIT_DEV=0
make EXEC=swe_cfdl.exe CONFIG_DIR=gnu CFDL=1 SFDL=0 SPLIT_DEV=0
make EXEC=swe_cfdl_sfdl.exe CONFIG_DIR=gnu CFDL=1 SFDL=1 SPLIT_DEV=0
