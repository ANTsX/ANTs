#!/bin/bash

## the intent of this file is to clean the files, make them KWSTYLE compliant,
## and then make sure that subsequent runs of this script do not modify `
## The goal of this script is to develop a consistent system of auto-formating
## and format testing such that if uncrustify is run, then the code will pass
## the KWStyle format checker

## http://www.itk.org/Wiki/ITKv4_StyleChangeProposal

UNCRUSTIFYBIN=/opt/uncrustify/bin/uncrustify
if [ ! -f ${UNCRUSTIFYBIN} ];then
   UNCRUSTIFYBIN=/usr/local/bin/uncrustify
fi

CONFIG=$(dirname $0)/uncrustify_itk_aggressive.cfg
if [ ! -f "${CONFIG}" ]; then
  CONFIG=/raid0/homes/johnsonhj/Dashboard/src/ITK/Utilities/Maintenance/uncrustify_itk_aggressive.cfg
fi
if [ ! -f "${CONFIG}" ]; then
  CONFIG=/Users/johnsonhj/Dashboards/src/ITK/Utilities/Maintenance/uncrustify_itk_aggressive.cfg
fi
if [ ! -f "${CONFIG}" ]; then
  echo "failed to find config file"
  exit -1
fi


for file in $@; do
    echo "======= ${file}"
    if [ ! -f ${KWERROS}_skip_uncrustify ]; then
        ### A possible bug in uncrustify breaks complicated macros by adding extra spaces surrounding "##" and "#"
        ### This is an interaction betweeen the macro directives and indenting of includes and defines.
        #echo ${UNCRUSTIFYBIN} -c ${CONFIG} -l CPP -f ${file} | sed 's/  *##  */##/g' | sed 's/#  */#/g'
        ${UNCRUSTIFYBIN} -c ${CONFIG} -l CPP -f ${file} | sed 's/  *##  */##/g' | sed 's/#  */#/g' |sed 's///g' > ${file}_uncrustify
        mv ${file}_uncrustify ${file}
    fi
done
