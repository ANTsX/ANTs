#!/bin/bash
# request Bourne shell as shell for job

ctest -D ExperimentalStart
ctest -D ExperimentalUpdate
ctest -D ExperimentalConfigure
ctest -D ExperimentalBuild
ctest -D ExperimentalSubmit
ctest -D ExperimentalTest
#ctest -D ExperimentalCoverage
ctest -D ExperimentalSubmit
#ctest -D ExperimentalMemCheck
#ctest -D ExperimentalSubmit
