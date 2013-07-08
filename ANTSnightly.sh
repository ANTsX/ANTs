#!/bin/bash
#
package=ANTS
repository=git@github.com:BRAINSia/ANTS.git
#
# when run by cron, the path variable is only /bin:/usr/bin
export PATH="/opt/cmake/bin:/usr/local/bin:/usr/sbin:$PATH"

#
# make the testing directory, based on current user name
#user=`who -m | sed -e 's/ .*$//'`
user=${LOGNAME}

ThisComputer=`hostname`


#
# the default is to use /brainsdev/kent -- which is
# appropriate on the b2dev VMs.
if [ $# = 0 ] ; then
    startdir=/scratch/kent/Testing
else
    startdir=$1
    shift
fi

#
# needed for ssl authentication for git
export GIT_SSL_NO_VERIFY=true

CXXFLAGS="${CXXFLAGS:-}"
CFLAGS="${CFLAGS:-}"
LDFLAGS="${LDFLAGS:-}"

# turn on coverage at command line
if [ $# = 0 ] ; then
    coverage=0
else
    if [ $1 = "coverage" ] ; then
	coverage=1
	shift
    fi
fi

OS=$(uname -s)
NPROCS=1

# if [ "${OS}" = "Linux" ] ; then
#     NPROCS=$(grep -c ^processor /proc/cpuinfo)
#     export CFLAGS="${CFLAGS} -fpic"
#     export CXXFLAGS="${CXXFLAGS} -fpic"
# else
#     NPROCS=$(system_profiler | awk '/Number Of Cores/{print $5}{next;}')
# fi

# create the testing directory if necessary
mkdir -p ${startdir}
if [ ! -d ${startdir} ] ; then
    echo ${startdir} cannot be created, exiting
    exit 1
fi

cd ${startdir}

mkdir -p ${startdir}/${ThisComputer}/${package}

cd ${startdir}/${ThisComputer}/${package}

top=`pwd`
echo WORKING IN $top

# check out package in a directory unique to each host -- this is unfortunately necessary
# because svn can't update a directory  checked out by a newer version of svn, so
# every host has their own copy of BRAINS3 so that it's compatible with the local svn version.
if [ -d ${package} ] ; then
    cd ${package}
    git pull
else
    git clone ${repository}
fi

if [ $? != 0 ]
then
    echo ${package} checkout failed, continuing with old version
fi


OsName=$(uname)
if [ "${OsName}" = "Darwin" ] ; then
    Compiler=clang-`clang -v 2>&1 | head -1 | awk '{print $4}'`
    Compiler=${Compiler}-`clang -v 2>&1 | tail -2 | head -1 | awk '{print $2}'`
    export CC=`which clang`
    export CXX=`which clang++`
else
    which gcc > /dev/null 2>&1
    if [ $? == 0 ] ; then
        Compiler=gcc-`gcc -dumpversion`-`gcc -dumpmachine`
    else
        Compiler=unknown
    fi
fi

echo "Compiler=${Compiler} CC=${CC} CXX=${CXX}"

for BUILD_TYPE in Debug Release
do
    BuildDir=${top}/${BUILD_TYPE}
    if [ "$BUILD_TYPE" = "Debug" -a "$coverage" = "1" ] ; then
	CXXFLAGS="${CXXFLAGS} -g -O0 -Wall -W -Wshadow -Wunused-variable \
	    -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers \
	    -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage"
	CFLAGS="${CFLAGS} -g -O0 -Wall -W -fprofile-arcs -ftest-coverage"
	LDFLAGS="${LDFLAGS} -fprofile-arcs -ftest-coverage"
    fi
    mkdir -p ${BuildDir}
    cd ${BuildDir}
    rm -f CMakeCache.txt
    # force reconfigure.
    find . -name '*-configure' | xargs rm -f
    #
    # the Build type
    cmake -DSITE:STRING=${ThisComputer} \
        -G "Unix Makefiles" \
	-DCMAKE_C_FLAGS:STRING="${CFLAGS}" \
	-DCMAKE_CXX_FLAGS:STRING="${CXXFLAGS}" \
	-DCMAKE_EXE_LINKER_FLAGS:STRING="${LDFLAGS}" \
	-DCMAKE_MODULE_LINKER_FLAGS:STRING="${LDFLAGS}" \
	-DCMAKE_SHARED_LINKER_FLAGS:STRING="${LDFLAGS}" \
	-DBUILDNAME:STRING="${OsName}-${Compiler}-${BUILD_TYPE}" \
        -DBUILD_SHARED_LIBS:BOOL=Off \
	-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE} \
        ${top}/${package}
    echo "Building in `pwd`"
    scriptname=`basename $0`
    make -j ${NPROCS}
    cd ${package}-build
    make clean
    if [ $scriptname = "nightly.sh" ] ; then
	ctest -j ${NPROCS} -D Nightly
    else
	ctest -j ${NPROCS} -D Experimental
    fi
    cd ..
done

cd ${top}
