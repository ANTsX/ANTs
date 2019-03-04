#!/bin/bash

VERSION="0.0.1"

function Usage {
    cat <<USAGE

Usage:

`basename $0` [-h ] [-r ] [-j nb_jobs ] command arg_list

Optional arguments:

     -h:  Shows help

     -r:  Replace asterix * in the command string with argument

     -j:  Number of cpu cores to use (default 2)

Examples:

    `basename $0` somecommand arg1 arg2 arg3

    `basename $0` -j 3 \"somecommand -r -p\" arg1 arg2 arg3

    `basename $0` -j 6 -r \"convert -scale 50% * small/small_*\" *.jpg"

In case you terminate this script prematurely by pressing CTRL + C, run
${here}/killme.sh to terminate any remaining processes.

USAGE
    exit 0
}

function Help {
    cat <<HELP

This is a simple wrapper for running processes in parallel. Tested both on
Mac (Darwin) and Linux (CentOS 5).

Usage:

`basename $0` [-h ] [-r ] [-j nb_jobs ] command arg_list

Optional arguments:

     -h:  Shows this help

     -r:  Replace asterix * in the command string with argument

     -j:  Number of cpu cores to use (default 2)

Examples:

    `basename $0` somecommand arg1 arg2 arg3

    `basename $0` -j 3 \"somecommand -r -p\" arg1 arg2 arg3

    `basename $0` -j 6 -r \"convert -scale 50% * small/small_*\" *.jpg"

In case you terminate this script prematurely by pressing CTRL + C, run
${here}/killme.sh to terminate any remaining processes.

--------------------------------------------------------------------------------------
Original script by Kawakamasu:
http://pebblesinthesand.wordpress.com/category/parallel-computing/

Script adapted by:
Brian Avants, Penn Image Computing And Science Laboratory
N.M. van Strien, http://www.mri-tutorial.com | NTNU MR-Center
--------------------------------------------------------------------------------------

HELP
    exit 0
}

function queue {
    QUEUE="$QUEUE $1"
    NUM=$(($NUM+1))
}

function regeneratequeuelinux {
    OLDREQUEUE=$QUEUE
    QUEUE=""
    NUM=0
    for PID in $OLDREQUEUE
    do
        if [ -d /proc/$PID ] ; then
            QUEUE="$QUEUE $PID"
            NUM=$(($NUM+1))
        fi
    done
}

function checkqueuelinux {
    OLDCHQUEUE=$QUEUE
    for PID in $OLDCHQUEUE
    do
        if [ ! -d /proc/$PID ] ; then
            regeneratequeuelinux # at least one PID has finished
            break
        fi
    done
}

function regeneratequeuemac {
    OLDREQUEUE=$QUEUE
    QUEUE=""
    NUM=0
    for PID in $OLDREQUEUE
    do
	whm=` whoami `
	num=` ps U $whm | grep -i "${PID} "  | wc -l `
        if [ $num = 2 ]    ;  then
            QUEUE="$QUEUE $PID"
            NUM=$(($NUM+1))
	fi

    done
}

function checkqueuemac {
    OLDCHQUEUE=$QUEUE
    for PID in $OLDCHQUEUE
    do
	whm=` whoami `
	num=` ps U $whm | grep -i "${PID} "  | wc -l `
        if [ $num = 1 ]  ; then
            regeneratequeuemac # at least one PID has finished
            break
	fi
    done
}

here=`pwd`
NUM=0
QUEUE=""
MAX_NPROC=2 # default
REPLACE_CMD=0 # no replacement by default

# parse command line
if [ $# -eq 0 ]; then #  must be at least one arg
    Usage >&2
fi

while getopts j:rh OPT; do # "j:" waits for an argument "h" doesnt
    case $OPT in
    h)  Help ;;
    j)  MAX_NPROC=$OPTARG ;;
    r)  REPLACE_CMD=1 ;;
    \?) Usage >&2 ;;
    esac
done

# Main program
echo Using max $MAX_NPROC parallel threads

if [ $MAX_NPROC -eq 1 ] ; then
echo " Dont use pexec to run 1 process at a time. "
echo " In this case, just run in series. "
exit
fi

shift `expr $OPTIND - 1` # shift input args, ignore processed args
COMMAND=$1
shift

# keep list of started processes
echo "#!/bin/bash" >> ${here}/killme.sh
chmod +x ${here}/killme.sh

for INS in $* # for the rest of the arguments
do
    # DEFINE COMMAND
    if [ $REPLACE_CMD -eq 1 ]; then
        CMD=${COMMAND//"*"/$INS}
    else
        CMD="$COMMAND $INS" #append args
    fi

    echo "Running $CMD"
    eval "$CMD &"
    # DEFINE COMMAND END

    PID=$!
    echo "kill $PID" >> ${here}/killme.sh
    queue $PID

    osmac=0
    osmac="` uname -a | grep Darwin  `"
    oslin=0
    oslin="`uname -a | grep Linux`"

    if [ ${#osmac} -ne 0 ]
    then
	while [ $NUM -ge $MAX_NPROC ]; do
            checkqueuemac
            sleep 0.5
	done
    elif [ ${#oslin} -ne 0 ]
    then
	while [ $NUM -ge $MAX_NPROC ]; do
            checkqueuelinux
            sleep 0.5
	done
    fi


done

wait # wait for all processes to finish before exit

rm ${here}/killme.sh

exit 0
