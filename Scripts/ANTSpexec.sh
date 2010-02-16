# !/bin/bash
NUM=0
QUEUE=""
MAX_NPROC=2 # default
REPLACE_CMD=0 # no replacement by default
USAGE="A simple wrapper for running processes in parallel.
Usage: `basename $0` [-h] [-r] [-j nb_jobs] command arg_list
    -h      Shows this help
    -r      Replace asterix * in the command string with argument
    -j nb_jobs  Set number of simultanious jobs [2]
 Examples:
    `basename $0` somecommand arg1 arg2 arg3
    `basename $0` -j 3 \"somecommand -r -p\" arg1 arg2 arg3
    `basename $0` -j 6 -r \"convert -scale 50% * small/small_*\" *.jpg"

function queue {
    QUEUE="$QUEUE $1"
    NUM=$(($NUM+1))
}

function regeneratequeue {
    OLDREQUEUE=$QUEUE
    QUEUE=""
    NUM=0
    for PID in $OLDREQUEUE
    do
	whm=` whoami `
	num=` ps U $whm | grep -i $PID  | wc -l `
        if [ $num = 1  ] ; then
            QUEUE="$QUEUE $PID"
            NUM=$(($NUM+1))
        fi
    done
}

function checkqueue {
    OLDCHQUEUE=$QUEUE
    for PID in $OLDCHQUEUE
    do
	whm=` whoami `
	num=` ps U $whm | grep -i $PID  | wc -l `
        if [ $num = 1   ] ; then
            regeneratequeue # at least one PID has finished
            break
        fi
    done
}

# parse command line
if [ $# -eq 0 ]; then #  must be at least one arg
    echo "$USAGE" >&2
    exit 1
fi

while getopts j:rh OPT; do # "j:" waits for an argument "h" doesnt
    case $OPT in
    h)  echo "$USAGE"
        exit 0 ;;
    j)  MAX_NPROC=$OPTARG ;;
    r)  REPLACE_CMD=1 ;;
    \?) # getopts issues an error message
        echo "$USAGE" >&2
        exit 1 ;;
    esac
done

# Main program
echo Using $MAX_NPROC parallel threads
shift `expr $OPTIND - 1` # shift input args, ignore processed args
COMMAND=$1
shift

for INS in $* # for the rest of the arguments
do
    # DEFINE COMMAND
    if [ $REPLACE_CMD -eq 1 ]; then
        CMD=${COMMAND//"*"/$INS}
    else
        CMD="$COMMAND $INS" #append args
    fi
    echo "Running $CMD"

    $CMD &
    # DEFINE COMMAND END

    PID=$!
    queue $PID

    while [ $NUM -ge $MAX_NPROC ]; do
        checkqueue
        sleep 0.4
    done
done
wait # wait for all processes to finish before exit
