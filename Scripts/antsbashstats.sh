#!/bin/bash

# sttdev.sh: Standard Deviation

# Original version obtained from: http://tldp.org/LDP/abs/html/contributed-scripts.html#STDDEV

# 2009-03-09: Panagiotis Kritikakos
# Function stat_median() added for calculating the median value
# of the data set

# ------------------------------------------------------------
#  The Standard Deviation indicates how consistent a set of data is.
#  It shows to what extent the individual data points deviate from the
#+ arithmetic mean, i.e., how much they "bounce around" (or cluster).
#  It is essentially the average deviation-distance of the
#+ data points from the mean.

# =========================================================== #
#    To calculate the Standard Deviation:
#
# 1  Find the arithmetic mean (average) of all the data points.
# 2  Subtract each data point from the arithmetic mean,
#    and square that difference.
# 3  Add all of the individual difference-squares in # 2.
# 4  Divide the sum in # 3 by the number of data points.
#    This is known as the "variance."
# 5  The square root of # 4 gives the Standard Deviation.
# =========================================================== #

count=0         # Number of data points; global.
SC=9            # Scale to be used by bc. Nine decimal places.
E_DATAFILE=90   # Data file error.

# ----------------- Set data file ---------------------
if [ ! -z "$1" ]  # Specify filename as cmd-line arg?
then
  datafile="$1" #  ASCII text file,
else            #+ one (numerical) data point per line!
  datafile=sample.dat
fi              #  See example data file, below.

if [ ! -e "$datafile" ]
then
  echo "\""$datafile"\" does not exist!"
  echo " Usage : "
  echo " sh $0 textdatafile.txt "
  exit $E_DATAFILE
fi
# -----------------------------------------------------
arith_mean ()
{
  local rt=0         # Running total.
  local am=0         # Arithmetic mean.
  local ct=0         # Number of data points.

  while read value   # Read one data point at a time.
  do
      echo $value
    rt=$(echo "scale=$SC; $rt + $value" | bc)
    (( ct++ ))
  done

  am=$(echo "scale=$SC; $rt / $ct" | bc)

  echo $am; return $ct   # This function "returns" TWO values!
  #  Caution: This little trick will not work if $ct > 255!
  #  To handle a larger number of data points,
  #+ simply comment out the "return $ct" above.
} <"$datafile"   # Feed in data file.

sd ()
{
  mean1=$1  # Arithmetic mean (passed to function).
  n=$2      # How many data points.
  sum2=0    # Sum of squared differences ("variance").
  avg2=0    # Average of $sum2.
  sdev=0    # Standard Deviation.

  while read value   # Read one line at a time.
  do
    diff=$(echo "scale=$SC; $mean1 - $value" | bc)
    # Difference between arith. mean and data point.
    dif2=$(echo "scale=$SC; $diff * $diff" | bc) # Squared.
    sum2=$(echo "scale=$SC; $sum2 + $dif2" | bc) # Sum of squares.
  done

    avg2=$(echo "scale=$SC; $sum2 / $n" | bc)  # Avg. of sum of squares.
    sdev=$(echo "scale=$SC; sqrt($avg2)" | bc) # Square root =
    echo $sdev                                 # Standard Deviation.

} <"$datafile"   # Rewinds data file.

stat_median() {
  NUMS=(`sort -n $1`)
  TOTALNUMS=${#NUMS[*]}
  MOD=$(($TOTALNUMS % 2))
  if [ $MOD -eq 0 ]; then
    ARRAYMIDDLE=$(echo "($TOTALNUMS / 2)-1" | bc)
    ARRAYNEXTMIDDLE=$(($ARRAYMIDDLE + 1))
    MEDIAN=$(echo "scale=$SC; ((${NUMS[$ARRAYMIDDLE]})+(${NUMS[$ARRAYNEXTMIDDLE]})) / 2" | bc)
  elif [ $MOD -eq 1 ]; then
    ARRAYMIDDLE=$(echo "($TOTALNUMS / 2)" | bc)
    MEDIAN=${NUMS[$ARRAYMIDDLE]}
  fi
  echo $MEDIAN
}

# ======================================================= #
mean=$(arith_mean); count=$?   # Two returns from function!
std_dev=$(sd $mean $count)
median=$(stat_median $1)

echo " Mean : $mean  SD: $std_dev "
#echo
#echo "Number of data points in \""$datafile"\" = $count"
#echo "Arithmetic mean (average) = $mean"
#echo "Standard Deviation = $std_dev"
#echo "Median number (middle) = $median"
#echo
# ======================================================= #

exit
