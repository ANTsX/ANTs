#!/usr/bin/perl -w
use File::Basename;

    if ($#ARGV >= 0) { $who = join(' ', $ARGV[0]); }
    if ($#ARGV >= 1) { $template = join(' ', $ARGV[1]); }
    if ($#ARGV >= 2) { $params = join(' ', $ARGV[2]); } else { $params=0; }
    if ($#ARGV >= 3) { $outnaming = join(' ', $ARGV[3]); } else { $outnaming=0; }
    if ($#ARGV >= 4) { $metric= join(' ', $ARGV[4]); } else { $metric="PR[".$template; }
    if ($#ARGV >= 5) { $metricpar= join(' ', $ARGV[5]); } else { $metricpar=",1,4]"; }

    opendir(DIR, $who);
    @filelist = glob("$who");
#    @filelist2 = glob("$who2");
    $ct=100;
print " Expecting Normalization params ".$params." \n \n \n doing \n $who \n ";
$dir = `pwd`;
$dir=substr($dir,0,length($dir)-1)."/";
$pathpre="";
$qsname="superfresh.sh";
    foreach $name (@filelist)
    {
      $pre=$outnaming.$name;
      $pre=~s/\.[^.]*$//;
      print $pre." \n  ";
      $prog=$params;
      $exe=$prog."  -m ".$metric.",".$name.",".$metricpar." -o ".$pre."  \n ";
        $exe2 = "/mnt/pkg/sge-root/bin/lx24-x86/qsub -q mac ".$qsname."  ".$dir." ".$exe." ";
        print $exe."\n";
     system $exe2;
# USE qstat TO VERIFY THE JOB IS SUBMITTED
	$ct=$ct+1;
    }

print " DONE WITH LOOP \n ";
$user=`whoami`;
chomp($user);
print " you are ".$user." = surf-king ";
$atcfn=0;

if ($atcfn)  { $waitcmd="vq -s | grep ".$user; }
else { $waitcmd=" qstat -u ".$user." | grep ".substr($qsname,0,4); }
      $continuewaiting=1;
if ($continuewaiting) { system "sleep 40";  }
      # Now Wait for the jobs to complete
      while($continuewaiting)
      {
	system $waitcmd;
	$output = `$waitcmd`;

	if ( $output eq "" )
	{
	    $continuewaiting=0;
	    print "Jobs all finished!\n";
	}
	else { print " wait ";system "sleep 30"; }
       }


    closedir(DIR);
