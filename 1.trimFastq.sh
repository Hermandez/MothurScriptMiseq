#!/bin/bash

alias trimmomatic='java -jar /var/scratch/ytchiechoua/analysis/trimmomatic-0.39.jar'

function usage() {
	printf "Usage: %s [ options ]\n" $(basename $0);
	echo -e """
		Trim fastq files

		Options:
		-i,--idlist    <str>    :File containing SRA run accessions. [default: NULL]
		-p,--fqpath    <str>    :Path to fastq files. 
					:NB: Enter '.' for current directory 
					:[default: `pwd`/]
		-l,--leadx     <int>    :Number of bases ot crop from the 5' end [default: 0]
		-t,--trailx    <int>    :Number of bases to crop from the 3' end [default: 0]
		-T,--threads   <int>    :Optional [default = 1]
		-h,--help               :Print this help message
	"""
}

if [ $# -lt 1 ]; then
    usage; 1>&2;
    exit 1;
fi

if [ $? != 0 ]; then
   echo "ERROR: Exiting..." 1>&2;
   exit 1;
fi

prog=`getopt -o "hi:p:l:t:T:" --long "help,idlist:,fqpath:,leadx:,trailx:,threads:" -- "$@"`

id=NULL
dname="`pwd`/"
leadx=0
trailx=0
t=1

eval set -- "$prog"

while true; do
    case "$1" in
      -i|--idlist) id="$2"; 
         if [[ "$2" == -* ]]; then 
            echo "ERROR: -i,--idlist must not begin with a '-'"; 1>&2; 
            exit 1; 
         fi; 
         shift 2 
         ;;
      -p|--fqpath) dname="$2";
         if [[ "$2" == -* ]]; then
            echo "ERROR: -p,--fqpath must not begin with a '-' "; 1>&2;
            exit 1;
         fi;
         shift 2
         ;;
      -l|--leadx) leadx="$2";
         if [[ "$2" == -* ]]; then
            echo "ERROR: -l,--leadx must not begin with a '-'"; 1>&2;
            exit 1;
         fi;
         shift 2
         ;;
      -t|--trailx) trailx="$2";
         if [[ "$2" == -* ]]; then
            echo "ERROR: -t,--trailx must not begin with a '-'"; 1>&2;
            exit 1;
         fi;
         shift 2
         ;;
      -T|--threads) t="$2";
         if [[ "$2" == -* ]]; then
            echo "ERROR: -T,--threads must not begin with a '-'"; 1>&2;
            exit 1;
         fi;
         shift 2
         ;;
      -h|--help) shift; usage; 1>&2; exit 1 ;;
      --) shift; break ;;
       *) shift; usage; 1>&2; exit 1 ;;
    esac
    continue
done

if [[ $id == NULL ]]; then
    echo "ERROR: -i,--idlist not provided! Terminating..."
    sleep 1;
    1>&2;
    exit 1;
else
echo -e """=========================
Option;		argument
idlist:		$id
fqpath:		$dname
leadx:		$leadx
trailx:		$trailx
threads:	$t
========================="""

which trimmomatic
if [[ $? != 0 ]]; then
   echo "Could not locate trimmomatic"
   echo "Please install by one of the following means"
   echo "'sudo apt-get install trimmomatic'"
   echo "'conda install trimmomatic -c bioconda'"
   echo "Add bioconda channel if it doesn't already exist: 'conda config --add chennels bioconda'"
   exit $?
else
   #id=$1
   #dname=$2 # strip trailing forward slash
   #leadx=$3
   #trailx=$4
   #t=$5
   mkdir -p ${dname}/paired
   mkdir -p ${dname}/unpaired
   pdr="${dname}/paired/"
   udr="${dname}/unpaired/"
   n="$((50/$t))"
   while read -r line; do
       trimmomatic PE -phred33 $line ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:$leadx TRAILING:$trailx SLIDINGWINDOW:4:15 MINLEN:36 -threads $t
   done < $id
   mv ${dname}/*fp.fq.gz ${dname}/*rp.fq.gz $pdr
   mv ${dname}/*fu.fq.gz ${dname}/*ru.fq.gz $udr
fi
#else
#        echo """
#	Usage: ./trimFastq.sh <idlist> <fpath> <leading> <trailing> <threads>
#	 idlist: File containing SRA run accessions. (NB: Must be in same path as fastq files)
#	  fpath: Path to fastq files. (NB: Must end with a forward '/' slash)
#	threads: Optional [default = 1]
#	leading: Number of bases ot crop from the 5' end
#       trailing: Number of bases to crop from the 3' end
#	e.g.: 1_S1_L001_R1_001.fastq.gz 1_S1_L001_R2_001.fastq.gz
#	
#    """
#fi
fi
