#!/bin/bash
#SBATCH --job-name=4pti_23052024
#SBATCH --output=4pti_23052024.out    # Standard output
#SBATCH --error=opt%j.err         # Standard error (stderr)
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8GB
#SBATCH --time=168:00:00
#SBATCH --partition=ompi30,ompi28,ompi27,ompi26,ompi24 #preference CPUs partitions
#SBATCH --mail-type=begin # send email when job begins
#SBATCH --mail-type=end # send email when job ends
#SBATCH --mail-user=sahinjanz@hotmail.de



# create a scratch directory

SCRATCHDIR=$TMPDIR

export JOBWORKDIR=/work/cingoez/cgui

export JOBDIR=`pwd`


#test if a scratch directory has been created

if [ "$SCRATCHDIR" == "" ] ;

then echo "Could not create a scratch directory. SCRATCHDIR is empty!"

     exit 1

fi


# things that are needed for parallel runs

ulimit -s unlimited


# some input data that should be reported in the output file of the queueing system

 echo "xxxxxxxxxxxxx DvS debug xxxxxxxxxxxxx"

 echo "Path is:       $PATH"

 echo "LdPath is:     $LD_LIBRARY_PATH"

 echo "tempDIR is:    $TMPDIR"

 echo "JobworkDIR is: /work/cingoez/cgui"

 echo "JobDIR is:     $JOBDIR"

 echo "ScratchDIR is: $SCRATCHDIR"

 echo "Master:        $HOSTNAME"

 echo "Slots:         $NSLOTS"

 echo "xxxxxxxxxxxxx DvS debug xxxxxxxxxxxxx"


# create copylist to copy all files except the queuing outputfile SHNAME.out to SCRATCHDIR

COPYLIST=$(ls | grep -v 4pti_23052024.out | awk '{printf $1" "}')

# echo -e "\n $COPYLIST \n"


# create scratch directory and copy everything into it

chmod go+rwx ${SCRATCHDIR}

cd ${SCRATCHDIR}


# copy all files except the queuing outputfile SHNAME.out to SCRATCHDIR

echo -e "\n copy all except 4pti_23052024.out file to ${SCRATCHDIR}! \n"

for i in $COPYLIST

do

       cp -p -r -v /work/cingoez/cgui/${i} ${SCRATCHDIR}/

done


# copy everything back and clean up scratch even if there is an unexpected termination

cleanup_copy (){

                echo -e "\n starting cleanup procedure!\n"

                cd ${SCRATCHDIR}


                if [ -d /work/cingoez/cgui ] ;

                then

                         if [ "$SCRATCHDIR" == "" ] ;

                         then echo "SCRATCHDIR is empty!"

                              exit 1

                         fi

                         cp -p -v -r ${SCRATCHDIR}/* /work/cingoez/cgui/

                else

                         echo "$(date) :"

                         echo "JOBWORKDIR /work/cingoez/cgui does not exist or has been renamed!"

                         echo "copying files to /work/cingoez/4pti_23052024/"


                         if [ $SCRATCHDIR == "" ] ;

                         then echo "SCRATCHDIR is empty!"

                              exit 1

                         fi


                         cp -p -v -r ${SCRATCHDIR}/* /work/cingoez/4pti_23052024/

                fi


                rm -f /work/cingoez/cgui/tempdir

}


# catch TERM/QUIT signal to copy output back after unexpected termination (include "#$ -notify" in head)

trap cleanup_copy USR1 USR2 0 KILL


# run the CHARMM job,


/projects/biomodeling/charmm-lite/c46b1/ser/bin/charmm < step5_production.inp > step5_production.logfile


# copy everything after normal termination

cleanup_copy






