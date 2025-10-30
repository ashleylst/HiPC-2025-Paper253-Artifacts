#!/bin/bash

#TODO: add 16^4 testing here.

./main.sh --purge=1
./main.sh --build=1 --compile=1

do_test(){
	./main.sh --infile="testing_input.ini" --run=1 &>> output.txt
}

check_output() {
	if ! grep "\[  PASSED  \]" output.txt
	then 
		printf "failed test with parameters: %d %d %d %d %d %s %d\n" $mpi $num_l $s $oe $si $cl $kcycle
		cp output.txt Failed-$mpi-$num_l-$s-$oe-$si-$cl-$kcycle.log
		cp testing_input.ini Failed-$mpi-$num_l-$s-$oe-$si-$cl-$kcycle.ini
	fi
	rm output.txt
}

turn_on_MPI () {
	sed "s/d0 local lattice: 8 8 8 8/d0 local lattice: 4 4 8 8/" -i testing_input.ini
}

set_num_level() {
	sed "s/number of levels: 3/number of levels: $1/" -i testing_input.ini
}

set_smoother() {
	sed "s/method: 4/method: $1/" -i testing_input.ini
}

set_odd_even() {
	sed "s/odd even preconditioning: 0/odd even preconditioning: $1/" -i testing_input.ini
}

set_setup_iters(){
	sed "s/d0 initial setup iter: 1/d0 initial setup iter: $1/" -i testing_input.ini
}

set_clover(){
	sed "s/csw: 0/csw: $1/" -i testing_input.ini
}

set_kcycle(){
	sed "s/kcycle: 0/kcycle: $1/" -i testing_input.ini
}



run_with_setting(){
	cp sample.ini testing_input.ini
	if [[ "$1" == "1" ]]; then
		echo "reached correct if-cond."
		turn_on_MPI
	fi
	set_num_level $2
	set_smoother $3
	set_odd_even $4
	set_setup_iters $5
	set_clover $6
	set_kcycle $7

	do_test
}


#comment on input: variable must be in the following order:
# MPI on/off, number of levels, smoothers, odd_even on/off, ....


#comment on MPI on/off
# 0 = off
# 1 = on

#comment number of levels
# just put the number, eg 3

#comment on smoothers
# "-1" - pure CGN (no AMG)              #NOT SUPPORTED
#    0 - pure GMRES (no AMG)            #SUPPORTED
#    1 - FGMRES + additive Schwarz      #SUPPORTED
#    2 - FGMRES + red-black Schwarz     #SUPPORTED
#    3 - FGMRES + 16 color Schwarz      #NOT SUPPORTED
#    4 - FGMRES + GMRES                 #SUPPORTED
#    5 - FGMRES + biCGstab (no AMG)	#NOT SUPPORTED

#comment on odd_even
# 0 = off 
# 1 = on

#comment on setup iters
# just put the number, eg. 2

#comment on clover
# just put the number, eg. "1.0"

#comment on kcycle
# just put the number, eg. 1

cases_t=648
current_c=0

for mpi in 0 1
do
	for num_l in 0 1 2
	do
		for s in 0 1 2
		do 
			for oe in 0 1
			do
				for si in 0 1 2
				do 
					for cl in "0.0" "1.0"
					do
						for kcycle in 0 1
						do
							current_c=$((current_c+1))
							printf "running case %d of %d\n" $current_c $cases_t
							printf "using settings: %d %d %d %d %d %s %d\n" $mpi $num_l $s $oe $si $cl $kcycle
							run_with_setting $mpi $num_l $s $oe $si $cl $kcycle
							check_output
						done
					done
				done
			done
		done
	done
done

echo "cases in total:"
echo "$cases"
