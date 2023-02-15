# /bin/bash !


read -p  "how many experiments do you intend to run: " nexp
read -p  "experiment series name (short), e.g. kspd : " sern


for i in `seq 1 $nexp`; do 
	fn1=$sern$i
	if [[ $nexp -ge 10 &&  $i -lt 10  ]] ; then
		fn1=$sern"0"$i
	fi
	mkdir ../$fn1
	cp DAC ../$fn1
	mkdir ../$fn1/RUN5
	mkdir ../$fn1/RESTART
	mkdir ../$fn1/ASCII
done


