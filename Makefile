include test_param.sh
all: X Y Z

X:
	mkdir out/$t/X -p
	bash X_local.sh 
Y:
	mkdir out/$t/Y -p
	bash Y_local.sh 
Z:
	mkdir out/$t/Z -p 
	bash Z_local.sh 
party1: X Z
party2: Y

