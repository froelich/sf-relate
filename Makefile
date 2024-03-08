include test_param.sh
all: X Y Z

X:
	bash X_local.sh 
Y:
	bash Y_local.sh 
Z:
	bash Z_local.sh 
party1: X Z
party2: Y
