include test_param.sh
all: X Y Z

compile:
	go get relativeMatch
	go build
	go test -c -o goParty

X: compile
	bash X_local.sh 
Y: compile
	bash Y_local.sh 
Z: compile
	bash Z_local.sh 
party1: X Z
party2: Y
