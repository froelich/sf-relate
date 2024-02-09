source test_param.sh
PID=0 /usr/local/go/bin/go test -run TestRelativeSearchProtocol -timeout 48h | tee /dev/tty > $OUT_FOLDER/Z/test.txt

