source test_param.sh
PID=1 /usr/local/go/bin/go test -run TestRelativeSearchProtocol -timeout 48h | tee /dev/tty > $OUT_FOLDER/X/test.txt
