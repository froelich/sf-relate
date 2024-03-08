source test_param.sh
PID=1 /usr/local/go/bin/go test -run TestRelativeSearchProtocol -timeout 48h | tee /dev/tty > $FOLDER/logs/X/test.txt
python3 notebooks/step3_post_process.py -PARTY 1 -FOLDER $FOLDER
