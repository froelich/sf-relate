source test_param.sh
PID=2 /usr/local/go/bin/go test -run TestRelativeSearchProtocol -timeout 48h | tee /dev/tty > $FOLDER/logs/Y/test.txt
python3 notebooks/step3_post_process.py -PARTY 2 -FOLDER $FOLDER
