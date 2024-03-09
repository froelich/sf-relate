source test_param.sh
PID=1 ./goParty | tee /dev/tty > $FOLDER/logs/X/test.txt
python3 notebooks/step3_post_process.py -PARTY 1 -FOLDER $FOLDER
