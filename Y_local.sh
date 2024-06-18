source test_param.sh
mkdir -p $FOLDER/logs/Y
mkdir -p $FOLDER/out/raw
PID=2 /usr/bin/time -o $FOLDER/logs/Y/time.txt -v ./goParty 1 > >(tee $FOLDER/logs/Y/stdout.txt ) 2> >(tee $FOLDER/logs/Y/stderr.txt >&2 )
python3 notebooks/step3_post_process.py -PARTY 2 -FOLDER $FOLDER
