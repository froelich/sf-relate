source test_param.sh
mkdir -p $FOLDER/logs/X
mkdir -p $FOLDER/out/raw
PID=1 /usr/bin/time -o $FOLDER/logs/X/time.txt -v ./goParty 1 > >(tee $FOLDER/logs/X/stdout.txt ) 2> >(tee $FOLDER/logs/X/stderr.txt >&2 )
python3 notebooks/step3_post_process.py -PARTY 1 -FOLDER $FOLDER
