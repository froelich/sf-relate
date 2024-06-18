source test_param.sh
mkdir -p $FOLDER/logs/Z
mkdir -p $FOLDER/out/raw
PID=0 /usr/bin/time -o $FOLDER/logs/Z/time.txt -v ./goParty 1 > >(tee $FOLDER/logs/Z/stdout.txt ) 2> >(tee $FOLDER/logs/Z/stderr.txt >&2 )
python3 notebooks/step3_post_process.py -PARTY 1 -FOLDER $FOLDER
