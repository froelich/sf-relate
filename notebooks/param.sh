# for 1KG datasets, using more lenient choices: smaller ell = 3, shorter splits of 3cM and longer encoding (len=160) works better.
n=1601
w=1
k=8
l=2
L=4
snpmin=160
target_len=160
samp_mode=4
len_seg_cM=4
rep=1
rid=2
LSH=Hamming
l_array=( $l )
tid=mode${samp_mode}cM${len_seg_cM}len${target_len}k${k}