include test_param.sh
all: X Y Z

full_pipeline:
	bash 0_prepare_1KG.sh
	bash 1_hashing.sh
	bash 2_sketch.sh
	bash 3_run_MHE.sh
	bash 4_verify_output.sh

X:
	mkdir out/$t/X -p
	bash X_local.sh 
Y:
	mkdir out/$t/Y -p
	bash Y_local.sh 
Z:
	mkdir out/$t/Z -p 
	bash Z_local.sh 
dev: Y Z
devi: X Z
devv: X Y

