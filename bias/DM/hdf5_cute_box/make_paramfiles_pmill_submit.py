import numpy as N
import os
import sys


boxsize=542.16

def run_all(snapnum):
	data_fn='/cosma5/data/jch/L800/Sample/data/snapshot_%i'%snapnum
	output_fn='/gpfs/data/nmccull/pmill/xi_real_sn%03i.txt'%snapnum
	param_fn='/cosma/home/nmccull/code/CUTE-1.2/CUTE_box/paramfiles/pmill/params_sn%i.txt'%snapnum
	batch_fn='/cosma/home/nmccull/code/scripts/pmill/batch_script_cute%03i.sh'%snapnum
	print batch_fn
	print param_fn
	make_param_file(snapnum, param_fn, data_fn, output_fn)
	make_batch_file(batch_fn, param_fn)
	command='bsub < %s'%batch_fn
	status = os.system(command)


def make_param_file(snapnum, param_fn, data_fn, output_fn):
	f=open(param_fn, 'w')
	f.write('#input information\n')
	f.write('data_filename= %s\n'%data_fn)
	f.write('num_lines= all\n')
	f.write('input_format= 5\n')
	f.write('output_filename= %s\n'%output_fn)
	f.write('box_size= %f\n'%boxsize)

	f.write('#tree stuff\n')
	f.write('use_tree= 0\n')
	f.write('max_tree_order= 6\n')
	f.write('max_tree_nparts= 100\n')

	f.write('#irrelevant yet\n')
	f.write('use_pm= 0\n')
	f.write('n_grid_side= 256\n')
	f.close()

def make_batch_file(batch_fn, param_fn):



	f=open(batch_fn, 'w')
	f.write('#!/bin/tcsh\n')
	f.write('#BSUB -L /bin/tcsh\n')
	f.write('#BSUB -n 12\n')
	f.write('#BSUB -J cutexi\n')
	f.write('#BSUB -oo log/cute_%J.out\n')
	f.write('#BSUB -eo log/cute_%J.err\n')
	f.write('#BSUB -q cordelia\n')
	f.write('#BSUB "ptile[hosts=1]"\n')
	f.write('#BSUB -P durham\n')
	f.write('#BSUB -W 200:00\n')
	f.write('setenv OMP_NUM_THREADS 12\n')
	f.write('./CUTE_box %s\n'%param_fn)
	f.close()









if __name__=="__main__":

	snapnum=N.int(sys.argv[1])
	run_all(snapnum)
