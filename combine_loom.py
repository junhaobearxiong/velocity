import loompy

input_dirname = 'velocyto_output/'
output_path = 'data/E1.loom'

files = [
	'finalmerged_E1CD1col1_Aligned_UYDM8.loom',
	'finalmerged_E1CD1col2_Aligned_O5OMP.loom',
	'finalmerged_E1CD1col3_Aligned_I9VO5.loom',
	'finalmerged_E1CD1col4_Aligned_8WF7Q.loom',
	'finalmerged_E1CD1col5_Aligned_4UFUT.loom',
	'finalmerged_E1CD1col6_Aligned_71KUA.loom',
	'finalmerged_E1CD2col1_Aligned_LTERJ.loom',
	'finalmerged_E1CD2col2_Aligned_1BCK5.loom',
	'finalmerged_E1CD2col3_Aligned_RAFAO.loom',
	'finalmerged_E1CD2col4_Aligned_9NBS5.loom',
	'finalmerged_E1CD2col5_Aligned_LH0JZ.loom',  
	'finalmerged_E1CD2col6_Aligned_NHHLR.loom',
	'finalmerged_E1CD3col1_Aligned_TSK6M.loom',  
	'finalmerged_E1CD3col2_Aligned_T2DC5.loom',
	'finalmerged_E1CD3col3_Aligned_JBAUL.loom',
	'finalmerged_E1CD3col4_Aligned_0DE7D.loom',
	'finalmerged_E1CD3col5_Aligned_J4BOP.loom',
	'finalmerged_E1CD3col6_Aligned_PFFAE.loom',
]

files = [input_dirname + f for f in files]


loompy.combine(files, output_path, key='Accession')