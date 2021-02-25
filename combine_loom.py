import loompy

input_dirname = 'velocyto_output/'
output_path = 'data/E2.loom'

'''
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
'''

files = [
    'finalmerged_E2CD1col1_Aligned_LQF2W.loom',
    'finalmerged_E2CD1col2_Aligned_3JT6M.loom',
    'finalmerged_E2CD1col3_Aligned_WKP9S.loom',
    'finalmerged_E2CD1col4_Aligned_UP0YP.loom',
    'finalmerged_E2CD1col5_Aligned_2DWQ1.loom',
    'finalmerged_E2CD1col6_Aligned_S80Q8.loom',
    'finalmerged_E2CD2col1_Aligned_0Y90Q.loom',
    'finalmerged_E2CD2col2_Aligned_QQEP0.loom',
    'finalmerged_E2CD2col3_Aligned_5673N.loom',
    'finalmerged_E2CD2col4_Aligned_3ZYCP.loom',
    'finalmerged_E2CD2col5_Aligned_14D5V.loom',
    'finalmerged_E2CD2col6_Aligned_KA7T4.loom',
    'finalmerged_E2CD3col1_Aligned_V5IPL.loom',
    'finalmerged_E2CD3col2_Aligned_46E83.loom',
    'finalmerged_E2CD3col3_Aligned_8KXJE.loom',
    'finalmerged_E2CD3col4_Aligned_U67XA.loom',
    'finalmerged_E2CD3col5_Aligned_ZLNP3.loom',
    'finalmerged_E2CD3col6_Aligned_4R8BF.loom'.
]

files = [input_dirname + f for f in files]


loompy.combine(files, output_path, key='Accession')