import os,sys
from pathlib import Path
r=Path(__file__).resolve().parent
env={k:v for k,v in os.environ.items() if not k.startswith(('SLURM_','PMI_','PMIX_'))}
env.update(OPAL_PREFIX=str(r),OMPI_MCA_mca_base_component_path=str(r/'openmpi'),OMPI_MCA_ess='singleton',OMPI_MCA_plm='isolated',OMPI_MCA_btl='self',PMIX_MCA_mca_base_component_path=str(r/'pmix'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
loader=str(r/'lib/ld-linux-x86-64.so.2')
os.execve(loader,[loader,'--library-path',str(r/'lib'),str(r/'bin/paper2_transport')]+sys.argv[1:],env)
