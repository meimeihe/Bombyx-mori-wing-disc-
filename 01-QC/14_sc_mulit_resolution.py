import scanpy as sc
import pandas as pd
import scanpy.external as sce
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter
import matplotlib.font_manager
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
font_size = 10
rc={'font.size': font_size, 'axes.labelsize': font_size, 'figure.dpi':400, 'axes.linewidth':1,
    'axes.titlesize': font_size, 'xtick.labelsize': font_size, 'ytick.labelsize': font_size} # 'figure.figsize':(11.7/1.5,8.27/1.5)



matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

plt.rcParams['axes.unicode_minus']=False # negative minus sign

centimeter = 1/2.54

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()

sc.settings.set_figure_params(dpi=200, dpi_save=400,facecolor='white')

os.chdir("/hsfscqjf1/ST_CQ/P23Z32300N0001/hemingmin/10.sc_merged1/results/2_scanpy_base_timepoint/")
data_path = "/hsfscqjf1/ST_CQ/P23Z32300N0001/hemingmin/10.sc_merged1/results/"
adata_filter_v2 = sc.read_h5ad(data_path+"seurat_filtered.h5ad")
adata_filter_v2.raw = adata_filter_v2

sc.pp.normalize_total(adata_filter_v2, target_sum=1e4,inplace=True)
sc.pp.log1p(adata_filter_v2)

# cell_cycle_genes = pd.read_csv("cell cycle gene list.csv",sep=',')

# s_genes=cell_cycle_genes['G1_S_gene']
# s_genes.dropna(inplace=True)
# s_genes = s_genes.to_list()

# G2M_genes = cell_cycle_genes['G2M_gene'].to_list()

# cell_cycle_genes = s_genes + G2M_genes
# cell_cycle_genes = [x for x in cell_cycle_genes if x in adata_filter_v2.var_names]

# sc.tl.score_genes_cell_cycle(adata_filter_v2, s_genes=s_genes, g2m_genes=G2M_genes)

# sc.pp.regress_out(adata_filter_v2, ['S_score', 'G2M_score'])

sc.pp.highly_variable_genes(adata_filter_v2,n_top_genes=2500)
adata_filter_v2 = adata_filter_v2[:, adata_filter_v2.var.highly_variable]

sc.pp.scale(adata_filter_v2)

sc.tl.pca(adata_filter_v2, n_comps=80, random_state=0, use_highly_variable=False)

sce.pp.harmony_integrate(adata_filter_v2,key='sample')
# sce.pp.harmony_integrate(adata_filter_v2,key='timepoint')

sc.pp.neighbors(adata_filter_v2,use_rep='X_pca_harmony',n_pcs=30,n_neighbors=30,random_state=0)
# sc.pp.neighbors(adata_filter,n_pcs=40,random_state=0)

sc.tl.tsne(adata_filter_v2,random_state=0)
sc.tl.umap(adata_filter_v2,random_state=0)

res=[0.4,0.5,0.6,0.8,1]
for i in res:
    sc.tl.leiden(adata_filter_v2, key_added="leiden_res"+str(i), resolution=i)


for j in res:
    sc.pl.umap(
        adata_filter_v2, 
        color="leiden_res"+str(j),
        size=1,
       # legend_loc='on data',
        palette=sc.pl.palettes.default_28,
        wspace=0.3,
        ncols=1,
        title="leiden_res"+str(j),
        save = "leiden_res"+str(j)+"leiden_umap.pdf",
        frameon=False,
        show = False)



adata_filter_v2 = adata_filter_v2.raw.to_adata()
adata_filter_v2.write("cluster_raw.h5ad",compression='gzip')

