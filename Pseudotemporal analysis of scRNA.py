### Cellrank
import cellrank as cr
import scanpy as sc
import os
import matplotlib.pyplot as plt
os.chdir("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2D/cellrank/candidate/GPCCA/PseudotimeKernel")
sc.settings.set_figure_params(frameon=False, dpi=300)
cr.settings.verbosity = 2
import warnings
warnings.simplefilter("ignore", category=UserWarning)
file_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/cellrank/MNFHC_cellrank.h5ad"
adata = sc.read_h5ad(file_path)
sc.pp.filter_genes(adata, min_cells=5)
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
from cellrank.kernels import PseudotimeKernel
pk = PseudotimeKernel(adata, time_key="velocity_pseudotime")
pk.compute_transition_matrix()
fig = pk.plot_random_walks(
    seed=0,
    n_sims=100,
    start_ixs={"leiden_mEPC_CAF": "mEPC 08"},
    basis="umap",  # 修改 basis 参数为 "umap"
    legend_loc="right",
    dpi=300,
)
if fig is None:
    fig = plt.gcf()

fig.savefig("plot_random_walks.png", format="png", dpi=300)
plt.close(fig)
from cellrank.estimators import GPCCA
g = GPCCA(pk)
g.fit(n_states=16, cluster_key="leiden_mEPC_CAF") # 因有16种细胞亚型因此选择16个状态
fig = g.plot_macrostates(which="all")
if fig is None:
    fig = plt.gcf()
output_path = os.path.join(os.getcwd(), "macrostates.png")
fig.savefig(output_path, format="png", dpi=300)
plt.close(fig)
g.predict_initial_states(n_states=6, allow_overlap=True)
fig = g.plot_macrostates(which="initial")
if fig is None:
    fig = plt.gcf()
output_path = os.path.join(os.getcwd(), "initial_macrostates.png")
fig.savefig(output_path, format="png", dpi=300)
plt.close(fig)

print(f"终末状态图形已保存至: {output_path}")
g.predict_terminal_states(method="top_n", n_states=6, allow_overlap=True)
fig = g.plot_macrostates(which="terminal")
if fig is None:
    fig = plt.gcf()
output_path = os.path.join(os.getcwd(), "terminal_macrostates.png")
fig.savefig(output_path, format="png", dpi=300)
plt.close(fig)
g.compute_fate_probabilities()
fig = g.plot_fate_probabilities(legend_loc="right")
if fig is None:
    fig = plt.gcf()

output_path = os.path.join(os.getcwd(), "fate_probabilities.png")
fig.savefig(output_path, format="png", dpi=300)
plt.close(fig)
fig = cr.pl.circular_projection(adata, keys="leiden_mEPC_CAF", legend_loc="right")
if fig is None:
    fig = plt.gcf()
output_path = os.path.join(os.getcwd(), "circular_projection.png")
fig.savefig(output_path, format="png", dpi=300)
plt.close(fig)
terminal_states = [
    'mEPC 04_1', 'mEPC 04_2', 'mEPC 04_3', 'mEPC 03_1', 
    'CAF 02', 'CAF 06_1', 'mEPC 04_4', 'CAF 04', 
    'CAF 06_2', 'mEPC 08', 'CAF 05', 'mEPC 06', 
    'CAF 01', 'mEPC 03_2', 'mEPC 05', 'mEPC 01'
]
g.set_terminal_states(terminal_states, allow_overlap=True)
g.compute_fate_probabilities()
fig = g.plot_fate_probabilities(legend_loc="right")
if fig is None:
    fig = plt.gcf()
mEPC_drivers = g.compute_lineage_drivers(lineages="mEPC 08")
mEPC_drivers.head(10)
model = cr.models.GAMR(adata)
fig = cr.pl.gene_trends(
    adata,
    model=model,
    data_key="counts",
    genes=["WDR54"],
    same_plot=True,
    ncols=2,
    time_key="velocity_pseudotime",
    hide_cells=True,
)
if fig is None:
    fig = plt.gcf()
output_path = os.path.join(os.getcwd(), "gene_trends.png")
fig.savefig(output_path, format="png", dpi=300)
plt.close(fig)
gmt_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/MSigDB_HALLMARK/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2024.1.Hs.gmt"
with open(gmt_file, "r") as f:
    line = f.readline().strip()

parts = line.split("\t")
gmt_genes = parts[2:]
matched_genes = [gene for gene in gmt_genes if gene in adata.var_names]
fig = cr.pl.heatmap(
    adata,
    model=model,
    data_key="counts",
    genes=matched_genes,
    lineages=["mEPC 08"],
    time_key="velocity_pseudotime",
    cbar=False,
    show_all_genes=True,
)
if fig is None:
    fig = plt.gcf()

output_path = os.path.join(os.getcwd(), "heatmap_emt.png")
fig.savefig(output_path, format="png", dpi=300)
plt.close(fig)
gmt_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/MSigDB_HALLMARK/HALLMARK_TGF_BETA_SIGNALING.v2024.1.Hs.gmt"
with open(gmt_file, "r") as f:
    line = f.readline().strip()

parts = line.split("\t")
gmt_genes = parts[2:]
matched_genes = [gene for gene in gmt_genes if gene in adata.var_names]
fig = cr.pl.heatmap(
    adata,
    model=model,
    data_key="counts",
    genes=matched_genes,
    lineages=["mEPC 08"],
    time_key="velocity_pseudotime",
    cbar=False,
    show_all_genes=True,
)
if fig is None:
    fig = plt.gcf()

output_path = os.path.join(os.getcwd(), "heatmap_TGFB.png")
fig.savefig(output_path, format="png", dpi=300)
plt.close(fig)

### mellon
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn.cluster import k_means
import palantir
import scanpy as sc
import warnings
from numba.core.errors import NumbaDeprecationWarning
warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
import mellon
import pandas as pd
import os
import palantir as pa
import scipy.sparse as sp
import seaborn as sns

work_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/mellon"
os.chdir(work_dir)
adata = sc.read("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/MNFHC_witch_counts.h5ad")

def override_with_raw_and_check_integers(adata):
    if adata.raw is None:
        return None
    if adata.shape == adata.raw.shape:
        adata.X = adata.raw.X.copy()
    else:
        current_genes = adata.var_names
        raw_genes = adata.raw.var_names
        common_genes = [g for g in current_genes if g in raw_genes]
        if len(common_genes) == 0:
            return None
        adata.X = adata.raw[:, common_genes].X.copy()
        adata._inplace_subset_var(common_genes)
    X = adata.X
    if sp.issparse(X):
        data = X.data
        is_all_integers = np.allclose(data, np.round(data))
    else:
        is_all_integers = np.allclose(X, np.round(X))
    return is_all_integers

result = override_with_raw_and_check_integers(adata)

dm_res = palantir.utils.run_diffusion_maps(adata, pca_key="X_umap", n_components=100)
ms_data = palantir.utils.determine_multiscale_space(adata)
imputed_X = palantir.utils.run_magic_imputation(adata)
start_cell = "Me_HN33_CACCAGGTCTTACCTA-1"
terminal_states = pd.Series(
    ["mEPC 08", "CAF 04", "CAF 06"],
    index=[
        "Pr_21240647_AACGCTTAAGTGGTCAGTCTGTCA",
        "Pr_22023021_ATCCTGTACTGGCATAATCATTCC",
        "Me_HN272_AGATTGCTCTTGACGA-12"
    ]
)
pr_res = pa.core.run_palantir(adata, start_cell, num_waypoints=1000, terminal_states=terminal_states, knn=400, n_jobs=16)
masks = palantir.presults.select_branch_cells(adata, q=0.1, eps=0.1)
model = mellon.DensityEstimator()
log_density = model.fit_predict(adata.obsm["DM_EigenVectors"])
predictor = model.predict
adata.obs["mellon_log_density"] = log_density
adata.obs["mellon_log_density_clipped"] = np.clip(log_density, *np.quantile(log_density, [0.05, 1]))

output_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/mellon/umap/branch_CAF06.png"
plt.figure()
palantir.plot.plot_branch(adata, branch_name="CAF 06", position="mellon_log_density", color="leiden_mEPC_CAF", s=10, linewidths=0.5, alpha=0.6)
plt.savefig(output_path, format="png", dpi=600, bbox_inches="tight")
plt.close()

output_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/mellon/umap/branch_CAF04.png"
plt.figure()
palantir.plot.plot_branch(adata, branch_name="CAF 04", position="mellon_log_density", color="leiden_mEPC_CAF", s=10, linewidths=0.5, alpha=0.6)
plt.savefig(output_path, format="png", dpi=600, bbox_inches="tight")
plt.close()

output_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/mellon/umap/branch_mEPC08.png"
plt.figure()
palantir.plot.plot_branch(adata, branch_name="mEPC 08", position="mellon_log_density", color="leiden_mEPC_CAF", s=10, linewidths=0.5, alpha=0.6)
plt.savefig(output_path, format="png", dpi=600, bbox_inches="tight")
plt.close()

### PAGA
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scanpy as sc
import os

sc.settings.verbosity = 3

adata = sc.read("/media/desk16/Essential_Data/Spatial_and_single_cell/MNFHC_witch_counts.h5ad")

def override_with_raw_and_check_integers(adata):
    if adata.raw is None:
        return None
    if adata.shape == adata.raw.shape:
        adata.X = adata.raw.X.copy()
    else:
        current_genes = adata.var_names
        raw_genes = adata.raw.var_names
        common_genes = [g for g in current_genes if g in raw_genes]
        if len(common_genes) == 0:
            return None
        adata.X = adata.raw[:, common_genes].X.copy()
        adata._inplace_subset_var(common_genes)
    X = adata.X
    if hasattr(X, "issparse") and X.issparse():
        data = X.data
        is_all_integers = np.allclose(data, np.round(data))
    else:
        is_all_integers = np.allclose(X, np.round(X))
    return is_all_integers

result = override_with_raw_and_check_integers(adata)

df = pd.read_csv("/media/desk16/Essential_Data/Spatial_and_single_cell/adata_partial_filtered_cellID.csv")
cell_index_list = df["cellID"].tolist()

adata = adata[cell_index_list].copy()
adata.X = adata.X.astype("float64")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.diffmap(adata)
adata.obsm["X_diffmap_"] = adata.obsm["X_diffmap"][:, 1:]
sc.pl.embedding(
    adata,
    "diffmap_",
    color=["leiden_mEPC_CAF"],
    show=False,
    size=80,
    alpha=0.6,
    ncols=1
)
save_path = os.path.join("diffmap_modified.png")
plt.savefig(save_path, dpi=300, bbox_inches="tight")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10, use_rep="X_diffmap_")
sc.tl.umap(adata)
sc.tl.paga(adata, groups="leiden_mEPC_CAF")
sc.pl.paga(
    adata,
    threshold=0,
    show=False,
    fontsize=0,
    node_size_scale=8,
    edge_width_scale=2
)
fig = plt.gcf()
fig.set_size_inches(16, 16)
save_path = "/media/desk16/Essential_Data/Spatial_and_single_cell/PAGA/paga.png"
plt.savefig(save_path, dpi=600)
plt.close()
sc.tl.draw_graph(adata, init_pos="paga")
sc.pl.paga_compare(
    adata,
    basis="X_umap",
    threshold=0,
    show=False,
    title="",
    edge_width_scale=2,
    legend_loc=None,
    node_size_scale=8,
    fontsize=0,
    alpha=0.5,
    size=240
)
fig = plt.gcf()
fig.set_size_inches(36, 16)
save_path = "/media/desk16/Essential_Data/Spatial_and_single_cell/PAGA/paga1.png"
plt.savefig(save_path, dpi=600)
plt.close()
adata.uns["iroot"] = np.flatnonzero(adata.obs["leiden_mEPC_CAF"] == "mEPC 04")[0]
sc.tl.dpt(adata)
adata_raw = sc.datasets.paul15()
sc.pp.log1p(adata_raw)
sc.pp.scale(adata_raw)
adata.raw = adata_raw
sc.pl.umap(
    adata,
    show=False,
    title="",
    color=["dpt_pseudotime"],
    legend_loc=None,
    alpha=0.5,
    size=240
)
fig = plt.gcf()
fig.set_size_inches(20, 16)
save_path = "/media/desk16/Essential_Data/Spatial_and_single_cell/PAGA/paga2.png"
plt.savefig(save_path, dpi=600)
plt.close()
sc.pl.paga_compare(
    adata,
    basis="X_umap",
    threshold=0,
    show=False,
    title="",
    color=["dpt_pseudotime"],
    edge_width_scale=2,
    legend_loc=None,
    node_size_scale=8,
    fontsize=0,
    alpha=0.5,
    size=240
)
fig = plt.gcf()
fig.set_size_inches(36, 16)
save_path = "/media/desk16/Essential_Data/Spatial_and_single_cell/PAGA/paga1_2.png"
plt.savefig(save_path, dpi=600)
plt.close()

### palantir
import palantir
import scanpy as sc
import numpy as np
import pandas as pd
import os
import palantir as pa
import scipy.sparse as sp
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

%matplotlib inline
sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [8, 8]
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['image.cmap'] = 'Spectral_r'
warnings.filterwarnings(action="ignore", module="matplotlib", message="findfont")
np.random.seed(5)

work_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/palantir"
os.chdir(work_dir)  

adata = sc.read("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/MNFHC_witch_counts.h5ad")

def override_with_raw_and_check_integers(adata):
    if adata.raw is None:
        return None
    
    if adata.shape == adata.raw.shape:
        adata.X = adata.raw.X.copy()
    else:
        current_genes = adata.var_names
        raw_genes = adata.raw.var_names
        common_genes = [g for g in current_genes if g in raw_genes]
        
        if len(common_genes) == 0:
            return None
        
        adata.X = adata.raw[:, common_genes].X.copy()
        adata._inplace_subset_var(common_genes)
    
    X = adata.X
    if sp.issparse(X):
        data = X.data
        is_all_integers = np.allclose(data, np.round(data))
    else:
        is_all_integers = np.allclose(X, np.round(X))
    
    return is_all_integers

result = override_with_raw_and_check_integers(adata)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.diffmap(adata)

fig, ax = plt.subplots(figsize=(6, 6))
sc.pl.embedding(adata, 'diffmap', color=['leiden_mEPC_CAF'], ax=ax, show=False, size=60, alpha=0.6)
save_path = os.path.join("diffmap_plot.png")
fig.savefig(save_path, dpi=300, bbox_inches='tight')

adata.obsm['X_diffmap_'] = adata.obsm['X_diffmap'][:, 1:]

sc.pl.embedding(
    adata,
    'diffmap_',
    color=['leiden_mEPC_CAF'],
    show=False,
    size=80,
    alpha=0.6,
    ncols=1
)
save_path = os.path.join("diffmap_modified.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10, use_rep='X_diffmap_')
sc.tl.umap(adata)
sc.pl.umap(
    adata, 
    color=['leiden_mEPC_CAF'],
    show=False,
    size=80,
    alpha=0.6,
    save="umap_leiden_mEPC_CAF.png"
)

plt.savefig(os.path.join(work_dir, "umap_leiden_mEPC_CAF.png"), dpi=300)

output_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/palantir"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

sc.tl.tsne(adata, use_rep='X_diffmap_', n_pcs=10)
sc.pl.tsne(
    adata,
    color=['leiden_mEPC_CAF'],
    show=False,
    size=80,
    alpha=0.6
)

fig = plt.gcf()
output_path = os.path.join(output_dir, "tsne_leiden_mEPC_CAF.png")
fig.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close(fig)

dm_res = palantir.utils.run_diffusion_maps(adata, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(adata)

imputed_X = palantir.utils.run_magic_imputation(adata)
palantir.plot.plot_diffusion_components(adata)

save_path = os.path.join(os.getcwd(), "diffusion_components.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

start_cell = "Me_HN33_CACCAGGTCTTACCTA-1"# 最幼稚亚群最幼稚细胞（根据scvelo_pseudotime）
terminal_states = pd.Series(
    ["mEPC 08", "CAF 04", "CAF 06"],
    index=[
        "Pr_21240647_AACGCTTAAGTGGTCAGTCTGTCA",
        "Pr_22023021_ATCCTGTACTGGCATAATCATTCC",
        "Me_HN272_AGATTGCTCTTGACGA-12"
    ]
) # 目标亚群的最成熟细胞（根据scvelo_pseudotime）
pr_res = pa.core.run_palantir(
    adata, start_cell, num_waypoints=500, terminal_states=terminal_states, knn=450, 
    n_jobs=16
)
pa.plot.plot_palantir_results(adata, s=10)
plt.savefig('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/palantir/palantir_results.png', format='png', dpi=300)
plt.close()
masks = palantir.presults.select_branch_cells(adata, q=0.0001, eps=0.0001)
pa.plot.plot_branch_selection(adata)
plt.savefig('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/palantir/branch_selection.png', format='png', dpi=300)
plt.close()

### sctour
import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy.sparse as sp

os.chdir("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/sctour")
adata = sc.read("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/MNFHC_witch_counts.h5ad")

def override_with_raw_and_check_integers(adata):
    if adata.raw is None:
        return None
    if adata.shape == adata.raw.shape:
        adata.X = adata.raw.X.copy()
    else:
        current_genes = adata.var_names
        raw_genes = adata.raw.var_names
        common_genes = [g for g in current_genes if g in raw_genes]
        if len(common_genes) == 0:
            return None
        adata.X = adata.raw[:, common_genes].X.copy()
        adata._inplace_subset_var(common_genes)
    X = adata.X
    if sp.issparse(X):
        data = X.data
        is_all_integers = np.allclose(data, np.round(data))
    else:
        is_all_integers = np.allclose(X, np.round(X))
    return is_all_integers

result = override_with_raw_and_check_integers(adata)

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=5000, subset=True)

if sp.issparse(adata.X):
    adata.X = adata.X.toarray()
adata.X = adata.X.astype('float32')

tnode = sct.train.Trainer(
    adata,
    loss_mode='nb',
    alpha_recon_lec=0.5,
    alpha_recon_lode=0.5,
    use_gpu=False
)
tnode.train()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10, use_rep='X_diffmap_')
sc.tl.umap(adata)

adata.obs['ptime'] = tnode.get_time()
mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
adata.obsm['X_TNODE'] = mix_zs
adata.obsm['X_VF'] = tnode.get_vector_field(
    adata.obs['ptime'].values,
    adata.obsm['X_TNODE']
)
adata = adata[np.argsort(adata.obs['ptime'].values), :]

fig1, ax1 = plt.subplots(figsize=(10, 10))
sc.pl.umap(
    adata,
    color='leiden_mEPC_CAF',
    ax=ax1,
    legend_loc='on data',
    show=False,
    frameon=False
)
output_path1 = os.path.join(os.getcwd(), "umap_leiden_mEPC_CAF.png")
plt.savefig(output_path1, dpi=300, bbox_inches='tight')
plt.close(fig1)

fig2, ax2 = plt.subplots(figsize=(10, 10))
sc.pl.umap(
    adata,
    color='Batch',
    ax=ax2,
    show=False,
    frameon=False
)
output_path2 = os.path.join(os.getcwd(), "umap_Batch.png")
plt.savefig(output_path2, dpi=300, bbox_inches='tight')
plt.close(fig2)

fig3, ax3 = plt.subplots(figsize=(8, 6))
sc.pl.umap(
    adata,
    color='ptime',
    ax=ax3,
    show=False,
    frameon=False,
    size=80,
    alpha=0.6
)
output_path3 = os.path.join(os.getcwd(), "umap_ptime.png")
plt.savefig(output_path3, dpi=300, bbox_inches='tight')
plt.close(fig3)

fig4, ax4 = plt.subplots(figsize=(9.5, 8))
sct.vf.plot_vector_field(
    adata,
    zs_key='X_TNODE',
    vf_key='X_VF',
    use_rep_neigh='X_TNODE',
    color='leiden_mEPC_CAF',
    show=False,
    ax=ax4,
    legend_loc='none',
    frameon=False,
    size=100,
    alpha=0.2,
    linewidth=3
)
output_path4 = os.path.join(os.getcwd(), "vector_field.png")
plt.savefig(output_path4, dpi=300, bbox_inches='tight')
plt.close(fig4)

### scvelo
import scvelo as scv
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import CytoSimplex as csx
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
file_path_h5ad = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/MNFHC_for_velo.h5ad'
adata = sc.read(file_path_h5ad)
print(adata)
scv.pl.proportions(adata, save='/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/figures/proportions.png', figsize=(10, 6), dpi=300)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata, n_jobs=20)
leiden_colors = {'CAF 01': '#1F77B4','CAF 02': '#FF7E0E','CAF 03': '#279A68','CAF 04': '#D62728','CAF 05': '#AA40FC','CAF 06': '#8C564B','mEPC 01': '#E377C2','mEPC 02': '#B5BD61','mEPC 03': '#17BEDF','mEPC 04': '#AEC7E8','mEPC 05': '#FFBB78','mEPC 06': '#98DF8A','mEPC 07': '#FF9896','mEPC 08': '#C5B0D5','mEPC 09': '#C49C94','mEPC 10': '#F7B6D2'}
scv.pl.velocity_embedding_stream(adata, basis='umap', save='/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/figures/velocity_embedding_stream.png', figsize=(10, 8), dpi=300, density=3, linewidth=2, color='leiden_mEPC_CAF', palette=leiden_colors, size=300)
scv.pl.velocity_embedding_stream(adata, basis='umap', save='/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/figures/velocity_embedding_stream_nolegend.png', figsize=(10, 8), dpi=300, density=3, linewidth=2, color='leiden_mEPC_CAF', palette=leiden_colors, size=300, legend_loc='none')
scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=200, dpi=300, save='/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/figures/velocity_embedding.png', figsize=(10, 8), color='leiden_mEPC_CAF', palette=leiden_colors, size=300, legend_loc='none')
scv.pl.velocity_embedding_grid(adata, basis='umap', arrow_length=1.5, arrow_size=10, figsize=(10, 8), color='leiden_mEPC_CAF', save='/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/figures/embedding_grid.png', title='', palette=leiden_colors, size=300, dpi=300)
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, size=100, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], figsize=(10, 8), save='/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/figures/velocity_confidence.png')
scv.pl.velocity_graph(adata, threshold=0.1, figsize=(10, 8), color='leiden_mEPC_CAF', palette=leiden_colors, size=20, dpi=300, legend_loc='none', save='/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/figures/velocity_transition.png')
scv.tl.velocity_pseudotime(adata)
velocity_pseudotime_df = adata.obs[['velocity_pseudotime']]
velocity_pseudotime_df.to_csv('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/velocity_pseudotime.csv')
loaded_df = pd.read_csv('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/velocity_pseudotime.csv')
print("前五行：")
print(loaded_df.head())
print("\n所有列名：")
print(loaded_df.columns)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', size=100, dpi=300, figsize=(10, 8), save='/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/figures/velocity_pseudotime.png')
scv.pl.velocity(adata, ['WDR54'], ncols=1, save='/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/figures/velocity_genes.png', dpi=300, figsize=(30,24))
vertices = ['mEPC 08', 'CAF 06', 'CAF 04']
vertex_colors = ['#C5B0D5', '#8C564B', '#D62728']
features = csx.select_top_features(
    adata, 
    cluster_var="leiden_mEPC_CAF", 
    vertices=vertices
)
cluster_color_map = {
    "CAF 01": "#1F77B4",
    "CAF 02": "#FF7E0E",
    "CAF 03": "#279A68",
    "CAF 04": "#D62728",
    "CAF 05": "#AA40FC",
    "CAF 06": "#8C564B",
    "mEPC 01": "#E377C2",
    "mEPC 02": "#B5BD61",
    "mEPC 03": "#17BEDF",
    "mEPC 04": "#AEC7E8",
    "mEPC 05": "#FFBB78",
    "mEPC 06": "#98DF8A",
    "mEPC 07": "#FF9896",
    "mEPC 08": "#C5B0D5",
    "mEPC 09": "#C49C94",
    "mEPC 10": "#F7B6D2"
}
dot_colors = adata.obs["leiden_mEPC_CAF"].map(cluster_color_map)
csx.plot_ternary(
    x=adata,
    cluster_var="leiden_mEPC_CAF",
    vertices=vertices,
    features=features,
    velo_graph="velocity_graph",
    save_fig=True,
    fig_path="/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2C/figures/ternary_plot_leiden_mEPC_CAF.png",
    fig_size=(24, 18),
    processed=False,
    method='euclidean',
    sigma=0.08,
    scale=True,
    title="Ternary Simplex Plot based on leiden_mEPC_CAF",
    split_cluster=False,
    dot_color=dot_colors,  
    n_velogrid=10,
    radius=0.08,
    dot_size=30,
    vertex_colors=vertex_colors,
    vertex_label_size=12,
    gridline_alpha=0.4,
    axis_text_show=True,
    arrow_linewidth=0.004
)

### stavia
import os
import scanpy as sc
import omicverse as ov
from omicverse.externel import VIA
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
np.int = int
os.chdir("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/StaVIA")
adata1 = sc.read_h5ad('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/MNFHC_witch_counts.h5ad')
df = pd.read_csv('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/adata_partial_filtered_cellID.csv')
cell_index_list = df['cellID'].tolist()
adata = adata1[cell_index_list].copy()
adata = ov.pp.preprocess(adata, mode='shiftlog|pearson', n_HVGs=2000)
hvgs = adata.var_names[adata.var['highly_variable']].tolist()
adata = adata[:, hvgs]
ov.pp.scale(adata)
ncomps = 30 
knn = 15
v0_random_seed = 4
memory = 10
dataset = ''
use_rep = 'X_harmony'
clusters = 'leiden_mEPC_CAF'
v0 = VIA.core.VIA(
    data=adata.obsm[use_rep][:, 0:ncomps],
    true_label=adata.obs[clusters],
    edgepruning_clustering_resolution=0.15,
    cluster_graph_pruning=0.15,
    knn=knn,
    resolution_parameter=1.5,
    dataset=dataset,
    random_seed=v0_random_seed,
    memory=memory,
    root_user=['mEPC 08']
)
v0.run_VIA()
color_mapping = {
    "CAF 01": "#1F77B4",
    "CAF 02": "#FF7E0E",
    "CAF 03": "#279A68",
    "CAF 04": "#D62728",
    "CAF 05": "#AA40FC",
    "CAF 06": "#8C564B",
    "mEPC 01": "#E377C2",
    "mEPC 02": "#B5BD61",
    "mEPC 03": "#17BEDF",
    "mEPC 04": "#AEC7E8",
    "mEPC 05": "#FFBB78",
    "mEPC 06": "#98DF8A",
    "mEPC 07": "#FF9896",
    "mEPC 08": "#C5B0D5",
    "mEPC 09": "#C49C94",
    "mEPC 10": "#F7B6D2",
}
if not pd.api.types.is_categorical_dtype(adata.obs["leiden_mEPC_CAF"]):
    adata.obs["leiden_mEPC_CAF"] = adata.obs["leiden_mEPC_CAF"].astype("category")
cluster_order = list(adata.obs["leiden_mEPC_CAF"].cat.categories)
color_list = [color_mapping.get(cl, "#000000") for cl in cluster_order]
adata.uns["leiden_mEPC_CAF_colors"] = color_list
fig, ax, ax1 = VIA.core.plot_piechart_viagraph_ov(
    adata,
    clusters="leiden_mEPC_CAF",
    dpi=300,
    via_object=v0,
    ax_text=False,
    show_legend=False
)
fig.set_size_inches(16, 8)
output_path = os.path.join(os.getcwd(), "viagraph_piechart.png")
fig.savefig(output_path, format='png', dpi=300)
plt.close(fig)

plt.ioff()
adata.obs['pt_via'] = v0.single_cell_pt_markov
ax = ov.pl.embedding(
    adata,
    basis='X_umap',
    color=['pt_via'],
    frameon='small',
    cmap='Reds',
    show=False
)
fig = ax.get_figure()
fig.set_size_inches(9.5, 8)
output_path = os.path.join(os.getcwd(), "embedding_plot.png")
fig.savefig(output_path, format='png', dpi=300)
plt.close(fig)
plt.ioff()
clusters = 'leiden_mEPC_CAF'
color_true_list = adata.uns['{}_colors'.format(clusters)]
fig, ax, ax1 = VIA.core.plot_trajectory_curves_ov(
    adata,
    clusters='leiden_mEPC_CAF',
    dpi=300,
    via_object=v0,
    embedding=adata.obsm['X_umap'],
    draw_all_curves=False
)
fig.set_size_inches(21, 8)
output_path = os.path.join(os.getcwd(), "trajectory_curves.png")
fig.savefig(output_path, format='png', dpi=300)
plt.close(fig)
plt.ioff()
v0.embedding = adata.obsm['X_umap']
fig, ax = VIA.core.plot_atlas_view(
    via_object=v0, 
    n_milestones=500, 
    sc_labels=adata.obs['leiden_mEPC_CAF'], 
    fontsize_title=12,
    fontsize_labels=12,
    dpi=300,
    extra_title_text='Atlas View colored by pseudotime'
)
fig.set_size_inches(9.5, 8)
output_path = os.path.join(os.getcwd(), "atlas_view_nocell.png")
fig.savefig(output_path, format='png', dpi=300)
plt.close(fig)
decay = 0.6
i_bw = 0.02
global_visual_pruning = 0.5
n_milestones = 500
v0.hammerbundle_milestone_dict = VIA.core.make_edgebundle_milestone(
    via_object=v0, 
    n_milestones=n_milestones, 
    decay=decay, 
    initial_bandwidth=i_bw,
    global_visual_pruning=global_visual_pruning
)

fig_atlas, ax_atlas = VIA.core.plot_atlas_view(
    via_object=v0,  
    add_sc_embedding=True, 
    sc_labels_expression=adata.obs['leiden_mEPC_CAF'], 
    cmap='jet',
    sc_labels=adata.obs['leiden_mEPC_CAF'], 
    text_labels=False, 
    extra_title_text='Atlas View by Cell type', 
    fontsize_labels=3,
    fontsize_title=3,
    dpi=300
)
fig_atlas.set_size_inches(9.5, 8)
output_path_atlas = os.path.join(os.getcwd(), "atlas_view.png")
fig_atlas.savefig(output_path_atlas, format='png', dpi=300)
plt.close(fig_atlas)
fig_stream, ax_stream = VIA.core.via_streamplot_ov(
    adata, 
    'leiden_mEPC_CAF',
    v0,  
    embedding=adata.obsm['X_umap'], 
    dpi=300,
    density_grid=0.8, 
    density_stream=3,
    arrow_size=1,
    scatter_size=30, 
    scatter_alpha=0.3, 
    linewidth=2
)
fig_stream.set_size_inches(9.5, 8)
for txt in ax_stream.texts:
    txt.set_fontsize(12)
    txt.set_ha('center')
    txt.set_va('center')
output_path_stream = os.path.join(os.getcwd(), "stream_plot.png")
fig_stream.savefig(output_path_stream, format='png', dpi=300)
plt.close(fig_stream)
fig, ax = VIA.core.via_streamplot_ov(
    adata, 
    'clusters',
    v0,
    density_grid=0.8, 
    scatter_size=30, 
    color_scheme='time', 
    linewidth=2, 
    min_mass=1, 
    cutoff_perc=5, 
    scatter_alpha=0.3, 
    marker_edgewidth=0.1, 
    density_stream=3, 
    smooth_transition=1, 
    smooth_grid=0.5, 
    arrow_size=1,
    dpi=300
)
fig.set_size_inches(9.5, 8)
output_path = os.path.join(os.getcwd(), "streamplot_hearmap.png")
fig.savefig(output_path, format="png", dpi=80)
plt.close(fig)
plt.ioff()
gene_list_magic = ['WDR54', 'C1QC']
df = adata.to_df()
df_magic = v0.do_impute(df, magic_steps=3, gene_list=gene_list_magic)
fig, axs = VIA.core.plot_viagraph(
    via_object=v0, 
    type_data='gene', 
    df_genes=df_magic, 
    gene_list=gene_list_magic, 
    arrow_head=0.1,
    label_text=False
)
fig.set_size_inches(32, 8)
output_path = os.path.join(os.getcwd(), "viagraph.png")
fig.savefig(output_path, format='png', dpi=300)
plt.close(fig)

### CytoTrace
library(data.table)
library(doParallel)
library(dplyr)
library(ggplot2)
library(magrittr)
library(Matrix)
library(plyr)
library(SeuratObject)
library(stringr)
library(HiClimR)
library(CytoTRACE2)
install_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/Rpackage"
library(Seurat, lib.loc = install_path)
packageVersion("Seurat")
seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2D/vector/MNFHC_for_vector.rds")
cytotrace2_result <- cytotrace2(
  input = seurat_obj,        
  species = "human",         
  slot_type = "counts",        
  is_seurat = TRUE           
)
output_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/cytotrace/cytotrace_meta_data.csv"
write.csv(cytotrace2_result@meta.data, file = output_path, row.names = TRUE)
file_path <- '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/cytotrace/cytotrace_meta_data.csv'
cytotrace_data <- read.csv(file_path)
colnames(cytotrace_data)[1] <- "CellID"
write.csv(cytotrace_data, file = file_path, row.names = FALSE)


### monocle
### monocle2
library(monocle)
library(Seurat)
library(argparse)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(patchwork)
set.seed(10)
## 第一部分——细胞全景
outputDir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/monocle/monocle2/Allcell"
if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
}
print(paste("输出目录已设置为:", outputDir))
seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/vector/MNFHC_for_vector.rds")
Idents(seurat_obj) <- seurat_obj@meta.data$leiden_mEPC_CAF
full_categories <- c("mEPC 08", "CAF 04", "CAF 06")
categories <- unique(seurat_obj@meta.data$leiden_mEPC_CAF)
selected_cells <- unlist(lapply(categories, function(cat) {
  cells_in_category <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$leiden_mEPC_CAF == cat, ])
    if (cat %in% full_categories) {
    return(cells_in_category)
  } else {
    return(sample(cells_in_category, size = length(cells_in_category) %/% 3))
  }
}))
seurat_obj <- subset(seurat_obj, cells = selected_cells) # 降低数据规模，避免monocle2的Bug 
DefaultAssay(seurat_obj) <- 'RNA'
pd <- new('AnnotatedDataFrame', data = seurat_obj@meta.data)
fData <- data.frame(gene_short_name = row.names(seurat_obj), row.names = row.names(seurat_obj))
fd <- new('AnnotatedDataFrame', data = fData)
gbm_cds <- newCellDataSet(
  as(as.matrix(seurat_obj@assays[["RNA"]]@counts), 'sparseMatrix'), 
  phenoData = pd,  
  featureData = fd,  
  lowerDetectionLimit = 0.5,  
  expressionFamily = negbinomial.size()  
)
gbm_cds <- estimateSizeFactors(gbm_cds)
gbm_cds <- estimateDispersions(gbm_cds)
expressed_genes <- row.names(subset(fData(gbm_cds)))
express_genes <- VariableFeatures(seurat_obj)
gbm_cds <- setOrderingFilter(gbm_cds[expressed_genes,], express_genes)
png(file = paste(outputDir, "all_ordering_genes.png", sep = '/'), width = 800, height = 600)
print(plot_ordering_genes(gbm_cds))  
dev.off()  
tryCatch({
  
  gbm_cds <- reduceDimension(gbm_cds, max_components = 2, method = 'DDRTree')
}, error = function(e) {
  print("Error in reduceDimension(). Try to apply auto_param_selection = F")
  gbm_cds <- reduceDimension(gbm_cds, max_components = 2, method = 'DDRTree', auto_param_selection = F)
})
gbm_cds <- orderCells(gbm_cds)  
gbm_cds <- orderCells(gbm_cds, root_state = 1)  
gbm_cds$barcode <- colnames(gbm_cds)
rep <- as.data.frame(pData(gbm_cds))
write.table(pData(gbm_cds), file = paste0(outputDir,"/","cell_Pseudotime.txt"), row.names = T, quote = F, sep = ',')
state_levels <- levels(pData(gbm_cds)$State)
if (length(state_levels) <= 2) {  
  widths_state = 8  
  heights_state = 5  
  nrows_state = 1
} else if(1 < length(state_levels) & length(state_levels) <= 8) {  
  widths_state = length(state_levels) * 1.5  
  heights_state = length(state_levels) * 1.5  
  nrows_state = 2
} else {  
  widths_state = length(state_levels) * 1.5  
  heights_state = length(state_levels) * 1.5  
  nrows_state = 3
}
celltype_levels <- levels(gbm_cds$celltype)  
if (length(celltype_levels) <= 2) {  
  widths_celltype = 8  
  heights_celltype = 5  
  nrows_celltype = 1
} else if (1 < length(celltype_levels) & length(celltype_levels) <= 8) {  
  widths_celltype = length(celltype_levels) * 1.5  
  heights_celltype = length(celltype_levels) * 1.5  
  nrows_celltype = 2
} else {  
  widths_celltype = length(celltype_levels) * 1.5  
  heights_celltype = length(celltype_levels) * 1.5  
  nrows_celltype = 3
}
sample_levels <- levels(gbm_cds$orig.ident)  
if (length(celltype_levels) <= 2) {  
  widths_orig = 8  
  heights_orig = 5  
  nrows_orig = 1
} else if (1 < length(celltype_levels) & length(celltype_levels) <= 8) {  
  widths_orig = length(celltype_levels) * 1.5  
  heights_orig = length(celltype_levels) * 1.5  
  nrows_orig = 2
} else {  
  widths_orig = length(celltype_levels) * 1.5  
  heights_orig = length(celltype_levels) * 1.5  
  nrows_orig = 3
}
p1 <- plot_cell_trajectory(gbm_cds, color_by = "Pseudotime", cell_size = 3) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "none")

png(file = paste0(outputDir, "/1_Pseudotime_trajectory.png"), width = 800, height = 600)
print(p1)  
dev.off()  
p1 <- plot_cell_trajectory(gbm_cds, color_by = "Pseudotime", cell_size = 3) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top")

png(file = paste0(outputDir, "/1_Pseudotime_trajectory_legend.png"), width = 800, height = 600)
print(p1)  
dev.off()  
p1 <- plot_cell_trajectory(gbm_cds, color_by = "Pseudotime", cell_size = 3) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top") +
      facet_wrap(~leiden_mEPC_CAF)
png(file = paste0(outputDir, "/2_Pseudotime_trajectory_splited_by_orig.png"), width = 1600, height = 1200)
print(p1)  
dev.off()  
p1 <- plot_cell_trajectory(gbm_cds, color_by = "Pseudotime") + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top") +
      facet_wrap(~group)
png(file = paste0(outputDir, "/3_Pseudotime_trajectory_splited_by_group.png"), width = 800, height = 600)
print(p1)  
dev.off()  
MYCOLOR <- c(
  "#1F77B4", "#FF7E0E", "#279A68", "#D62728", "#AA40FC", "#8C564B",  
  "#E377C2", "#B5BD61", "#17BEDF", "#AEC7E8", "#FFBB78", "#98DF8A",  
  "#FF9896", "#C5B0D5", "#C49C94", "#F7B6D2"                       
)
p3 <- plot_cell_trajectory(gbm_cds, color_by = 'leiden_mEPC_CAF', cell_size = 3) + 
  scale_color_manual(values = scales::alpha(MYCOLOR, 0.6)) +  
  ggtitle('celltype') + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
png(file = paste0(outputDir, "/4_Celltype_trajectory.png"), width = 800, height = 600)
print(p3)  
dev.off()  
p3 <- plot_cell_trajectory(gbm_cds, color_by = 'leiden_mEPC_CAF', cell_size = 0.5, cell_link_size = 0.5) + 
  facet_wrap(~leiden_mEPC_CAF, nrow = nrows_celltype, scales = "free") + 
  scale_color_manual(values = MYCOLOR) + 
  ggtitle('leiden_mEPC_CAF') + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
png(file = paste0(outputDir, "/5_Celltype_trajectory_splited_by_celltype.png"), 
    width = 3000, height = 300)
print(p3)  
dev.off()  
p3 <- ggplot(rep, aes(Pseudotime, fill = leiden_mEPC_CAF)) + 
  geom_density() + 
  theme_bw() + 
  RotatedAxis() + 
  theme(
    strip.text = element_blank(),
    strip.background = element_rect(color = "white", fill = "white"),
    panel.grid = element_blank()
  ) + 
  scale_fill_manual(values = MYCOLOR)
png(file = paste0(outputDir, "/6_Density_trajectory.png"), 
    height = 600, width = 800)
print(p3)  
dev.off() 
outputDir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/monocle/monocle2/3cell"
if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)  
}
print(paste("输出目录已设置为:", outputDir))
seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2D/vector/MNFHC_for_vector.rds")
Idents(seurat_obj) <- seurat_obj@meta.data$leiden_mEPC_CAF  
full_categories <- c("mEPC 08", "CAF 04", "CAF 06")
seurat_obj <- subset(seurat_obj, subset = leiden_mEPC_CAF %in% full_categories)
DefaultAssay(seurat_obj) <- 'RNA'
pd <- new('AnnotatedDataFrame', data = seurat_obj@meta.data)
fData <- data.frame(gene_short_name = row.names(seurat_obj), row.names = row.names(seurat_obj))
fd <- new('AnnotatedDataFrame', data = fData)
gbm_cds <- newCellDataSet(
  as(as.matrix(seurat_obj@assays[["RNA"]]@counts), 'sparseMatrix'),  # 将Seurat中的RNA计数数据转换为稀疏矩阵
  phenoData = pd, 
  featureData = fd,  
  lowerDetectionLimit = 0.5,  
  expressionFamily = negbinomial.size()  
)
gbm_cds <- estimateSizeFactors(gbm_cds)
gbm_cds <- estimateDispersions(gbm_cds)
expressed_genes <- row.names(subset(fData(gbm_cds)))
express_genes <- VariableFeatures(seurat_obj)
gbm_cds <- setOrderingFilter(gbm_cds[expressed_genes,], express_genes)
png(file = paste(outputDir, "all_ordering_genes.png", sep = '/'), width = 800, height = 600)
print(plot_ordering_genes(gbm_cds))  
dev.off()  
tryCatch({
  gbm_cds <- reduceDimension(gbm_cds, max_components = 2, method = 'DDRTree')
}, error = function(e) {
  print("Error in reduceDimension(). Try to apply auto_param_selection = F")
  gbm_cds <- reduceDimension(gbm_cds, max_components = 2, method = 'DDRTree', auto_param_selection = F)
})
gbm_cds <- orderCells(gbm_cds)  
gbm_cds <- orderCells(gbm_cds, root_state = 1)  
gbm_cds$barcode <- colnames(gbm_cds)
rep <- as.data.frame(pData(gbm_cds))
write.table(pData(gbm_cds), file = paste0(outputDir,"/","cell_Pseudotime.txt"), row.names = T, quote = F, sep = ',')
state_levels <- levels(pData(gbm_cds)$State)
if (length(state_levels) <= 2) {  
  widths_state = 8  
  heights_state = 5  
  nrows_state = 1
} else if(1 < length(state_levels) & length(state_levels) <= 8) {  
  widths_state = length(state_levels) * 1.5  
  heights_state = length(state_levels) * 1.5  
  nrows_state = 2
} else {  
  widths_state = length(state_levels) * 1.5  
  heights_state = length(state_levels) * 1.5  
  nrows_state = 3
}
celltype_levels <- levels(gbm_cds$celltype)  
if (length(celltype_levels) <= 2) {  
  widths_celltype = 8  
  heights_celltype = 5  
  nrows_celltype = 1
} else if (1 < length(celltype_levels) & length(celltype_levels) <= 8) {  
  widths_celltype = length(celltype_levels) * 1.5  
  heights_celltype = length(celltype_levels) * 1.5  
  nrows_celltype = 2
} else {  
  widths_celltype = length(celltype_levels) * 1.5  
  heights_celltype = length(celltype_levels) * 1.5  
  nrows_celltype = 3
}
sample_levels <- levels(gbm_cds$orig.ident)  
if (length(celltype_levels) <= 2) {  
  widths_orig = 8  
  heights_orig = 5  
  nrows_orig = 1
} else if (1 < length(celltype_levels) & length(celltype_levels) <= 8) {  
  widths_orig = length(celltype_levels) * 1.5  
  heights_orig = length(celltype_levels) * 1.5  
  nrows_orig = 2
} else {  
  widths_orig = length(celltype_levels) * 1.5  
  heights_orig = length(celltype_levels) * 1.5  
  nrows_orig = 3
}
p1 <- plot_cell_trajectory(gbm_cds, color_by = "Pseudotime", cell_size = 3) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "none")

png(file = paste0(outputDir, "/7_Pseudotime_trajectory_3cell.png"), width = 800, height = 600)
print(p1)  
dev.off()  
p1 <- plot_cell_trajectory(gbm_cds, color_by = "Pseudotime", cell_size = 3) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top")


png(file = paste0(outputDir, "/7_Pseudotime_trajectory_3cell_legend.png"), width = 800, height = 600)
print(p1)  
dev.off()  
p1 <- plot_cell_trajectory(gbm_cds, color_by = "Pseudotime", cell_size = 3) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top") +
      facet_wrap(~leiden_mEPC_CAF)

png(file = paste0(outputDir, "/8_Pseudotime_trajectory_splited_by_orig_3cell.png"), width = 1600, height = 1200)
print(p1)  
dev.off()  
MYCOLOR <- c(
  "#D62728", "#8C564B",
  "#C5B0D5"                       
)
p3 <- plot_cell_trajectory(gbm_cds, color_by = 'leiden_mEPC_CAF', cell_size = 3) + 
  scale_color_manual(values = scales::alpha(MYCOLOR, 0.6)) +  
  ggtitle('celltype') + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
png(file = paste0(outputDir, "/9_Celltype_trajectory.png"), width = 800, height = 600)
print(p3)  
dev.off()  
p3 <- plot_cell_trajectory(gbm_cds, color_by = 'leiden_mEPC_CAF', cell_size = 3) + 
  scale_color_manual(values = scales::alpha(MYCOLOR, 0.6)) +  
  ggtitle('celltype') + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"  
  )
png(file = paste0(outputDir, "/9_Celltype_trajectory_nolegend.png"), width = 800, height = 600)
print(p3)  
dev.off()  
p3 <- ggplot(rep, aes(Pseudotime, fill = leiden_mEPC_CAF)) + 
  geom_density() + 
  theme_bw() + 
  RotatedAxis() + 
  theme(
    strip.text = element_blank(),
    strip.background = element_rect(color = "white", fill = "white"),
    panel.grid = element_blank()
  ) + 
  scale_fill_manual(values = MYCOLOR)
png(file = paste0(outputDir, "/10_Density_trajectory.png"), 
    height = 600, width = 800)
print(p3)  
dev.off()  
## 第三部分——WDR54表达
outputDir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/monocle/monocle2/WDR54"
if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)  
}
print(paste("输出目录已设置为:", outputDir))
set.seed(123)
seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/vector/MNFHC_for_vector.rds")
Idents(seurat_obj) <- seurat_obj@meta.data$leiden_mEPC_CAF  
wdr54_expression <- seurat_obj[["RNA"]]@counts["WDR54"]
selected_cells <- names(wdr54_expression[wdr54_expression > 0])
cat("提取的细胞数量：", length(selected_cells), "\n")
seurat_obj <- subset(seurat_obj, cells = selected_cells)
DefaultAssay(seurat_obj) <- 'RNA'
pd <- new('AnnotatedDataFrame', data = seurat_obj@meta.data)
fData <- data.frame(gene_short_name = row.names(seurat_obj), row.names = row.names(seurat_obj))
fd <- new('AnnotatedDataFrame', data = fData)
gbm_cds <- newCellDataSet(
  as(as.matrix(seurat_obj@assays[["RNA"]]@data), 'sparseMatrix'),  
  phenoData = pd,  
  featureData = fd,  
  lowerDetectionLimit = 0.5,  
  expressionFamily = negbinomial.size()  
)
gbm_cds <- estimateSizeFactors(gbm_cds)
gbm_cds <- estimateDispersions(gbm_cds)
expressed_genes <- row.names(subset(fData(gbm_cds)))
disp_table <- dispersionTable(gbm_cds)
ordering_genes_temp <- subset(disp_table, mean_expression >= 0.005 & dispersion_empirical >= 1 * dispersion_fit)
ordering_genes <- ordering_genes_temp$gene_id
gbm_cds <- setOrderingFilter(gbm_cds, ordering_genes)
write.table(ordering_genes, file = paste0(outputDir, "/", "all_dispersion.ordering_genes.tsv"), 
            row.names = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)
png(file = paste(outputDir, "all_ordering_genes.png", sep = '/'), width = 800, height = 600)
print(plot_ordering_genes(gbm_cds))  
dev.off()  
tryCatch({
  gbm_cds <- reduceDimension(gbm_cds, max_components = 2, method = 'DDRTree')
}, error = function(e) {
  print("Error in reduceDimension(). Try to apply auto_param_selection = F")
  gbm_cds <- reduceDimension(gbm_cds, max_components = 2, method = 'DDRTree', auto_param_selection = F)
})
gbm_cds <- orderCells(gbm_cds)  
gbm_cds <- orderCells(gbm_cds, root_state = 1)  
gbm_cds$barcode <- colnames(gbm_cds)
rep <- as.data.frame(pData(gbm_cds))
write.table(pData(gbm_cds), file = paste0(outputDir,"/","cell_Pseudotime.txt"), row.names = T, quote = F, sep = ',')
state_levels <- levels(pData(gbm_cds)$State)
if (length(state_levels) <= 2) {  
  widths_state = 8  
  heights_state = 5  
  nrows_state = 1
} else if(1 < length(state_levels) & length(state_levels) <= 8) {  
  widths_state = length(state_levels) * 1.5  
  heights_state = length(state_levels) * 1.5  
  nrows_state = 2
} else {  
  widths_state = length(state_levels) * 1.5  
  heights_state = length(state_levels) * 1.5  
  nrows_state = 3
}
celltype_levels <- levels(gbm_cds$celltype)  
if (length(celltype_levels) <= 2) {  
  widths_celltype = 8  
  heights_celltype = 5  
  nrows_celltype = 1
} else if (1 < length(celltype_levels) & length(celltype_levels) <= 8) {  
  widths_celltype = length(celltype_levels) * 1.5  
  heights_celltype = length(celltype_levels) * 1.5  
  nrows_celltype = 2
} else {  
  widths_celltype = length(celltype_levels) * 1.5  
  heights_celltype = length(celltype_levels) * 1.5  
  nrows_celltype = 3
}
sample_levels <- levels(gbm_cds$orig.ident)  
if (length(celltype_levels) <= 2) {  
  widths_orig = 8  
  heights_orig = 5  
  nrows_orig = 1
} else if (1 < length(celltype_levels) & length(celltype_levels) <= 8) {  
  widths_orig = length(celltype_levels) * 1.5  
  heights_orig = length(celltype_levels) * 1.5  
  nrows_orig = 2
} else {  
  widths_orig = length(celltype_levels) * 1.5  
  heights_orig = length(celltype_levels) * 1.5  
  nrows_orig = 3
}
color_mapping <- c(
  "CAF 01" = "#1F77B4",
  "CAF 02" = "#FF7E0E",
  "CAF 03" = "#279A68",
  "CAF 04" = "#D62728",
  "CAF 05" = "#AA40FC",
  "CAF 06" = "#8C564B",
  "mEPC 01" = "#E377C2",
  "mEPC 02" = "#B5BD61",
  "mEPC 03" = "#17BEDF",
  "mEPC 04" = "#AEC7E8",
  "mEPC 05" = "#FFBB78",
  "mEPC 06" = "#98DF8A",
  "mEPC 07" = "#FF9896",
  "mEPC 08" = "#C5B0D5",
  "mEPC 09" = "#C49C94",
  "mEPC 10" = "#F7B6D2"
)
genes_of_interest <- row.names(subset(fData(gbm_cds), gene_short_name %in% c("WDR54", "CD8A"))) 
if(length(genes_of_interest) == 0) {
    stop("没有找到相关基因。请检查基因名称是否正确。")
}
outputFile <- file.path(outputDir, "gene_expression_branched_pseudotime_all.png")
gbm_cds$leiden_mEPC_CAF <- factor(gbm_cds$leiden_mEPC_CAF, levels = names(color_mapping))
png(outputFile, width = 600, height = 600)
plot_genes_in_pseudotime(gbm_cds[genes_of_interest,], 
                         color_by = "leiden_mEPC_CAF", 
                         ncol = 1, cell_size = 3) +
  scale_color_manual(values = color_mapping)  
dev.off() 
mEPC_cells <- gbm_cds[, grepl("mEPC", pData(gbm_cds)$leiden_mEPC_CAF)]
head(pData(mEPC_cells))
color_mapping <- c(
  "mEPC 01" = "#E377C2",
  "mEPC 02" = "#B5BD61",
  "mEPC 03" = "#17BEDF",
  "mEPC 04" = "#AEC7E8",
  "mEPC 05" = "#FFBB78",
  "mEPC 06" = "#98DF8A",
  "mEPC 07" = "#FF9896",
  "mEPC 08" = "#C5B0D5",
  "mEPC 09" = "#C49C94",
  "mEPC 10" = "#F7B6D2"
)
genes_of_interest <- row.names(subset(fData(mEPC_cells), gene_short_name %in% c("WDR54", "CD8A"))) ## 必须得两个，不然报错，拉一个凑数
if(length(genes_of_interest) == 0) {
    stop("没有找到相关基因。请检查基因名称是否正确。")
}
outputFile <- file.path(outputDir, "gene_expression_branched_pseudotime_EPC.png")
mEPC_cells$leiden_mEPC_CAF <- factor(mEPC_cells$leiden_mEPC_CAF, levels = names(color_mapping))
png(outputFile, width = 600, height = 600)
plot_genes_in_pseudotime(mEPC_cells[genes_of_interest,], 
                         color_by = "leiden_mEPC_CAF", 
                         ncol = 1, cell_size = 3) +
  scale_color_manual(values = color_mapping)  
dev.off()  
saveRDS(gbm_cds, "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/monocle/monocle2/WDR54/gbm_cds.rds")
### monocle3
library(monocle3)
library(Matrix)
library(dplyr)
install_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/Rpackage"
library(Seurat, lib.loc = install_path)
packageVersion("Seurat")
setwd("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/monocle/monocle3")
seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/vector/MNFHC_for_vector.rds")
Idents(seurat_obj) <- seurat_obj@meta.data$leiden_mEPC_CAF 
file_path <- '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2D/cytotrace/cytotrace_meta_data.csv'
cytotrace_data <- read.csv(file_path)
seurat_cell_ids <- colnames(seurat_obj)
csv_cell_ids <- cytotrace_data$CellID
cytotrace_data_sorted <- cytotrace_data[match(seurat_cell_ids, csv_cell_ids), ]
seurat_obj@meta.data$CytoTRACE2_Score <- as.numeric(cytotrace_data_sorted$CytoTRACE2_Score)
seurat_obj@meta.data$CytoTRACE2_Relative <- as.numeric(cytotrace_data_sorted$CytoTRACE2_Relative)
seurat_obj@meta.data$CytoTRACE2_Potency <- as.character(cytotrace_data_sorted$CytoTRACE2_Potency)
seurat_obj@meta.data$preKNN_CytoTRACE2_Potency <- as.character(cytotrace_data_sorted$preKNN_CytoTRACE2_Potency)
seurat_obj@meta.data$preKNN_CytoTRACE2_Score <- as.numeric(cytotrace_data_sorted$preKNN_CytoTRACE2_Score)
cellnames <- read.csv("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/vector/cellname.csv", header = TRUE, stringsAsFactors = FALSE)
cell_index_list <- cellnames$cell_index
seurat_obj <- subset(seurat_obj, cells = cell_index_list)
seurat_obj$celltype = seurat_obj@meta.data$leiden_mEPC_CAF
expression_matrix <- GetAssayData(seurat_obj, assay = 'RNA', slot = 'counts')
cell_metadata <- seurat_obj@meta.data 
gene_annotation <- data.frame(gene_short_name = rownames(seurat_obj))
rownames(gene_annotation) <- rownames(seurat_obj)
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)

output_path <- "pca_variance_explained.png"
png(output_path)
plot_pc_variance_explained(cds)
dev.off()
cds <- reduce_dimension(cds, 
                        core = 8, 
                        umap.fast_sgd = TRUE)

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat_obj, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
output_path <- "umap.png"
png(output_path, width = 1200, height = 900)
plot_cells(cds, 
           color_cells_by = "leiden_mEPC_CAF", 
           reduction_method = "UMAP",  
           cell_size = 2,  
           alpha = 0.8,  
           cell_stroke = I(3)  
)
dev.off()
marker_test_res <- top_markers(cds, group_cells_by="leiden_mEPC_CAF", 
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>% 
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

output_path <- "top_specific_marker_genes_by_group.png"
png(output_path)
plot_genes_by_group(cds, 
                    top_specific_marker_ids,
                    group_cells_by = "leiden_mEPC_CAF",
                    ordering_type = "maximal_on_diag",
                    max.size = 3)

dev.off()
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
output_path <- "cells_by_cell_type.png"
png(output_path, width = 1200, height = 900)

plot_cells(cds,
           color_cells_by = "CytoTRACE2_Relative",  
           label_groups_by_cluster = FALSE,  
           label_leaves = TRUE,  
           label_branch_points = TRUE,  
           label_roots = TRUE,  
           cell_size = 2,  
           alpha = 0.8,  
           cell_stroke = I(3),  
           show_trajectory_graph = TRUE,  
           trajectory_graph_color = "#1A1A1A",  
           trajectory_graph_segment_size = 4,  
           graph_label_size = 6  
)

dev.off()
get_earliest_principal_node <- function(cds) {
  time_threshold <- quantile(colData(cds)[, "CytoTRACE2_Relative"], 0.75)
  cell_ids <- which(colData(cds)[, "CytoTRACE2_Relative"] >= time_threshold)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))
  ]
  
  return(root_pr_nodes)
}
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))
output_path <- "cells_by_cell_type_ordered.png"
png(output_path, width = 1200, height = 900)
plot_cells(cds,
           color_cells_by = "CytoTRACE2_Relative",  
           label_groups_by_cluster = FALSE,  
           label_leaves = TRUE,  
           label_branch_points = TRUE,  
           label_roots = TRUE,  
           cell_size = 2,  
           alpha = 0.8,  
           cell_stroke = I(3),  
           show_trajectory_graph = TRUE,  
           trajectory_graph_color = "#1A1A1A",  
           trajectory_graph_segment_size = 4,  
           graph_label_size = 6  
)
dev.off()

### slingshot
library(slingshot)
library(tidyverse)
library(uwot)
library(viridisLite)
library(RColorBrewer)
library(mclust)
library(RColorBrewer)
library(grDevices)
library(SingleCellExperiment)
install_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/Rpackage"
library(Seurat, lib.loc = install_path)
packageVersion("Seurat")
setwd("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/Slingshot")
seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/vector/MNFHC_for_vector.rds")
print(seurat_obj)
sce <- as.SingleCellExperiment(seurat_obj)
library(SingleCellExperiment)
umap_coords <- read.csv("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/Slingshot/X_umap.csv", header = TRUE, stringsAsFactors = FALSE)
umap_coords_ordered <- umap_coords[match(colnames(sce), umap_coords[["CellID"]]), ]
colnames(umap_coords_ordered)[colnames(umap_coords_ordered) == "UMAP1"] <- "UMAP_0"
colnames(umap_coords_ordered)[colnames(umap_coords_ordered) == "UMAP2"] <- "UMAP_1"
reducedDims(sce)$UMAP <- as.matrix(umap_coords_ordered[, c("UMAP_0", "UMAP_1")])
head(reducedDims(sce)$UMAP, 5)
sce <- slingshot(sce, clusterLabels = 'leiden_mEPC_CAF', reducedDim = 'UMAP', approx_points = FALSE)
sce
summary(sce$slingPseudotime_1)
colData(sce)$slingshot
my_colors <- c(
  "CAF 01"  = "#1F77B4",
  "CAF 02"  = "#FF7E0E",
  "CAF 03"  = "#279A68",
  "CAF 04"  = "#D62728",
  "CAF 05"  = "#AA40FC",
  "CAF 06"  = "#8C564B",
  "mEPC 01" = "#E377C2",
  "mEPC 02" = "#B5BD61",
  "mEPC 03" = "#17BEDF",
  "mEPC 04" = "#AEC7E8",
  "mEPC 05" = "#FFBB78",
  "mEPC 06" = "#98DF8A",
  "mEPC 07" = "#FF9896",
  "mEPC 08" = "#C5B0D5",
  "mEPC 09" = "#C49C94",
  "mEPC 10" = "#F7B6D2"
)
png(filename = "custom_colors_plot_with_border.png", width = 2850, height = 2400, res = 300)
group_info <- colData(sce)$leiden_mEPC_CAF
plot_colors <- adjustcolor(my_colors[group_info], alpha.f = 0.6)
plot(
  reducedDims(sce)$UMAP,
  pch = 21,
  bg = plot_colors,
  col = plot_colors,
  cex = 0.5,
  asp = 8.6 / 9.5,
  main = "UMAP with Custom Color Mapping & Matching Borders"
)
lines(
  SlingshotDataSet(sce),
  lwd = 2,
  type = "lineages",
  col = "black"
)
dev.off()
rd <- reducedDim(sce, "UMAP")
cl <- colData(sce)$leiden_mEPC_CAF
lin1 <- getLineages(rd, cl)
png(filename = "custom_colors_plot.png", width = 2850, height = 2400, res = 300)
plot(rd,
     pch = 21,
     bg = plot_colors,
     col = plot_colors,
     cex = 0.5,
     asp = 10 / 9.5,
     main = "UMAP with Custom Color Mapping & Matching Borders"
)
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')
dev.off()
crv1 <- getCurves(lin1)
png(filename = "curve.png", width = 2850, height = 2400, res = 300)
plot(rd,
     pch = 21,
     bg = plot_colors,
     col = plot_colors,
     cex = 0.5,
     asp = 8.6 / 9.5,
     main = "UMAP with Custom Color Mapping"
)
line_colors <- brewer.pal(9, "Set1")
lines(SlingshotDataSet(crv1), lwd = 3, col = line_colors)
dev.off()
pt2 <- assay(crv1, "pseudotime")[, "Lineage2"]
cell_colors <- rep("lightgrey", length(pt2))
is_lineage2 <- !is.na(pt2)
plasma_n <- 100
plasma_pal <- plasma(plasma_n)
temp_rank <- rank(pt2[is_lineage2], ties.method = "first")
color_idx <- ceiling(plasma_n * temp_rank / max(temp_rank))
cell_colors[is_lineage2] <- plasma_pal[color_idx]
output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/Slingshot/curves_rainbow2.png"
png(filename = "curve111.png", width = 2850, height = 2400, res = 300)
plot(rd,
     pch = 16,
     bg = plot_colors,
     col = cell_colors,
     cex = 0.5,
     asp = 8.6 / 9.5,
     main = "UMAP with Custom Color Mapping"
)
line_colors <- brewer.pal(9, "Set1")
lines(SlingshotDataSet(crv1), lwd = 5, col = line_colors)
dev.off()

### vector
packages <- c("circlize", "gatepoints", "stringr", "igraph", "gmodels")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
setwd("C:/R/C_Project/vector")
unzip_dir <- "C:/R/C_Project/vector/Vector-master"
source(file.path(unzip_dir, "Vector.R"))
seurat_obj <- readRDS("C:/R/C_Project/vector/MNFHC_for_vector.rds")
cellnames <- read.csv("C:/R/C_Project/vector/cellname.csv", header = TRUE, stringsAsFactors = FALSE)
cell_index_list <- cellnames$cell_index
seurat_obj <- subset(seurat_obj, cells = cell_index_list)
VEC <- seurat_obj@reductions$umap@cell.embeddings
rownames(VEC) <- colnames(seurat_obj)
PCA <- seurat_obj@reductions$pca@cell.embeddings
PCA <- vector.rankPCA(PCA)
OUT <- vector.buildGrid(VEC, N = 30, SHOW = TRUE)
OUT <- vector.buildNet(OUT, CUT = 1, SHOW = TRUE)
OUT <- vector.getValue(OUT, PCA, SHOW = TRUE)
OUT <- vector.gridValue(OUT, SHOW = TRUE)
OUT <- vector.autoCenter(OUT, UP = 0.90, SHOW = TRUE)
OUT <- vector.drawArrow(OUT, P = 0.10, SHOW = TRUE, COL = OUT$COL, SHOW.SUMMIT = FALSE, AL = 100)
OUT <- vector.buildGrid(VEC, N = 30, SHOW = TRUE)
OUT <- vector.buildNet(OUT, CUT = 1, SHOW = TRUE)
OUT <- vector.getValue(OUT, PCA, SHOW = TRUE)
OUT <- vector.gridValue(OUT, SHOW = TRUE)
OUT <- vector.selectCenter(OUT)
OUT <- vector.drawArrow(OUT, P = 0.10, SHOW = TRUE, COL = OUT$COL, SHOW.SUMMIT = FALSE, AL = 100)
source_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2D/vector/Rplot06.png"
destination_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/figures"
destination_file <- file.path(destination_dir, "F2D_1_cytotrace.PNG")
if (!dir.exists(destination_dir)) {
  dir.create(destination_dir, recursive = TRUE)
}
file.copy(source_file, destination_file)
source_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2D/vector/Rplot08.png"
destination_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/figures"
destination_file <- file.path(destination_dir, "F2D_2_cytotrace.PNG")
if (!dir.exists(destination_dir)) {
  dir.create(destination_dir, recursive = TRUE)
}
file.copy(source_file, destination_file)