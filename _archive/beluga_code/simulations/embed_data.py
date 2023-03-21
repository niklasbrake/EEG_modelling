import sys
import numpy as np
import umap
from scipy.sparse import coo_array


def load_spike_correlations(filename):
    # Get pairwise correlation of spike trains from file
    C = np.loadtxt(filename,delimiter=",")
    idcs2,idcs = np.unique(C[:,:2],return_inverse=True)
    N = len(idcs2)
    newCoords = idcs.reshape([-1,2])
    preIdcs = idcs2.reshape([-1,1])
    dist_mat = coo_array((1-C[:,2], (newCoords[:,0], newCoords[:,1])), shape=(N,N))
    dist_mat = dist_mat + dist_mat.T
    return dist_mat,preIdcs

def sphere_embed(d):
    # Run UMAP on functional correlation matrix and plot
    sphere_mapper = umap.UMAP(output_metric='haversine', random_state=42, low_memory=True).fit(d)
    x = np.sin(sphere_mapper.embedding_[:, 0]) * np.cos(sphere_mapper.embedding_[:, 1])
    y = np.sin(sphere_mapper.embedding_[:, 0]) * np.sin(sphere_mapper.embedding_[:, 1])
    z = np.cos(sphere_mapper.embedding_[:, 0])
    x = np.arctan2(x, y)
    y = -np.arccos(z)
    theta = np.reshape(x,[-1,1])
    phi = np.reshape(y,[-1,1])
    return theta, phi

if __name__ == "__main__":
    pars = sys.argv
    file = pars[1]
    file2 = pars[2]
    d,preIdcs = load_spike_correlations(file)
    theta,phi = sphere_embed(d)
    X = np.concatenate((preIdcs,theta,phi),1)
    fmt = '%d', '%.2f', '%.2f'
    np.savetxt(file2,X,delimiter=",",fmt=fmt)
