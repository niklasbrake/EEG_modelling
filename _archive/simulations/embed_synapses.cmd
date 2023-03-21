set folder1=%1\presynaptic_network\correlations.csv
set folder1=%folder1:\=/%
set folder2=%1\presynaptic_network\UMAP_embedding.csv
set folder2=%folder2:\=/%

python ./functions/embed_data.py %folder1% %folder2%