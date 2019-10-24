from scipy.special import binom

k_mt_min = 1
k_nuc_min = 1
k_mt_max = 1
k_nuc_max = 4

mt_num = 0
nuc_num = 0

for k_mt in range(k_mt_min, k_mt_max + 1):
    mt_num += binom(k_mt_max, k_mt)
for k_nuc in range(k_nuc_min, k_nuc_max + 1):
    nuc_num += binom(k_nuc_max, k_nuc)

num_combinations = mt_num * nuc_num
print(str(int(num_combinations)))
