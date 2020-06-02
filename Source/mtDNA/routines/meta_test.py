from scipy.stats import norm
import math

pval = [0.084503005, 0.00599414, 0.017554205, 0.261152011]
ni = [656, 729, 1803, 2670]
zi = [norm.ppf(x / 2.0) for x in pval]
wi = [math.sqrt(x) for x in ni]
z_num = 0
z_den = 0
for i in range(0, len(pval)):
    z_num += zi[i] * wi[i]
    z_den += wi[i] * wi[i]

z = z_num / math.sqrt(z_den)
p = 2 * norm.cdf(z)
print(p)
