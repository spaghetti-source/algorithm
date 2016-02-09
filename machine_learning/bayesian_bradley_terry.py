#
# Bayesian version of Bradley-Terry model
# 
# Reference: 
#   Ruby C. Weng and Chih-Jen Lin (2011):
#   A Bayesian approximation method for online ranking.
#   Jornal on Machine Learning Research, vol.12, pp.267--300
#   (no-team version of Algorithm 1)
#

import math
import random
from collections import defaultdict
import matplotlib.pyplot as plt

beta = 25.0/6.0
kappa = 0.0001
mu = defaultdict(lambda: 25.0)
sigma = defaultdict(lambda: 25.0/3.0)

# exp(a) / (exp(a) + exp(b))
def logit(a, b):
    return 1.0 / (1.0 + math.exp(b - a))
def c(i,q):
    return (sigma[i]**2 + sigma[q]**2 + 2.0 * beta**2)**0.5
def p(i,q):
    return logit(mu[i]/c(i,q), mu[q]/c(i,q))
def gamma(i,q):
    return sigma[i] / c(i,q)

# result: {player} -> Int (smaller is better) 
def update(result):
    for i in result:
        Omega = 0
        Delta = 0
        for q in result:
            if i == q: continue
            if result[i] <  result[q]: s = 1.0
            if result[i] >  result[q]: s = 0.0
            if result[i] == result[q]: s = 0.5
            Omega += sigma[i]**2 / c(i,q) * (s - p(i,q))
            Delta += gamma(i,q) * sigma[i]**2 / c(i,q)**2 * p(i,q) * p(q,i)
        mu[i]     += Omega
        sigma[i]  *= max(1.0 - Delta, kappa)

# verification
out = []
for iter in range(10000):
    if random.random() < 0.25:
        update({'a': 1, 'b': 2})
    else:
        update({'a': 2, 'b': 1})
    out.append(p('a', 'b'))
plt.plot(out)
plt.show()
