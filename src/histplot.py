import sys
import matplotlib.pyplot as plt
import numpy as np

c = np.array([])
k = sys.argv[1]
l = eval(sys.argv[2])
r = eval(sys.argv[3])

f = open(k + "mer_kmc_hist.txt", "r")
for x in f:
  c = np.append(c, eval(x))

plt.plot(c, color='r', label='solid')
plt.plot(l, c[l], '-o', color='orange')
plt.plot(r, c[r], '-o', color='orange')

plt.legend()
plt.xlim([2, 150])
plt.ylim([0, 1.5*1e6])
plt.show()
plt.savefig(k + "mer.png")

# delta = np.array([])
# delta = (c[2:-1] - c[3:]) / c[2:-1]
# plt.xlim([0, 100])
# plt.plot(delta)
# plt.plot(l, delta[l], "-o")
# plt.savefig(k + "_delta.png")


