import numpy as np
import os
from pathlib import Path
print(os.getcwd())
datas = np.load('.\datas\data.npz')

Z = datas['Z']
W = datas['W']

Z0 = Z[0]
Z1 = Z[1]
Z2 = Z[2]
Z3 = Z[3]
Z4 = Z[4]
Z5 = Z[5]
Z6 = Z[6]
Z7 = Z[7]
Z8 = Z[8]
Z9 = Z[9]
W0 = W[0]
W1 = W[1]
W2 = W[2]
W3 = W[3]
W4 = W[4]
W5 = W[5]
W6 = W[6]
W7 = W[7]
W8 = W[8]
W9 = W[9]

np.savez('.\datas\distribution_n{}_y{}.npz'.format(1, None), Z=Z0, W=W0)

np.savez('.\datas\distribution_n{}_y{}.npz'.format(2, None), Z=Z1, W=W1)

np.savez('.\datas\distribution_n{}_y{}.npz'.format(3, None), Z=Z2, W=W2)

np.savez('.\datas\distribution_n{}_y{}.npz'.format(4, None), Z=Z3, W=W3)

np.savez('.\datas\distribution_n{}_y{}.npz'.format(5, None), Z=Z4, W=W4)

np.savez('.\datas\distribution_n{}_y{}.npz'.format(6, None), Z=Z5, W=W5)

np.savez('.\datas\distribution_n{}_y{}.npz'.format(7, None), Z=Z6, W=W6)

np.savez('.\datas\distribution_n{}_y{}.npz'.format(8, None), Z=Z7, W=W7)

np.savez('.\datas\distribution_n{}_y{}.npz'.format(9, None), Z=Z8, W=W8)

np.savez('.\datas\distribution_n{}_y{}.npz'.format(10, None), Z=Z9, W=W9)
