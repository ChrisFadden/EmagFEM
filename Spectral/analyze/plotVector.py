import struct
import numpy as np
import matplotlib.pyplot as plt

f = open('../build/test.dat','rb').read()

#7 comes from 8-byte longs used for size of vectors
vecSize = int.from_bytes(f[0:7],byteorder='little')

byteSize = int.from_bytes(f[8:15],byteorder='little')

strForm = 'f'
if(byteSize == 8):
    strForm = 'd'

x = np.asarray(struct.unpack(str(vecSize) + strForm,f[16:]))

plt.plot(x)
plt.show()
