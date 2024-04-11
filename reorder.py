import numpy as np

data1 = np.genfromtxt("example/save_data.txt",
                      usecols=(1, 2, 3), unpack=True).T
data2 = np.genfromtxt("build/test", usecols=(1, 2, 3), unpack=True).T

col1 = np.genfromtxt("example/save_data.txt",
                     usecols=(0), unpack=True)
col2 = np.genfromtxt("build/test", usecols=(0), unpack=True)


d1dict = {}
d2dict = {}
for i, data in enumerate(col1):
    d1dict[int(col1[i])] = data1[i]

for i, data in enumerate(col2):
    d2dict[int(col2[i])] = data2[i]

sum_field1 = np.zeros(3)
sum_field2 = np.zeros(3)

for key in d1dict:
    sum_field1 += d1dict[key]
    sum_field2 += d2dict[key]
angperau = 0.52917721092
print(sum_field1*angperau**2)
print(sum_field2*angperau**2)
