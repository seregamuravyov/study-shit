import math
import numpy
import sys
import pylab
from matplotlib import mlab
n = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
d = [2, 2, 3, 4, 4, 4, 6, 6, 6, 6, 7, 8]
result = []

j = 0
def fastPow (t, k): 
    res = 1
    while k:
        if (k & 1):
            res *= t
        k = k >> 1
        if k == 0:
            break
        t *= t
 
    return res

def bernulli(n, k, p0):
    koef = math.factorial(n) / ((math.factorial(k))*(math.factorial(n - k)))
    return koef * math.pow(p0, k) * math.pow(1 - p0, (n - k))

while j < len(n):
   i = math.floor((d[j] - 1)/2)
   s = i + 1
   summ = 0
   while s < n[j]:
       summ += bernulli(n[j], s, p0 = 0.01)
       s += 1
   result.append(summ)    
   j += 1
   
print (result)

xlist = mlab.frange (2, 24, 2)
ylist = result
pylab.plot (xlist, ylist)
pylab.show()
