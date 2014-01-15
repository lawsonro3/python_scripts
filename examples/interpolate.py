from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 10, 10)
y = np.exp(-x/3.0)
f = interp1d(x, y)
f2 = interp1d(x, y, kind='cubic')
print f
print type(f)
xnew = np.linspace(0, 10, 40)
print xnew
plt.plot(x,y,'o',xnew,f(xnew),'-x', xnew, f2(xnew),'--o')
plt.legend(['data', 'linear', 'cubic'], loc='best')
plt.show()
raw_input("Press ENTER to exit") 

