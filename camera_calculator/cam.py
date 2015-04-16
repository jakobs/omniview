import numpy as np
import matplotlib.pyplot as plt


# draw mirror
mirror_diameter = 30
x1 = np.linspace( -mirror_diameter/2, mirror_diameter/2 )
k = 10.0
c = 20.0
a = c/2*np.sqrt((k-2)/k)
b = c/2*np.sqrt(2/k)
y1 = -np.sqrt((1 + (x1**2) / (b**2))*(a**2)) + c/2

print 'a = ' + str(a)
print 'b = ' + str(b)
line, = plt.plot(x1, y1, '-', linewidth=2)

# draw helper axis
line, = plt.plot([-mirror_diameter/2, mirror_diameter/2], [0,0], '--', linewidth=1)

plt.show()
