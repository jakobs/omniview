import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import math

class OmniCam:
    def __init__( self, k, c, f, mirror_diameter ):
        self.k = k
        self.c = c
        self.f = f
        self.mirror_diameter = mirror_diameter
        
        # derived values 
        self.a = c/2*np.sqrt((k-2)/k)
        self.b = c/2*np.sqrt(2/k)

    def y( self, x ):
        return -np.sqrt((1 + (x**2) / (self.b**2))*(self.a**2)) + self.c/2

    def p_u( self, x, y ):
        tana = y/x
        x1 = (self.b**2*self.c*tana- \
                self.a*self.b*np.sqrt(4*self.b**2*tana**2+self.c**2-4*self.a**2))/ \
                (2*(self.b**2*tana**2-self.a**2))
        y1 = x1 * tana 
        def ev(x2):
            y2 = self.y(x2)
            theta = 2.0*(np.arctan((y2-y)/(x-x2))+np.arctan(((self.c/2-y2)*self.b**2)/(self.a**2*x2)))
            l1sq = x2**2+(self.c+self.f-y2)**2
            l2sq = (x-x2)**2+(y-y2)**2
            l3sq = x**2 + (self.c+self.f-y)**2
            return l1sq + l2sq - l3sq - 2.0*np.sqrt(l1sq)*np.sqrt(l2sq)*np.cos(theta) 

#        def find_x2:
#            i = 1.0
#            while 
#
#        xd = np.linspace( x1-5.0, x1+5.0 )
#        plt.plot( xd, ev(xd) )
#        plt.show()
        intv = 0.0
        while np.sign( ev(x1-intv) ) == np.sign( ev(x1+intv) ):
            intv += 1.0

        x2 = scipy.optimize.brentq( ev, x1 - intv, x1 + intv )
        y2 = self.y( x2 )
        v = x1*self.f*(self.c+self.f-y2) / (x2*(self.c-y1))
        u = self.f*v/(v-self.f)
        xu = x2*u/(self.c+self.f-y2)
        yu = self.c+self.f-u

        return xu, yu

#        plt.plot([x,x1],[y,y1], 'k-')
#        plt.plot([x,x2],[y,y2], 'k-')
#        plt.plot([x1,0],[y1,self.c], 'k-')
#        plt.plot([x2,0],[y2,self.c+self.f], 'k-')
#        
#        plt.plot([xu],[yu], 'ro')

    def plot_mirror( self, plt = plt ):
        x1 = np.linspace( -self.mirror_diameter/2, self.mirror_diameter/2 )
        line, = plt.plot(x1, self.y(x1), '-', linewidth=2)
        # draw helper axis
        line, = plt.plot([-self.mirror_diameter/2, self.mirror_diameter/2], [0,0], '--', linewidth=1)

    def virtual_points( self, x, y ):
        xu = []
        yu = []
        for xi, yi in zip(x,y):
            xui, yui = self.p_u( xi, yi )
            xu.append(xui)
            yu.append(yui)
        return xu, yu

    def plot_circle( self ):
        alpha = np.linspace( 1.0, 3.0 )
        x, y = self.virtual_points( np.sin(alpha)*100.0, np.cos(alpha)*100.0 )
        plt.plot(x,y, 'ro')

oc = OmniCam( k=10.0, c=20.0, f=5.0, mirror_diameter=30.0)
oc.plot_mirror()
oc.plot_circle()

plt.show()

