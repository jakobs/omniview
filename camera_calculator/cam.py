import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import math

def annotate_dim(ax,xyfrom,xyto,text=None):
    if text is None:
        text = str(np.sqrt( (xyfrom[0]-xyto[0])**2 + (xyfrom[1]-xyto[1])**2 ))
    ax.annotate("",xyfrom,xyto,arrowprops=dict(arrowstyle='<->'))
    ax.text((xyto[0]+xyfrom[0])/2,(xyto[1]+xyfrom[1])/2,text,fontsize=8)

class OmniCam:
    def __init__( self, k, c, f, mirror_diameter, fnumber, pixel_spacing = 0.0055, pixels = 2048):
        self.k = k
        self.c = c
        self.f = f
        self.mirror_diameter = mirror_diameter

        # derived values 
        self.a = c/2.0*np.sqrt((k-2.0)/k)
        self.b = c/2.0*np.sqrt(2.0/k)

        # max angle
        xe = mirror_diameter/2.0 
        ye = self.y( xe )
        self.max_angle = math.atan2( xe, ye )

        # calc volume
        yc = self.y( 0 )
        yi = np.arange( ye, yc, 0.1 )
        self.volume = np.trapz( self.x( yi )**2 * math.pi, yi )
        self.mirror_height = yc - ye

        # sensor 
        self.pixel_spacing = pixel_spacing
        self.pixels = pixels
        self.sensor_width = pixel_spacing * pixels
        self.edge_view = ( c - ye ) * self.sensor_width / f
        self.corner_view = self.edge_view * np.sqrt(2.0) 

        # DOF
        circle_of_confusion = pixel_spacing
        subject_distance = f+c
        self.fnumber = fnumber
        self.dof = 2*fnumber*circle_of_confusion*f**2*(subject_distance)**2/ \
            (f**4-fnumber**2*circle_of_confusion**2*subject_distance**2)

    def print_stats(self):
        print "max angle: %d deg (%d over horizon)" % (self.max_angle / math.pi * 180, self.max_angle / math.pi * 180 - 90)
        print "volume: %d ccm" % ( self.volume * 1e-3 ) 
        print "edge view (corner view): %d mm (%d mm)" % ( self.edge_view, self.corner_view )
        print "mirror diameter: %d mm" % ( self.mirror_diameter )
        print "mirror height: %d mm" % ( self.mirror_height )
        print "DOF: %d mm @ f/%.1f" % ( self.dof, self.fnumber )
        print "mirror parameters:"
        print "curve k=%.1f, focal point distance c=%.1f" % (self.k, self.c)
        print "a=%f b=%f" % (self.a, self.b)

    def y( self, x ):
        return -np.sqrt((1.0 + (x**2) / (self.b**2))*(self.a**2)) + self.c/2.0

    def x( self, y ):
        return np.sqrt( ((y-self.c/2.0)**2/self.a**2-1.0)*self.b**2 )
    
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
        y = self.y(self.mirror_diameter/2)
        line, = plt.plot([-self.mirror_diameter/2, self.mirror_diameter/2], [y,y], '-', linewidth=1)

        text = "$y=-a\\sqrt{\\frac{x^2}{b^2}+1}$\n$a=%.5f$\n$b=%.5f$" % (self.a, self.b)
        plt.annotate( text, xy=(-self.mirror_diameter/2, y), xytext=(-self.mirror_diameter/2 - 100, y), \
                arrowprops=dict(facecolor='black', shrink=0.01) )

        annotate_dim( plt.gca(), [-self.mirror_diameter/2, y - 15], [self.mirror_diameter/2, y - 15] )
        annotate_dim( plt.gca(), [self.mirror_diameter/2 + 15, y], [self.mirror_diameter/2 + 15, 0] )
        annotate_dim( plt.gca(), [self.mirror_diameter/2 + 15, 0], [self.mirror_diameter/2 + 15, self.y(0)] )
        annotate_dim( plt.gca(), [-(self.mirror_diameter/2 + 15), 0], [-(self.mirror_diameter/2 + 15), self.c] )
        
        # draw helper axis
        line, = plt.plot([-self.mirror_diameter/2, self.mirror_diameter/2], [0,0], '--', linewidth=1)
        line, = plt.plot([0,0], [0,self.c+self.f], '--', linewidth=1)

        # draw max lines
        plt.plot( [0,self.mirror_diameter*2], [0, y*4] )
        plt.plot( [0,-self.mirror_diameter*2], [0, y*4] )

    def virtual_points( self, x, y ):
        xu = []
        yu = []
        for xi, yi in zip(x,y):
            xui, yui = self.p_u( xi, yi )
            xu.append(xui)
            yu.append(yui)
        return xu, yu

    def plot_circle( self ):
        alpha = np.linspace( 0.4, self.max_angle )
        for i in range(1,5):
            px, py = np.sin(alpha)*100.0*i, np.cos(alpha)*100.0*i 
#            plt.plot( px, py, 'go' )
            x, y = self.virtual_points( px, py )
            plt.plot(x,y, 'ro')
#        x, y = self.virtual_points( np.sin(alpha)*1000.0, -np.cos(alpha)*1000.0 )
#        plt.plot(x,y, 'ro')

# Lensagon C5M1618GS
# 16mm F1.8, 1.1" Sensor, min obj distance 97.7mm
#oc = OmniCam( k=20.0, c=140.0, f=16.0, mirror_diameter=58.0, fnumber=32.0)

# Lensagon C5M1618GS 
#oc = OmniCam( k=15.0, c=33.0, f=16.0, mirror_diameter=58.0, fnumber=16.0)

# Lensagon C5M3514GS 
# 35mm F1.4, 1.1" Sensor, min obj distance 110mm
print "Lensagon C5M3514GS" 
oc = OmniCam( k=17.0, c=155.0, f=35.0, mirror_diameter=58.0, fnumber=16.0)

# Lensagon C5M2514GS
#print "Lensagon C5M2514GS"
#oc = OmniCam( k=15.0, c=90.0, f=25.0, mirror_diameter=58.0, fnumber=16.0)

oc.print_stats()

oc.plot_mirror()
#oc.plot_circle()

plt.axes().set_aspect('equal')
#plt.add_subplot(111,aspect='equal')
plt.show()

