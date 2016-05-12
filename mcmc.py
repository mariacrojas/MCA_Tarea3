from matplotlib import rcParams
from scipy.stats import norm
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
import matplotlib.mlab as mlab

rcParams['savefig.dpi'] = 100

data=np.loadtxt("salida1.dat")
x=data[:,0]
y=data[:,1]
def straight_line_log_likelihood(x, y, sigmay, a, b,c,d):
    '''
    Returns the log-likelihood of drawing data values *y* at
    known values *x* given Gaussian measurement noise with standard
    deviation with known *sigmay*, where the "true" y values are
    *y_t = m * x + b*

    x: list of x coordinates
    y: list of y coordinates
    sigmay: list of y uncertainties
    m: scalar slope
    b: scalar line intercept

    Returns: scalar log likelihood
    '''
#    funcion=np.log10((10**a)/(((10**x/10**b)**c)((1+10**x/10**b)**d)))
    funcion=np.log10((10**a)/(((10**x/10**b)**c)*((1+10**x/10**b)**d)))
    return (np.sum(np.log(1./(np.sqrt(2.*np.pi) * sigmay))) +
            np.sum(-0.5 * (y - funcion)**2 / sigmay**2))
    
def straight_line_log_prior(m, b):
    return 0.
    
def straight_line_log_posterior(x,y,sigmay, m,b):
    return (straight_line_log_likelihood(x,y,sigmay, a,b,c,d) +
            straight_line_log_prior(m, b))

def straight_line_posterior(x, y, sigmay, m, b):
    return np.exp(straight_line_log_posterior(x, y, sigmay, a, b,c,d))

# initial m, b
a,b,c,d = 1.0, 0.0,0.0,3.0

# step sizes
astep, bstep,cstep, dstep = 1.5,0.5,1,0.1
        
# how many steps?
nsteps = 1000
    
chain = []
probs = []
naccept = 0
    
print 'Running MH for', nsteps, 'steps'
sigmay=1
# First point:
L_old    = straight_line_log_likelihood(x, y, sigmay, a, b,c,d)
p_old    = straight_line_log_prior(a, b)
#p_old    = 0
prob_old = np.exp(L_old + p_old)

for i in range(nsteps):
    # step
    anew = a + np.random.normal() * astep
    bnew = b + np.random.normal() * bstep
    cnew = c + np.random.normal() * cstep
    dnew = d + np.random.normal() * dstep
    # evaluate probabilities
    # prob_new = straight_line_posterior(x, y, sigmay, mnew, bnew)

    L_new    = straight_line_log_likelihood(x, y, sigmay, anew, bnew,cnew, dnew)
    p_new    = straight_line_log_prior(anew, bnew)
#    p_new    = 0
    prob_new = np.exp(L_new + p_new)

    if (prob_new / prob_old > np.random.uniform()):
        # accept
        a = anew
        b = bnew
        c = cnew
        d = dnew
        L_old = L_new
        p_old = p_new
        prob_old = prob_new
        naccept += 1
    else:
        # Stay where we are; m,b stay the same, and we append them
        # to the chain below.
        pass

    chain.append((a,b,c,d))
    probs.append((L_old,p_old))
aa = [a for a,b,c,d in chain]
bb = [b for a,b,c,d in chain]
cc = [c for a,b,c,d in chain]
dd = [d for a,b,c,d in chain]

plt.figure(1)
plt.plot(aa)

plt.figure(2)
plt.plot(bb)

plt.figure(3)
plt.plot(cc)

plt.figure(4)
plt.plot(dd)
plt.title('beta')
plt.figure(5)
na,binsa,patchesa=plt.hist(aa,10,normed=True)
(mua,sigmaa)=norm.fit(aa)
ya=norm.pdf(binsa,mua,sigmaa)
plt.plot(binsa,ya,'r--')
plt.title('rho_0')
plt.figure(6)
nb,binsb,patchesb=plt.hist(bb,10,normed=True)
(mub,sigmab)=norm.fit(bb)
yb=norm.pdf(binsb,mub,sigmab)
plt.plot(binsb,yb,'r--')
plt.title('r_c')
plt.figure(7)
nc,binsc,patchesc=plt.hist(cc,bins=10,normed=True)
(muc,sigmac)=norm.fit(cc)
yc=norm.pdf(binsc,muc,sigmac)
plt.plot(binsc,yc,'r--')
plt.title('alpha')
plt.figure(8)
nd,binsd,patchesd=plt.hist(dd,10,normed=True)
(mud,sigmad)=norm.fit(dd)
yd=norm.pdf(binsd,mud,sigmad)
plt.plot(binsd,yd,'r--')
plt.title('beta')
plt.figure(0)
plt.scatter(x,y)
x1=np.linspace(-1.5,2.5,100)
y1=np.log10(10**mua/(((10**x1/10**(mub))**(muc))*((1+10**x1/10**(mub))**(mud))))
plt.scatter(x1,y1)

plt.show()

print 'Incertidumbres'
print 'rho_0'
print sigmaa
print 'r_c'
print sigmab
print 'alpha'
print sigmac
print 'beta'
print sigmad


