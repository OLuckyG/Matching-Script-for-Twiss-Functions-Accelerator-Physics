##  Optimizer for matching Twiss functions and magnetic elements:
##  This optimizer uses minize function to match twiss functions(beta,alhpa) functions to a lattice. 
##  It uses the linear transformation (3x3) matrix that transforms twiss functions in the linear elements
##  The elements that are included are:
##      - Quadrupoles
##      - Dipoles 
##      - Solenoids
##      - Drifts
##  This script uses the transfer matrix of elements and can sepcify the conditions you want to match.
##  For instance: zero dispersions, matched beta functions, upper and lower bound on beta functions, specific phase advances,
## Author: OluckyG
## Date: 07/24/2024
import numpy as np
import matplotlib.pyplot as plt
import random as random
import time
from numpy import sqrt, pi, cos, sin, arctan, sinh, cosh
from scipy.optimize import minimize

start_time = time.time()  ## For recording how long it took!

## Drift Class:
class Drift:
    def __init__(self,l):
        self.l = l
    
    def t_matrix(self):
        M=np.array([[1,self.l,0,0],[0,1,0,0],[0,0,1,self.l],[0,0,0,1]])
        return M
    
    def lreturn(self):
        return self.l

## Quadrupole Class:
class Quad:
    def __init__(self,l,k):
        self.l = l
        self.k = k
    
    def t_matrix(self):
        if (self.k>0):
            sqrtk = sqrt(self.k)
            M=np.array([[cos(sqrtk*self.l),(1/sqrtk)*sin(sqrtk*self.l),0,0],[-sqrtk*sin(sqrtk*self.l),cos(sqrtk*self.l),0,0],[0,0,cosh(sqrtk*self.l),(1/sqrtk)*sinh(sqrtk*self.l)],[0,0,sqrtk*sinh(sqrtk*self.l),cosh(sqrtk*self.l)]])
        else:
            sqrtk = sqrt(np.abs(self.k))
            M=np.array([[cosh(sqrtk*self.l),(1/sqrtk)*sinh(sqrtk*self.l),0,0],[sqrtk*sinh(sqrtk*self.l),cosh(sqrtk*self.l),0,0],[0,0,cos(sqrtk*self.l),(1/sqrtk)*sin(sqrtk*self.l)],[0,0,-sqrtk*sin(sqrtk*self.l),cos(sqrtk*self.l)]])
        
        return M
    
    def lreturn(self):
        return self.l

## Dipole Class
class Dipole:
    def __init__(self,l,rho):
        self.l = l
        self.rho = rho
    
    def t_matrix(self):
        M = np.array([[np.cos(self.l/self.rho),self.rho*np.sin(self.l/self.rho),0.0,0.0],[-(1/self.rho)*sin(self.l/self.rho),np.cos(self.l/self.rho),0.0,0.0],[0.0,0.0,1.0,self.l],[0.0,0.0,0.0,1.0]])
        return M
    
    def lreturn(self):
        return self.l
## Solenoid Class: Assuming Larmor plane(decoupled in both planes)
class solenoid:
    def __init__(self,l,K):
        self.l = l
        self.K = K
    
    def t_matrix(self):
        theta = 0.5*self.K*self.l
        # M=np.array([[(cos(theta)**2),(2/self.K)*sin(theta)*cos(theta),sin(theta)*cos(theta),(2/self.K)*(sin(theta)**2)],[-(self.K/2)*sin(theta)*cos(theta),(cos(theta)**2),-(self.K/2)*(sin(theta)**2),sin(theta)*cos(theta)],[-sin(theta)*cos(theta),-(2/self.K)*(sin(theta)**2),(cos(theta)**2),(2/self.K)*sin(theta)*cos(theta)],[(self.K/2)*(sin(theta)**2),-sin(theta)*cos(theta),-(self.K/2)*sin(theta)*cos(theta),(cos(theta)**2)]])
        M=np.array([[cos(theta),(2/self.K)*sin(theta),0.0,0.0],[-(self.K/2)*sin(theta),cos(theta),0.0,0.0],[0.0,0.0,cos(theta),(2/self.K)*sin(theta)],[0.0,0.0,-(self.K/2)*sin(theta),cos(theta)]])
        return M
    
    def lreturn(self):
        return self.l

## Solenoid Class, with rotation:
class solenoidfull:
    def __init__(self,l,k):
        self.l = l
        self.k = k
    
    def t_matrix(self):
        c = cos(self.l * self.k)
        s = sin(self.l * self.k)
        M = np.array([[(1+c)/2 , s/self.k, s/2, (1-c)/self.k],[-self.k*s*0.25,0.5*(1+c),-0.25*self.k*(1-c),s/2],[-s/2,-(1/self.k)*(1-c),0.5*(1+c),s/self.k],[(self.k/4)*(1-c),-s/2,-(self.k/4)*s,0.5*(1+c)]])
        return M
    def lreturn(self):
        return self.l

## Callback function to print the parameters that are being tried onto screen
def callback(xk):
    print(f"Current parameters: {xk}")

## Courant-Snyder 3x3 transfer matrix for linear optic functions for both x and y
def CStransfermat(betax0,betay0,alfx0,alfy0,gamx0,gamy0,Mmatrix):
    betaxf = (Mmatrix[0][0]**2)*betax0 -2*(Mmatrix[0][0]*Mmatrix[0][1])*alfx0 + (Mmatrix[0][1]**2)*gamx0
    alfaxf = (-Mmatrix[0][0]*Mmatrix[1][0])*betax0 + (Mmatrix[0][0]*Mmatrix[1][1] + Mmatrix[0][1]*Mmatrix[1][0])*alfx0 - (Mmatrix[0][1]*Mmatrix[1][1])*gamx0
    gammaxf = (Mmatrix[1][0]**2)*betax0 -2*(Mmatrix[1][1]*Mmatrix[1][0])*alfx0 + (Mmatrix[1][1]**2)*gamx0
    betayf = (Mmatrix[2][2]**2)*betay0 -2*(Mmatrix[2][2]*Mmatrix[2][3])*alfy0 + (Mmatrix[2][3]**2)*gamy0
    alfayf = (-Mmatrix[2][2]*Mmatrix[3][2])*betay0 + (Mmatrix[2][2]*Mmatrix[3][3] + Mmatrix[2][3]*Mmatrix[3][2])*alfy0 - (Mmatrix[2][3]*Mmatrix[3][3])*gamy0
    gammayf = (Mmatrix[3][2]**2)*betay0 -2*(Mmatrix[3][3]*Mmatrix[3][2])*alfy0 + (Mmatrix[3][3]**2)*gamy0
    return (betaxf,alfaxf,gammaxf,betayf,alfayf,gammayf)

## Dispersion function at dipoles, included dispersion in y as well
def dispersiontransportdipole(dx0,dpx0,dy0,dpy0,Mmatrix,rhodipole,ldipole):
    dxf = Mmatrix[0][0]*dx0 + Mmatrix[0][1]*dpx0 + np.sqrt(2)*rhodipole*(1-np.cos(ldipole/(np.sqrt(2)*rhodipole)))
    dpxf = Mmatrix[1][0]*dx0 + Mmatrix[1][1]*dpx0 + np.sin(ldipole/(np.sqrt(2)*rhodipole))
    dyf = Mmatrix[2][2]*dy0 + Mmatrix[2][3]*dpy0
    dpyf = Mmatrix[3][2]*dy0 + Mmatrix[3][3]*dpy0
    return (dxf,dpxf,dyf,dpyf)
## Dispersion at normal elements
def dispersiontransport(dx0,dpx0,dy0,dpy0,Mmatrix):
    dxf = Mmatrix[0][0]*dx0 + Mmatrix[0][1]*dpx0 + Mmatrix[0][2]*dy0 + Mmatrix[0][3]*dpy0
    dpxf = Mmatrix[1][0]*dx0 + Mmatrix[1][1]*dpx0 + Mmatrix[1][2]*dy0 + Mmatrix[1][3]*dpy0
    dyf = Mmatrix[2][0]*dx0 + Mmatrix[2][1]*dpx0 + Mmatrix[2][2]*dy0 + Mmatrix[2][3]*dpy0
    dpyf = Mmatrix[3][0]*dx0 + Mmatrix[3][1]*dpx0 + Mmatrix[3][2]*dy0 + Mmatrix[3][3]*dpy0
    return (dxf,dpxf,dyf,dpyf)

## Phase advance transportation through integration
## Need to use big step size for this, since it is kind of integrattion.
def mutransport(muI0,muII0,betx0,betay0,l):
    muIf = muI0 + l/betx0
    muIIf = muII0 + l/betay0
    return muIf , muIIf

## Objective function to be minimized
def objfunc(x):
    betax = x[5] # optic functions
    betay = x[5]
    alfax = 0.0
    alfay = 0.0
    gammax = (1+alfax**2)/betax
    gammay = (1+alfay**2)/betay
    betax0 = x[5]
    betay0 = x[5]
    alfax0 = 0.0
    alfay0 = 0.0
    gammax0 = (1+alfax0**2)/betax0,
    gammay0 = (1+alfay0**2)/betay0
    muxf = 0.0 # Phase advances
    muyf = 0.0
    muxfdesired = np.pi/2
    muyfdesired = np.pi/2
    dx = 0.0 # Dispersions
    dx0 = 0.0
    dy = 0.0
    dpy = 0.0
    dy0 = 0.0
    dpy0=0.0
    dpx = 0.0
    dpx0 = 0.0
    stepsize = 3 # stepsizes
    # Element definitions:
    # Lengths are divided by stepsize for integration purposes
    ## The use of this is to give an array of things to optimize, so array positions are focusing strenghts or drift lengths, it can also be beta functions
    D00 = Drift(l=x[4]/stepsize)
    Q00 = Quad(l=0.2/stepsize,k=x[0])
    Q01 = Quad(l=0.2/stepsize,k=x[1])
    Q02 = Quad(l=0.2/stepsize,k=x[2])
    Q03 = Quad(l=0.2/stepsize,k=x[3])
    MBI = Dipole(l=1.0/stepsize,rho=5.729577951308232) ## 10 degrees
    TEST = [D00,Q00,D00,Q01,D00,D00,Q02,D00,Q03,D00]
    extended = []
    ## This loop adds that many elements as there is stepsize!
    for element in TEST:
        for _ in range(stepsize):
            extended.append(element)
    ## This loop basically transforms the optic functions
    for element in extended:
        mat = element.t_matrix()
        lelement = element.lreturn()
        betax, alfax, gammax, betay, alfay, gammay = CStransfermat(betax,betay,alfax,alfay,gammax,gammay,mat)
        muxf,muyf = mutransport(muxf,muyf,betax,betay,lelement)
    ## This loop calculates the dispersion functions
    for element in extended:
        mat = element.t_matrix()
        if element==MBI:
            ## Here specify the bending radius of the dipole correctly and the length of the dipole as well!
            dx,dpx,dy,dpy = dispersiontransportdipole(dx,dpx,dy,dpy,mat,5.729577951308232,1.0/stepsize)
        else:
            dx,dpx,dy,dpy = dispersiontransport(dx,dpx,dy,dpy,mat)
    ## The following differenes are the things we wanna minimize
    diff1 = betax0 - betax
    diff2 = betay0 - betay
    diff3 = alfax0 - alfax
    diff4 = alfay0 - alfay
    diff5 = dx0 - dx
    diff6 = dpx0 - dpx
    diff7 = gammax0 - gammax
    diff8 = gammay0 - gammay
    diff9 =  muxfdesired - muxf
    diff10 =  muyfdesired - muyf
    objective = np.sum(diff1**2 + diff2**2 + diff3**2 + diff4**2 + diff7**2 + diff8**2) 
    return objective

## This function checks and later used to plot the optic functions
def testing(x):
    betax = x[5]
    betay = x[5]
    alfax = 0.0
    alfay = 0.0
    gammax = (1+alfax**2)/betax
    gammay = (1+ alfay**2)/betay
    dx = 0.0 
    dpx = 0.0
    dy = 0.0
    dpy = 0.0
    muxf = 0.0
    muyf = 0.0
    betaxarray = []
    betayarray = []
    dxarray = []
    dyarray = []
    larray = []
    stepsize = 3000
    D00 = Drift(l=x[4]/stepsize)
    Q00 = Quad(l=0.2/stepsize,k=x[0])
    Q01 = Quad(l=0.2/stepsize,k=x[1])
    Q02 = Quad(l=0.2/stepsize,k=x[2])
    Q03 = Quad(l=0.2/stepsize,k=x[3])
    MBI = Dipole(l=1.0/stepsize,rho=5.729577951308232) ## 10 degrees
    TEST = [D00,Q00,D00,Q01,D00,D00,Q02,D00,Q03,D00]
    Nturns = 1 ## Can specify number of repeated cells
    extended = []
    extendedturns = []
    betaxarray.append(betax)
    betayarray.append(betay)
    larray.append(0.0)
    dxarray.append(dx)
    dyarray.append(dy)
    for _ in range(Nturns):
        for element in TEST:
            extendedturns.append(element)
    for element in extendedturns:
        for _ in range(stepsize):
            extended.append(element)
    muxarray = [0]*(len(extended)+1)
    muyarray = [0]*(len(extended)+1)
    for element in extended:
        mat = element.t_matrix()
        lelement = element.lreturn()
        betax, alfax, gammax, betay, alfay, gammay = CStransfermat(betax,betay,alfax,alfay,gammax,gammay,mat)
        betaxarray.append(betax)
        betayarray.append(betay)
        larray.append(element.lreturn())
        muxf,muyf = mutransport(muxf,muyf,betax,betay,lelement)
    for element in extended:
        mat = element.t_matrix()
        if element==MBI:
            ## Here specify the bending radius of the dipole correctly and the length of the dipole as well!
            dx,dpx,dy,dpy = dispersiontransportdipole(dx,dpx,dy,dpy,mat,5.729577951308232,1.0/stepsize)
        else:
            dx,dpx,dy,dpy = dispersiontransport(dx,dpx,dy,dpy,mat)
        dxarray.append(dx)
        dyarray.append(dy)
    for i in range(len(betaxarray)):
        if i==0:
            muxarray[i] = 0.0
            muyarray[i] = 0.0
        else:
            betxinv = 1/betaxarray[i]
            betyinv = 1/betayarray[i]
            muxarray[i] = muxarray[i-1] + larray[i]*(betxinv)
            muyarray[i] = muyarray[i-1] + larray[i]*(betyinv)
    return larray, betaxarray, betayarray, dxarray, muxarray, muyarray,dyarray

## This function can be used to put a constraint on maximum and minimum boundaries through lattice:
def max_beta_constraint(x, max_beta=8.0, min_beta = 1.0):
    betax = 5.0
    betay = 5.0
    alfax = 0.0
    alfay = 0.0
    gammax = (1+alfax**2)/betax
    gammay = (1+ alfay**2)/betay
    stepsize = 1
    D00 = Drift(l=x[4]/stepsize)
    Q00 = Quad(l=0.2/stepsize,k=x[0])
    Q01 = Quad(l=0.2/stepsize,k=x[1])
    Q02 = Quad(l=0.2/stepsize,k=x[2])
    Q03 = Quad(l=0.2/stepsize,k=x[3])
    MBI = Dipole(l=1.0/stepsize,rho=5.729577951308232) ## 10 degrees
    TEST = [D00,Q00,D00,Q01,D00,D00,Q02,D00,Q03,D00]

    extended = []
    for element in TEST:
        for _ in range(stepsize):
            extended.append(element)
    for element in extended:
        mat = element.t_matrix()
        betax, alfax, gammax, betay, alfay, gammay = CStransfermat(betax,betay,alfax,alfay,gammax,gammay,mat)
        if betax > max_beta or betay > max_beta:
            return max_beta - max(betax, betay)  # return a negative value if constraint is violated
        if betax < min_beta or betay < min_beta:
            return min(betax,betay) - min_beta
    # for element in extended:
        # mat = element.t_matrix()
        # betax, alfax, gammax, betay, alfay, gammay = CStransfermat(betax, betay, alfax, alfay, gammax, gammay, mat)
        # if betax > max_beta or betay > max_beta:
            # return max_beta - max(betax, betay)  # return a negative value if constraint is violated

    return max_beta  # return positive value if constraint is satisfied


## One can find other optimizers on scipy.minimize function online!
## For constraint matching: 'SLSQP' 
## For normal matching: 'Nelder-Mead' || Simplex algo
## Matching function
def matchingfunc():
    numberofquads = 4
    numberofdrifts = 1
    lower_boundaryquad = -15.0
    upper_boundaryquad = 15.0
    constraints = [{'type': 'ineq', 'fun': max_beta_constraint, 'args': (1.8,1.3)}]
    x0 = [random.gauss(4.0,5.0) for _ in range(numberofquads)] + [0.2]*numberofdrifts + [5.0]
    result = minimize(objfunc, x0, tol=1e-10, method='Nelder-Mead',callback=callback, bounds=[(lower_boundaryquad, upper_boundaryquad)]*numberofquads + [(0.01,1.0)]*numberofdrifts + [(1.0,10.0)], options={'maxiter':50000}) 
    counter = 0
    while result.fun > 0.00000000001:
        x0 = [random.gauss(4.0,5.0) for _ in range(numberofquads)] + [0.2]*numberofdrifts + [5.0]
        result = minimize(objfunc, x0, tol=1e-10, method='Nelder-Mead',callback=callback, bounds=[(lower_boundaryquad, upper_boundaryquad)]*numberofquads + [(0.01,1.0)]*numberofdrifts + [(1.0,10.0)], options={'maxiter':50000}) 
        counter += 1
        if (counter % 100 ==0):
            print("===============================================")
            print("counter",counter)
            print("error func:",result.fun)
            print("===============================================")
            time.sleep(5)

    return result.fun, result.x

fun, x = matchingfunc()
print("=================== Results:")
print(fun)
print(x)
##======================================================
## Plotting:
ltot, betxarray, betyarray, dxarray, muxarr, muyarr, dyarray = testing(x)

ltotarray = [0]*len(ltot)
sumt = 0.0
for i in range(len(ltot)):
    sumt += ltot[i]
    ltotarray[i] = sumt


plt.plot(ltotarray,betxarray,color="red",label=r'$\beta_{x}$')
plt.plot(ltotarray,betyarray,color="blue",label=r'$\beta_{y}$')
plt.ylabel(r'$\beta_{x,y}$ m')
plt.xlabel("s[m]")
plt.legend()
plt.savefig("betafunctions.png")
plt.clf()

plt.plot(ltotarray,dxarray,color="red")
plt.plot(ltotarray,dyarray,color="green")
plt.ylabel(r'$D_{x,y}$ m')
plt.xlabel("s[m]")
plt.savefig("dispfunctions.png")
plt.clf()

plt.plot(ltotarray,muxarr,color="black",label=r'$\mu_{x}$')
plt.plot(ltotarray,muyarr,color="red",label=r'$\mu_{y}$')
plt.legend()
plt.ylabel(r'$\mu_{x,y}$ [rad]')
plt.xlabel("s [m]")
plt.savefig("phaseadvance.png")
plt.clf()


end_time = time.time()

print("Time took to finish in [secs]:", end_time - start_time)


