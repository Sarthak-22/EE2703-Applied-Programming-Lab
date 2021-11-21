
# EE19B140 - Sarthak Vora
# Assignment 7
# Run the following on Anaconda prompt - python EE2703_ASSIGN7_EE19B140.py


import scipy.signal as sp
import numpy as np
import scipy 
import matplotlib.pyplot as plt
import sympy
sympy.init_session


#plotting helper function
def plotter(x,y,title,xlabel,ylabel):
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.plot(x,y,'r')
    plt.grid()
    plt.show()

#low pass filter definition function
def lowpass(R1,R2,C1,C2,G,Vi):
    s = sympy.symbols("s")
    A = sympy.Matrix([[0,0,1,-1/G],\
            [-1/(1+s*R2*C2),1,0,0],\
            [0,-G,G,1],\
            [-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
    b = sympy.Matrix([0,0,0,-Vi/R1])
    V = A.inv()*b
    return A,b,V


#high pass filter definition function
def highpass(R1,R3,C1,C2,G,Vi):
    s = sympy.symbols("s")
    A = sympy.Matrix([[0,-1,0,1/G],
        [s*C2*R3/(s*C2*R3+1),0,-1,0],
        [0,G,-G,1],
        [-s*C2-1/R1-s*C1,0,s*C2,1/R1]])
    b = sympy.Matrix([0,0,0,-Vi*s*C1])
    V = A.inv()*b
    return (A,b,V)


#this is a function that converts a sympy function to an expression that is understood by scipy.signals 
def symToTransferFn(Y):
    Y = sympy.expand(sympy.simplify(Y))
    n,d = sympy.fraction(Y)
    n,d = sympy.Poly(n,s), sympy.Poly(d,s)
    num,den = n.all_coeffs(), d.all_coeffs()
    num,den = [float(f) for f in num], [float(f) for f in den]
    return num,den


#plots step response
def stepresponse(Y,title):
    num,den = symToTransferFn(Y)
    den.append(0)
    H = sp.lti(num,den)
    t,y=sp.impulse(H,T = np.linspace(0,1e-3,10000)) 
    plotter(t,y,title,r'$t\rightarrow$',r'$V_{o}(t)\rightarrow$')
    return


#sum of sinusoids
def inputs(t):
    return (np.sin(2000*np.pi*t)+np.cos(2e6*np.pi*t))



#function that calculates response for arbitrary function
def inp_response(Y,title,inp=inputs,tlim=1e-3): 
    num,den = symToTransferFn(Y)
    H = sp.lti(num,den)
    t = np.linspace(0,tlim,100000)
    t,y,svec = sp.lsim(H,inp(t),t)
    plotter(t,y,title,r'$t\rightarrow$',r'$V_{o}(t)\rightarrow$')
    return



#High frequency damped sinusoids
def damped1(t,decay=3e3,freq=2e5*(np.pi)):
    return np.sin(freq*t)*np.exp(-decay*t) * (t>0)



#Low frequency damped sinusoids
def damped2(t,decay=1e1,freq=200*(np.pi)):
    return np.sin(freq*t)*np.exp(-decay*t) * (t>0)



#first we deal with lowpass filters
s =  sympy.symbols("s")


#defining the low pass transfer function
A,b,V=lowpass(10000,10000,1e-9,1e-9,1.586,1)
H=V[3]
print('The low-pass filter transfer function is -  ', end='')
print(H)
w = np.logspace(0,8,801)
ss = 1j*w
hf=sympy.lambdify(s,H,"numpy")
v=hf(ss)

#plotting low pass Magnitude response
plt.title("Bode plot of Low-pass filter transfer function")
plt.xlabel(r'$\omega \rightarrow$')
plt.ylabel(r'$|H(j\omega)|\rightarrow$')
plt.loglog(w,abs(v),'r',lw=2)
plt.grid(True)
plt.show()


#Problem 1, finding step response
stepresponse(H,r'Step response $V_{o}(t)$ for Low pass filter')



#Problem 2
t = np.linspace(0,1e-3,1000000)
plotter(t,inputs(t),r'$V_i(t)=(\sin(2x10^3\pi t)+\cos(2x10^6\pi t))u(t)$',r'$t\rightarrow$',r'$V_{i}(t)\rightarrow$')

inp_response(H,r'$V_o(t)$ for $V_i(t)=(\sin(2x10^3\pi t)+\cos(2x10^6\pi t))u(t)$ for Lowpass filter',inputs)




#High pass filters
A,b,V=highpass(10000,10000,1e-9,1e-9,1.586,1)
H=V[3]
print('The High-pass filter transfer function is - ',end='')
print(H)
w=np.logspace(0,8,801)
ss=1j*w
hf=sympy.lambdify(s,H,"numpy")
v=hf(ss)


plt.title("Bode plot of High-pass filter transfer function")
plt.xlabel(r'$\omega \rightarrow$')
plt.ylabel(r'$|H(j\omega)|\rightarrow$')
plt.loglog(w,abs(v),'r',lw=2)
plt.grid(True)
plt.show()

#plotting step response
stepresponse(H,r'Step response $V_{o}(t)$ for High-pass filter')


#plotting sum of sinusoids 
inp_response(H,r'$V_{o}(t)$ for $V_{i}(t)=(sin(2x10^3\pi t)+cos(2x10^6\pi t))u(t)$ for Highpass filter',inputs,tlim= 1e-5)


#plotting response to damped sinusoids
inp_response(H,r'$V_{o}(t)$ for $V_{i}(t)=sin(2x10^5\pi t)e^{-3000t}u(t)$ for Highpass filter',damped1)


inp_response(H,r'$V_{o}(t)$ for $V_{i}(t)=sin(200\pi t)e^{-1000t}u(t)$ for Highpass filter',damped2,tlim = .5)


#Plotting high damped sinusoid
t = np.linspace(0,1e-3,100000)
plt.title(r'High frequency damped sinusoid: $V_{i}(t)=sin(2x10^5\pi t)e^{-3000t}u(t)$')
plt.xlabel(r'$t\rightarrow$')
plt.ylabel(r'$v_{i}(t)\rightarrow$')
plt.plot(t,damped1(t,decay=3e3,freq=2e5*(np.pi)),'r')
plt.grid()
plt.show()

#Plotting low damped sinusoid
t = np.linspace(0,0.5,100000)
plt.title(r'Low frequency damped sinusoid: $V_{i}(t)=sin(200\pi t)e^{-1000t}u(t)$')
plt.xlabel(r'$t\rightarrow$')
plt.ylabel(r'$v_{i}(t)\rightarrow$')
plt.plot(t,damped2(t),'r')
plt.grid()
plt.show()