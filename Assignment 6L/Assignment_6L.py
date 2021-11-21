# EE2703
# Assignment 6L - Sarthak Vora (EE19B140)

# Run the code : python EE2703_ASSGN6L_EE19B140.py


import numpy as np
import scipy.signal as sp
import matplotlib.pyplot as plt





# Question 1
Num = np.poly1d([1, 0.5])
denom1 = np.poly1d([1, 0, 2.25])
denom2 = np.poly1d([1, 1, 2.5])
Den = np.polymul(denom1, denom2)
H = sp.lti(Num, Den)
t,x = sp.impulse(H,None, np.linspace(0,50,501))

# Plotting the signal
plt.plot(t,x,'r', label='decay=0.5')
plt.grid()
plt.legend(loc='upper right')
plt.xlabel(r"t$\rightarrow$")
plt.ylabel(r"x(t)$\rightarrow$")
plt.title('Solution for decay=0.5')
plt.show()






# Question 2
Num1 = np.poly1d([1, 0.05])
denom1 = np.poly1d([1, 0, 2.25])
denom2 = np.poly1d([1, 0.1, 2.2525])
Den1 = np.polymul(denom1, denom2)
H1 = sp.lti(Num1, Den1)
t1,x1 = sp.impulse(H1,None, np.linspace(0,50,501))

# Plotting the signal
plt.plot(t1,x1,'r', label='decay=0.05')
plt.grid()
plt.legend(loc='upper left')
plt.xlabel(r"t$\rightarrow$")
plt.ylabel(r"x(t)$\rightarrow$")
plt.title('Solution for decay=0.05')
plt.show()







# Question 3
for f in np.arange(1.4,1.6,0.05):
    Num2 = np.poly1d([1])
    Den2 = np.poly1d([1, 0, 2.25])
    Hf = sp.lti(Num2, Den2)
    t = np.linspace(0,50,501)
    u = np.ones_like(t)
    f_t = (np.cos(f*t))*(np.exp((-0.05)*t))*u
    tout, yout, xout = sp.lsim(Hf, f_t, t)
    
    # Plotting all the signals on the same plot
    plt.plot(tout, yout, label='w = '+str(f)+'rad/s')
    plt.xlabel(r"t$\rightarrow$")
    plt.ylabel(r"x(t)$\rightarrow$")
    plt.grid()
    plt.legend(loc='lower left')
    plt.title('Response of LTI system to various frequencies')
        
plt.show()







# Question 4
t = np.linspace(0,20,1000)
H_x = sp.lti(np.poly1d([1,0,2]), np.poly1d([1,0,3,0]))
t_out, x_t = sp.impulse(H_x,None, t)

# Plotting x(t)
plt.plot(t_out, x_t, 'b', label='x(t)')


H_y = sp.lti(np.poly1d([2]), np.poly1d([1,0,3,0]))
t_out, y_t = sp.impulse(H_y,None, t)

# Plotting y(t)
plt.plot(t_out, y_t, 'r', label='y(t)')
plt.xlabel(r"t$\rightarrow$")
plt.ylabel(r"Continous time response$\rightarrow$")
plt.title('Solution for coupled spring problem')
plt.legend()
plt.grid()
plt.show()









# Question 5
R = 100
L = 1e-6
C = 1e-6

H_num = np.poly1d([1])
H_den = np.poly1d([L*C, R*C, 1])
H_RLC = sp.lti(H_num, H_den)
w,S,phi = H_RLC.bode()

# Plotting the magnitude plot
plt.semilogx(w,S,'r')
plt.xlabel(r"$\omega$$\rightarrow$")
plt.ylabel(r"|H(jw)|$\rightarrow$")
plt.grid()
plt.title('Magnitude plot of the RLC network')
plt.show()

# Plotting the phase plot
plt.semilogx(w,phi,'r')
plt.xlabel(r"$\omega$$\rightarrow$")
plt.ylabel(r"$\angle$H(jw)(in $^{\circ}$)$\rightarrow$")
plt.grid()
plt.title('Phase plot of the RLC network')
plt.show()  








# Question 6
t_rlc_in = np.linspace(0, 30e-6, 301)
u_rlc = np.ones_like(t_rlc_in)
v_rlc = ((np.cos(1e3*t_rlc_in))*u_rlc) - ((np.cos(1e6*t_rlc_in))*u_rlc)
t_rlc_out, y_rlc, x_rlc = sp.lsim(H_RLC, v_rlc, t_rlc_in)

# Plotting vo(t)
plt.plot(t_rlc_out, y_rlc,'r', label=r'$v_{0}$(t)')
plt.xlabel(r"t$\rightarrow$")
plt.ylabel(r"$v_{0}(t)$$\rightarrow$")
plt.title(r"$v_{o}$(t) for t<30$\mu$s")
plt.legend()
plt.grid()
plt.show()

# Changing time range
t_rlc_in = np.linspace(0, 0.01, 100000)
u_rlc = np.ones_like(t_rlc_in)
v_rlc = ((np.cos(1e3*t_rlc_in))*u_rlc) - ((np.cos(1e6*t_rlc_in))*u_rlc)
t_rlc_out, y_rlc, x_rlc = sp.lsim(H_RLC, v_rlc, t_rlc_in)

# Plotting vo(t)
plt.plot(t_rlc_out, y_rlc,'r', label=r'$v_{0}$(t)')
plt.xlabel(r"t$\rightarrow$")
plt.ylabel(r"$v_{0}(t)$$\rightarrow$")
plt.title(r"$v_{o}$(t) for 0<t<10ms")
plt.legend()
plt.grid()
plt.show()

# Plotting the magnified plot for vo(t)
plt.plot(t_rlc_out, y_rlc,'r', label=r'$v_{0}$(t)')
plt.xlim([0.0062, 0.00668])
plt.ylim([0.960, 1.006])
plt.xlabel(r"t$\rightarrow$")
plt.ylabel(r"$v_{0}(t)$$\rightarrow$")
plt.title(r"Magnified plot for $v_{o}$(t): 0<t<10ms")
plt.legend()
plt.grid()
plt.show()



