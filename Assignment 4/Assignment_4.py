
from pylab import *
from scipy.integrate import quad
from scipy.linalg import lstsq


# Question 1: Define the functions and plotting them

def exponential(x): # Defining exponential function
    out = exp(x) 
    return out

# Defining cos(cos(x)) function
def coscos(x):
    out_i = cos(x)
    out_f = cos(out_i)
    return out_f



x_range = linspace(-2*pi, 4*pi, 1201)

figure('Figure 1')
# Actual function: exp(x)
semilogy(x_range, exponential(x_range), label='Actual plot: exp(x)')
# Periodic extension
semilogy(x_range, exponential(x_range%(2*pi)), 'r--', label='periodic extension')
title('Figure 1: exp(x)')
xlabel(r'x$\rightarrow$')
ylabel(r'y$\rightarrow$')
legend(loc='best')
grid()
show()




figure('Figure2')
# Actual function: coscos(x)
plot(x_range, coscos(x_range), label='Actualplot: cos(cos(x))')
# Periodic extension
plot(x_range, (coscos(x_range)%(2*pi)),'r--', label='periodic extension')
ylim([0.5, 1.1])
title('Figure 2: cos(cos(x))')
xlabel(r'x$\rightarrow$')
ylabel(r'y$\rightarrow$')
legend(loc='upper right')
grid()
show()




# Question2: Obtaining 51 Fourier series coefficients for both the functions

def u(x,k,f): # defining functions for integration
    out_1 = f(x)*cos(k*x)
    return out_1

def v(x,k,f):
    out_2 = f(x)*sin(k*x)
    return out_2




# Defining a function that returns 51 coeffs
def get51coeffs(f):
    acoeffs = list()
    bcoeffs = list()
    acoeffs.append(quad(u,0,2*pi,args=(0,f))[0]/(2*pi))
    bcoeffs.append(0)
    for i in linspace(1,25,25):
        acoeffs.append(quad(u,0,2*pi,args=(i, f))[0]/pi)
        bcoeffs.append(quad(v,0,2*pi,args=(i, f))[0]/pi)
        
    coeffs = [0 for i in linspace(0,50,51)]
    coeffs[0] = acoeffs[0]
    coeffs[1::2] = acoeffs[1:]
    coeffs[2::2] = bcoeffs[1:]
    return coeffs

exp_coeff = get51coeffs(exponential)
coscos_coeff = get51coeffs(coscos)
exp_coeff = [abs(a) for a in exp_coeff] # Taking modulus of all the elements





# Semilogy for Fourier series coeffs of exp(x)
figure('Figure3: semilog plot: $e^x$')
semilogy(0,exp_coeff[0],'ko',label='$a_{0}$')
semilogy(linspace(1,25,25),exp_coeff[1::2],'ro',label='$a_{n}$')
semilogy(linspace(1,25,25),exp_coeff[2::2],'bo',label='$b_{n}$')
#annotate('Exact location',(0,aexp[0]))
title('Figure 3: Semilogy plot for $e^x$')
xlabel(r'n$\rightarrow$')
ylabel(r'semilog(y)$\rightarrow$')
grid()
legend()
show()



# loglog plot for Fourier series coeffs of exp(x)
figure('Figure4: loglog plot: $e^x$')
loglog(0,exp_coeff[0],'ko',label='$a_{0}$')
loglog(linspace(1,25,25),exp_coeff[1::2],'ro',label='$a_{n}$')
loglog(linspace(1,25,25),exp_coeff[2::2],'bo',label='$b_{n}$')
title('Figure 4: loglog plot for $e^x$')
xlabel(r'log(n)$\rightarrow$')
ylabel(r'log(y)$\rightarrow$')
grid()
legend()
show()



# Semilogy for Fourier series coeffs of coscos(x)
figure('Figure5: semilog plot: $cos(cos(x))$')
semilogy(0,coscos_coeff[0],'ko',label='$a_{0}$')
semilogy(linspace(1,25,25),coscos_coeff[1::2],'ro',label='$a_{n}$')
semilogy(linspace(1,25,25),coscos_coeff[2::2],'bo',label='$b_{n}$')
title('Figure 5: Semilogy plot for $cos(cos(x))$')
xlabel(r'n$\rightarrow$')
ylabel(r'semilog(y)$\rightarrow$')
grid()
legend()
show()



# loglog plot for Fourier series coeffs of coscos(x)
figure('Figure6: loglog plot: $cos(cos(x))$')
loglog(0,coscos_coeff[0],'ko',label='$a_{0}$')
loglog(linspace(1,25,25),coscos_coeff[1::2],'ro',label='$a_{n}$')
loglog(linspace(1,25,25),coscos_coeff[2::2],'bo',label='$b_{n}$')
title('Figure 6: loglog plot for $cos(cos(x))$')
xlabel(r'log(n)$\rightarrow$')
ylabel(r'log(y)$\rightarrow$')
grid()
legend()
show()



# Question4: Least Squares Approach of the given functions

x = linspace(0,2*pi,401) # Vector with 401 steps
x = x[:-1] # dropping the last term for periodic integral

b1 = exponential(x) # vector entries for exponential(x)
b2 = coscos(x) # vector entries for cos(cos(x))

A = zeros((400,51)) 
A[:,0] = 1 # 1st col is all 1

range_col = [int(el) for el in linspace(1,25,25)]

for k in range_col :
    A[:,2*k-1] = cos(k*x) #cos(kx) column
    A[:,2*k] = sin(k*x) #sin(kx) column

c1 = lstsq(A,b1)[0] # best fit coeffeicients for exp(x)
c2 = lstsq(A,b2)[0] # best fit coefficients for cos(cos(x))



# Comparing the semilog plot of exp(x) with best fit and integration
figure('Figure3.1: Comparing semilog plot for FS Coeff: $e^x$')
semilogy(linspace(1,25,25),exp_coeff[1::2],'ro',label='$a_{n}$ by integration')
semilogy(linspace(1,25,25),exp_coeff[2::2],'bo',label='$b_{n}$ by integration')
semilogy(linspace(1,25,25),abs(c1[1::2]),'go',label='$a_{n}$ by best fit')
semilogy(linspace(1,25,25),abs(c1[2::2]),'yo',label='$b_{n}$ by best fit')
title('Figure3.1: Comparing Semilogy plots: $e^x$ FS coefficients')
xlabel(r'n$\rightarrow$')
ylabel(r'semilog(y)$\rightarrow$')
grid(True)
legend()
show()



# Comparing the loglog plot of exp(x) with best fit and integration
figure('Figure4.1: Comparing loglog plot for FS Coeff: $e^x$')
loglog(linspace(1,25,25),exp_coeff[1::2],'ro',label='$a_{n}$ by integration')
loglog(linspace(1,25,25),exp_coeff[2::2],'bo',label='$b_{n}$ by integration')
loglog(linspace(1,25,25),abs(c1[1::2]),'go',label='$a_{n}$ by best fit')
loglog(linspace(1,25,25),abs(c1[2::2]),'yo',label='$b_{n}$ by best fit')
title('Figure4.1: Comparing Loglog plots: $e^x$ FS coefficients')
xlabel(r'log(n)$\rightarrow$')
ylabel(r'log(y)$\rightarrow$')
grid(True)
legend()
show()



# Comparing the Semilogy plots of cos(cos(x)) with best fit and integration
figure('Figure5.1: Comparing semilog plots for FS coeff: $cos(cos(x))$')
semilogy(0,coscos_coeff[0],'ko')
semilogy(linspace(1,25,25),coscos_coeff[1::2],'ro',label='$a_{n}$ by integration')
semilogy(linspace(1,25,25),coscos_coeff[2::2],'bo',label='$b_{n}$ by integration')
semilogy(linspace(1,25,25),abs(c2[1::2]),'go',label='$a_{n}$ by best fit')
semilogy(linspace(1,25,25),abs(c2[2::2]),'yo',label='$b_{n}$ by best fit')
title('Figure5.1: Comparing Semilogy plots: $cos(cos(x))$ FS coefficients')
xlabel(r'n$\rightarrow$')
ylabel(r'semilog(y)$\rightarrow$')
grid()
legend()
show()



# Comparing the loglog plots of cos(cos(x)) with best fit and integration
figure('Figure6.1: Comparing loglog plots for FS coeff: $cos(cos(x))$')
loglog(0,coscos_coeff[0],'ko')
loglog(linspace(1,25,25),coscos_coeff[1::2],'ro',label='$a_{n}$ by integration')
loglog(linspace(1,25,25),coscos_coeff[2::2],'bo',label='$b_{n}$ by integration')
loglog(linspace(1,25,25),abs(c2[1::2]),'go',label='$a_{n}$ by best fit')
loglog(linspace(1,25,25),abs(c2[2::2]),'yo',label='$b_{n}$ by best fit')
title('Figure 6.1: Comparing loglog plots: $cos(cos(x))$ FS coefficients')
xlabel(r'log(n)$\rightarrow$')
ylabel(r'log(y)$\rightarrow$')
grid()
legend()
show()



# Question 6: Finding error(deviation) in FS Coeffs by integration and best fit
abserr_coscos = max([abs(coscos_coeff[int(i)] - c2[int(i)]) for i in linspace(0,50,51)])
abserr_exp = max([abs(exp_coeff[int(i)] - c1[int(i)]) for i in linspace(0,50,51)])
print("The max deviation for cos(cos(x)) FS coefficients is: " + str(abserr_coscos))
print("The max deviation for exp(x) FS coefficients is: " + str(abserr_exp))



# Question 7: Plotting the function using FS coefficients by integration and best fit and comparing them
new_x = linspace(-2*pi,4*pi,1201)
new_x = new_x[:-1] # Removing the last elements to include periodic elements only

A_new = zeros((1200,51)) 
A_new[:,0] = 1 # 1st col is all 1

range_col = [int(el) for el in linspace(1,25,25)]

for k in range_col :
    A_new[:,2*k-1] = cos(k*new_x) #cos(kx) column
    A_new[:,2*k] = sin(k*new_x) #sin(kx) column

# Plotting exp(x) from its FS coefficients by integration and best fit
figure('Figure 7: Plotting exp(x) from FS coefficients')
plot(new_x, exponential(new_x%(2*pi)), 'go', label='Integration coeffs')
plot(new_x, A_new@c1, 'r', label='Best fit coeffs')
legend(loc='best')
title('Figure 7: Plotting $e^x$ from FS coefficients')
xlabel(r'x$\rightarrow$')
ylabel(r'y$\rightarrow$')
grid(True)
show()




# Plotting cos(cos(x)) from its FS coefficients by integration and best fit
figure('Figure 8: Plotting cos(cos(x)) from FS coefficients')
plot(new_x, dot(A_new,coscos_coeff), 'go', label='Integration coeffs')
plot(new_x, dot(A_new,c2), 'r', label='Best fit coeffs')
legend(loc='upper right')
title('Figure 8: Plotting $cos(cos(x))$ from FS coefficients')
xlabel(r'x$\rightarrow$')
ylabel(r'y$\rightarrow$')
grid(True)
show()

