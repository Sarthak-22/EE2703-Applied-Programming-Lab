



from pylab import *
from numpy import *
import scipy.special as sp
from scipy.linalg import *
import sys


#Run the following command in terminal : 'python generate_data.py'


#Loading the dataframe
try:
    df = loadtxt(r"fitting.dat", delimiter=' ')
    df = transpose(df) # 1st row is the time row i.e df[0]
    sigma=logspace(-1,-3,9) # Define the standard deviation.
except:
    sys.exit('fitting.dat not found. Run \'python generate_data.py\' and try again')
    


#Define the 2nd order Bessel Function
def g(t,A,B): 
    return (A*sp.jv(2,t)) + B*t



#Question 3: Code for plotting Figure 0
for i in range(1,10):
    plot(df[0],df[i], label="$\sigma$"+str(i)+"="+str(round(sigma[i-1],3)))
plot(df[0], g(df[0],1.05,-0.105), label='True value', color='black') # Plot for Question 4: True value function.
legend(loc='upper right')
title('Noisy plots v/s True plots')
xlabel('t')
ylabel('f(t) + noise')
grid(b=True, which='major', color='#666666', linestyle='-')
show()



#Question 5: Error bar plot
errorbar(df[0][::5],df[1][::5],sigma[0],fmt='ro', label='errorbar') # Plotting at an interval of 5 values
plot(df[0], g(df[0],1.05,-0.105), label='True value', color='black')
legend()
xlabel(r't $\rightarrow$')
ylabel('f(t)')
title('Errorbar Plot')
grid(b=True, which='major', color='#666666', linestyle='-')
show()




#Question 6: Proving that two vectors are equal
a, b = 1.05, -0.105
M = c_[sp.jv(2,df[0]), df[0]] # concatenating(stacking) the two vectors vertically.
p = hstack((a,b))
x = dot(M,p) # Multiplying 2 matrices
y =  where(x!=g(df[0],a,b)) # Find where are the 2 values not equal
print("The number of positions at which the 2 vectors are not equal is " + str(len(y)-1))  # empty numpy array has one element




#Question 7: Finding the mean squared error
A = linspace(0,2,21) # Defining the coefficients
B = linspace(-0.2,0,21)
epsilon = zeros((len(A),len(B))) # Initializing the error
for i in range(0,len(A)):
    for j in range(0,len(B)):
        epsilon[i,j] = sum(square(df[1] - g(df[0],A[i],B[j])))
        epsilon[i,j] = epsilon[i,j]/101




#Question 8: Contour plot of MSerror
C = contour(A,B,epsilon, levels=20)
clabel(C, [0.   , 0.025, 0.05, 0.1])
scatter(1.05, -0.105, color='red', label='Exact location')
annotate("Exact location", (1.05,-0.105))
xlabel(r'A $\rightarrow$')
ylabel(r'B $\rightarrow$')
title(r'Contours of $\epsilon_{ij}$')
show()




#Question 9: Least squares program 
p,*rest = lstsq(M, x) # p contains the fitting(optimum) value of A and B inside a list




#Question 10: Plotting the estimation error
Aerr = zeros((9,1)) # Initialize the error values
Berr = zeros((9,1))
for i in range(1, 10):
    predoutput,*rest = lstsq(M, df[i])
    Aerr[i-1] = abs(predoutput[0] - p[0]) # # Calculating Aerr and Berr for each lstsq estimation
    Berr[i-1] = abs(predoutput[1] - p[1]) 

plot(sigma, Aerr, 'o--', label='Aerr')
plot(sigma, Berr, 'o--', label='Berr')
xlabel(r'Noise standard deviation $\rightarrow$')
ylabel(r'MS error $\rightarrow$')
legend()
grid()
title('Variation of error with noise')
show()




#Question 11: loglog plot for error estimation
loglog(sigma, Aerr,'ro', label='Aerr')
loglog(sigma, Berr, 'bo', label='Berr')
legend(loc='upper right')
grid(b=True)
errorbar(sigma, Aerr, std(Aerr), fmt='ro')
errorbar(sigma, Berr, std(Berr), fmt='bo')
xlabel(r"$\sigma_n$")
ylabel(r"MSerror $\rightarrow$")
title('Variation of error with noise')
show()




