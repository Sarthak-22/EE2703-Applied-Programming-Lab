from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import scipy.linalg
import sys


# Enter the command as : `python EE19B140_assgn5.py Nx Ny radius Niter`
# Only if all 4 parameters are mentioned in the command, custom values would be considered. For user input params less than 4, defualt values would be considered

# Note: Ensure giving proportional values of Nx,Ny and radius to ensure a constant radius of 0.35 in all cases


# get user inputs
if(len(sys.argv)==5):
    Nx = int(sys.argv[1])
    Ny = int(sys.argv[2])
    radius = int(sys.argv[3])
    Niter = int(sys.argv[4])
    print("Using user provided params")
    
else:
    # Initializing the default parameters
    Nx=25
    Ny=25
    radius = 8
    Niter = 1500 
    print("Using default values. Input all 4 params to use custom values")
    





# Initializing the Electric Potential
phi = zeros((Ny,Nx))



# Defining the coordinates and scale
x = linspace(-0.5,0.5,Nx)
y = linspace(0.5,-0.5,Ny)
Y,X = meshgrid(x,y)



# Finding the desired coordinates
ii = where(square(X)+square(Y)<=(0.35)**2) #round((radius/(Nx-1))**2,2))
phi[ii]=1.0



# Plotting Contour Plot
figure('Figure 1: Contour Plot of potential')
contourf(phi,cmap=cm.jet)
title('Figure 1: Contour plot of potential $\phi$')
xlabel(r'x $\rightarrow$')
ylabel(r'y $\rightarrow$')
show()




# Initialing the error and finding the error
errors = zeros((Niter))
for i in range(Niter):
    oldphi = phi.copy()
    
    #Updating potential
    phi[1:-1,1:-1] = 0.25*(phi[1:-1,0:-2] + phi[1:-1,2:] + phi[0:-2,1:-1] + phi[2:,1:-1]) 
    
    # Boundary conditions for potential
    phi[1:-1,0] = phi[1:-1,1] # Left Edge 
    phi[1:-1,Nx-1] = phi[1:-1,Nx-2] # Right edge
    phi[0,:] = phi[1,:] # Top Edge
    # Bottom Edge is already grounded
    
    phi[ii] = 1.0 # Assigning 1V to electrode region
    
    # Calculating error
    errors[i] = (abs(phi-oldphi)).max()

iterations = arange(0,Niter,1)




# Plotting loglog plot
figure('Figure 2: loglog plot of Error')
loglog(iterations,errors,'r')
ylabel(r'log(Error) $\rightarrow$')
xlabel(r'log(iterations) $\rightarrow$')
title('Figure 2: Loglog plot of error')
grid()
show()



# Plotting semilog plot
figure('Figure 3: Semilog plot of error')
semilogy(iterations,errors,'r')
ylabel(r'log(Error) $\rightarrow$')
xlabel(r'No. of iterations $\rightarrow$')
title('Figure 3: Semilog plot of error')
grid()
show()




# Defining the best fit function
def fit(iterations,errors):
    coeff = zeros((len(errors),2))
    const = zeros((len(errors),1))
    coeff[:,0] = 1
    coeff[:,1] = iterations
    const = log(errors)
    
    fit = scipy.linalg.lstsq(coeff, const)[0]
    estimate = dot(coeff,fit)
    
    return fit, estimate




# Finding the fit and estimate of the fit for all values and the last 500 values
fitAll, estimate1 = fit(iterations, errors)
fitAfter500, estimate2 = fit(iterations[501:], errors[501:])
print('A = '+ str(exp(fitAll[0]))+' B = '+str(fitAll[1]) + ' from Fit1')
print('A = '+ str(exp(fitAfter500[0]))+' B = '+str(fitAfter500[1]) + ' from Fit2')




# Plotting the comaprision between actual error and fits in semilog 
figure('Figure 4: Comparision of actual error and fits (Semilog)')
semilogy(iterations,errors,'g',markersize=6,label='Actual Error')
semilogy(iterations[::50],exp(estimate1[::50]),'ro',markersize=8,mfc=None,label='Total fit')
semilogy(iterations[501::50], exp(estimate2[::50]),'bo',markersize=5.5,mfc=None,label='Fit after 500')
ylabel(r'log(Error) $\rightarrow$')
xlabel(r'Iterations $\rightarrow$')
title('Figure 4: Comparision of actual Error v/s Fits (Semilog)')
legend()
grid()
show()



# Defining stopping condition error function
def cumerror(error,N,A,B):
    return -(A/B)*exp(B*(N+0.5))



# Defining a function to find stopping condition
def findStopCondn(errors, Niter, error_tol):
    cum_Err=[]
    for n in range(Niter):
        cum_Err.append(cumerror(errors[n],n, exp(fitAll[0]), fitAll[1]))
        if(cum_Err[n-1] <= error_tol):
            print("last per-iteration change in the error is ", (np.abs(cum_Err[-1]-cum_Err[-2])))
            return cum_Err[n-1], n
    print("last per-iteration change in the error is ", (np.abs(cum_Err[-1]-cum_Err[-2])))
    return cum_Err[-1], Niter




# Defining error tolerance and finding the stopping condition and error for N=1500(default)
errorTol = 10e-8
cum_Err, nStop = findStopCondn(errors, Niter, errorTol)
print("Stopping Condition ----> N: %g and Error: %g" % (nStop, cum_Err))




# Plotting the 3-D Surface Potential
fig5=figure('Figure 5: 3-D plot of the potential') # open a new figure
ax=p3.Axes3D(fig5) # Axes3D is the means to do a surface plot
title('Figure 5: 3-D surface plot of Potential $\phi$')
surf = ax.plot_surface(Y, X, phi, rstride=1,cstride=1, cmap=cm.jet)
xlabel(r'$\leftarrow $X $\rightarrow$')
ylabel(r'$\leftarrow $Y $\rightarrow$')
ax.set_zlabel('$Z$')
show()




# Plotting the contour plot of potential phi
figure('Figure 6: Contour Plot of Potential $\phi$')
plot(x[ii[0]],y[ii[1]],'ro', label='Central lead: 1V region')
contour(Y,X,phi)
xlabel(r'X $\rightarrow$')
ylabel(r'Y $\rightarrow$')
title('Figure 6: Contour plot of potential $\phi$')
legend()
show()




# Initializing and defining the current densities
Jx = zeros((Ny,Nx))
Jy = zeros((Ny,Nx))
Jx[1:-1,1:-1] = 0.5*(phi[1:-1,0:-2]-phi[1:-1,2:])
Jy[1:-1,1:-1] = 0.5*(phi[2:,1:-1]-phi[0:-2,1:-1])



# Plotting the vector plot of current flow
figure('Figure 7: Vector Plot of the Current Flow')
plot(x[ii[0]],y[ii[1]],'ro', label='Central lead: 1V region')
quiver(x,y,Jx,Jy, label='Current density')
xlabel(r'X $\rightarrow$')
ylabel(r'Y $\rightarrow$')
title('Figure 7: Vector Plot Of Current flow')
legend(loc='upper right')
show()




