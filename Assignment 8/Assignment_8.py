# EE2703 : Assignment 8

# Sarthak Vora - EE19B140

# Run the code: python EE2703_ASSIGN8_EE19B140.py



from pylab import *


# Corrected spectrum of sin(5t)
t=linspace(0,2*pi,129);t=t[:-1]
y=sin(5*t)
Y=fftshift(fft(y))/128.0
w=linspace(-64,63,128)

figure()
subplot(2,1,1)
plot(w,abs(Y),'b',lw=2)
xlim([-10,10])
ylabel(r'|Y|$\rightarrow$')
title(r'Spectrum of $\sin(5t)$')
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-10,10])
ylabel(r'$\angle Y\rightarrow$')
xlabel(r'$\omega\rightarrow$')
grid(True)
show()






# Corrected spectrum of (1 + 0.1cos(t))cos(10t)
t=linspace(-4*pi,4*pi,513);t=t[:-1]
y=(1+0.1*cos(t))*cos(10*t)
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513); w = w[:-1]

figure()
subplot(2,1,1)
plot(w,abs(Y),'b',lw=2)
xlim([-15,15])
ylabel(r'|Y|$\rightarrow$')
title(r'Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$')
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-15,15])
ylabel(r'$\angle Y\rightarrow$')
xlabel(r'$\omega\rightarrow$')
grid(True)
show()






# Spectrum of sin^3(t)
t=linspace(-4*pi,4*pi,513);t=t[:-1]
y=pow(sin(t),3)
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513); w = w[:-1]

figure()
subplot(2,1,1)
plot(w,abs(Y),'b',lw=2)
xlim([-15,15])
ylabel(r'|Y|$\rightarrow$')
title(r'Spectrum of $\sin^{3}(t)$')
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-15,15])
ylabel(r'$\angle Y\rightarrow$')
xlabel(r'$\omega\rightarrow$')
grid(True)
show()






# Spectrum of cos^3(t)
t=linspace(-4*pi,4*pi,513);t=t[:-1]
y=pow(cos(t),3)
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513); w = w[:-1]

figure()
subplot(2,1,1)
plot(w,abs(Y),'b',lw=2)
xlim([-15,15])
ylabel(r'|Y|$\rightarrow$')
title(r'Spectrum of $\cos^{3}(t)$')
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-15,15])
ylabel(r'$\angle Y\rightarrow$')
xlabel(r'$\omega\rightarrow$')
grid(True)
show()





# Spectrum of cos(20t + 5cos(t))
t=linspace(-4*pi,4*pi,513);t=t[:-1]
y=cos(20*t + 5*cos(t))
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513); w = w[:-1]

figure()
subplot(2,1,1)
plot(w,abs(Y),'b',lw=2)
xlim([-40,40])
ylabel(r'|Y|$\rightarrow$')
title(r'Spectrum of $\cos(20t + 5\cos(t))$')
grid(True)
subplot(2,1,2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-40,40])
ylabel(r'$\angle Y\rightarrow$')
xlabel(r'$\omega\rightarrow$')
grid(True)
show()





# Estimated Spectrum of Gaussian
t=linspace(-4*pi,4*pi,513);t=t[:-1]
y=exp(-(t**2)/2)
Y=fftshift(abs(fft(ifftshift(y))))/512.0
Y=Y*sqrt(2*pi)/max(Y)
Y_ = Y.copy()
w=linspace(-64,64,513); w = w[:-1]

figure()
subplot(2,1,1)
plot(w,abs(Y),'b',lw=2)
xlim([-10,10])
yticks([0,0.5,1,1.5,2,2.5])
ylabel(r'|Y|$\rightarrow$')
title(r'Estimated Spectrum of $\exp(-t^{2}/2)$')
grid(True)
subplot(2,1,2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-10,10])
ylabel(r'$\angle Y\rightarrow$')
xlabel(r'$\omega\rightarrow$')
grid(True)
show()




# True spectrum of Gaussian
t=linspace(-4*pi,4*pi,513);t=t[:-1]
y=exp(-(t**2)/2)
Y=fftshift(abs(fft(y)))/512.0
Y=Y*sqrt(2*pi)/max(Y)
w=linspace(-64,64,513); w = w[:-1]

figure()
subplot(2,1,1)
plot(w,abs(Y),'b',lw=2)
xlim([-10,10])
yticks([0,0.5,1,1.5,2,2.5])
ylabel(r'|Y|$\rightarrow$')
title(r'True Spectrum of $\exp(-t^{2}/2)$')
grid(True)
subplot(2,1,2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-10,10])
ylabel(r'$\angle Y\rightarrow$')
xlabel(r'$\omega\rightarrow$')
grid(True)
show()

error = max(abs(Y_ - Y))




t=linspace(-8*pi,8*pi,513);t=t[:-1]
y=exp(-(t**2)/2)
Y=fftshift(abs(fft(y)))/512.0
Y=Y*sqrt(2*pi)/max(Y)
w=linspace(-64,64,513); w = w[:-1]

figure()
subplot(2,1,1)
plot(w,abs(Y),'b',lw=2)
xlim([-10,10])
yticks([0,0.5,1,1.5,2,2.5])
ylabel(r'|Y|$\rightarrow$')
title(r'Spectrum of $\exp(-t^{2}/2)$ in the interval ($-8\pi$,$8\pi$)')
grid(True)
subplot(2,1,2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-10,10])
ylabel(r'$\angle Y\rightarrow$')
xlabel(r'$\omega\rightarrow$')
grid(True)
show()
