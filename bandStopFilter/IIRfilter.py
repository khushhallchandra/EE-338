#' % Bandstop IIR filter design
#' % Khushhall Chandra Mahajan
#' % 7th April 2016
#' M = 141

import numpy as np
import pylab
import scipy.signal as sg
import matplotlib.pyplot as plt
import pylab as pl


#' #Given specifications
M = 141
delta_p = 0.15
delta_s = 0.15
f_sample = 100000

omega_p1 = 9.4*1000
omega_s1 = 11.4*1000
omega_s2 = 21.4*1000
omega_p2 = 23.4*1000

#' # Conversion of the given frequency
#' 1. Normalize the given frequency
#' 2. Convert the frequency from digital domain to analog domain
#' 3. Bandpass to Lowpass conversion
omega = np.array([omega_p1,omega_s1,omega_s2,omega_p2],dtype='f')

# Normalized digital frequencies
omega_normalized=(omega/f_sample)*2*np.pi
omega_analog = np.tan(omega_normalized/2)

omega_o = omega_analog[0] * omega_analog[3]
B = omega_analog[3] - omega_analog[0]

#Frequency Transformation
omega_l = ( B*omega_analog)/(( omega_analog**2)-omega_o)
omega_lpf = min(-1*omega_l[1],omega_l[2])
W_s = omega_lpf
W_p = omega_l[3]
#print omega_l
#' # Designing the low pass filter
D1 = (1/(1-delta_p)**2)-1
D2=(1/delta_s**2)-1

N = np.log10(D2/D1)/np.log10(W_s/W_p)
N = np.ceil(N/2)
W_c = W_p / np.power(D1,(0.5/N))

#Finding the filter poles
# Total 2N poles, out of which take N 
left_poles = np.zeros([N],dtype='complex64')
for k in range(int(N)):
	left_poles[k] = W_c*np.exp(1.j*np.pi*(2*k + 1 + N)/(2*N))

# Plotting poles of analog low pass filter

plt.figure(1)
plt.grid(True)	
plt.scatter(left_poles.real,left_poles.imag,s=50,marker='x')
plt.title('Pole zero plot of Low Pass Filter')
plt.ylabel('Imaginary')
plt.xlabel('Real')
#gain calculation
gain = np.power(W_c,N)

#' # Convert back to analog bandpass domain
#-------------------------------------------------------------
numerator = gain
denominator = 1+0.j
for p in left_poles:
    numerator = np.poly1d([1,0,omega_o ],r=0)*numerator 
    denominator = np.poly1d([-p,B, -p*omega_o],r=0)*denominator
z,p,k = sg.tf2zpk(numerator,denominator)


plt.figure(2)
plt.grid(True)
plt.scatter(p.real,p.imag,s=50,c='b',marker='x')
plt.scatter(z.real,z.imag,s=50,c='b',marker='o')
plt.title('Pole Zero plot of Analog Bandstop filter')
plt.ylabel('Imaginary')
plt.xlabel('Real')

#' # Convert back to digital domain
#-------------------------------------------------------------
numerator = gain
denominator = 1+0.j
for p in left_poles:
    numerator = np.poly1d([1+omega_o, -2+(2*omega_o), 1+omega_o],r=0)*numerator 
    denominator = np.poly1d([-p+B-(p*omega_o), (2*p)-(2*p*omega_o) , -p-B-(omega_o*p)],r=0)*denominator
z,p,k = sg.tf2zpk(numerator,denominator)
#print "Check",z,p,k
plt.figure(3)
plt.grid(True)
plt.scatter(p.real,p.imag,s=50,marker='x')
plt.scatter(z.real,z.imag,s=50,marker='o')
plt.title('Pole Zero plot of Digital Bandstop filter')
plt.ylabel('Imaginary')
plt.xlabel('Real')

#Plot Frequency response
plt.figure(4)
plt.clf()
plt.grid(True)
w,h= sg.freqz(numerator,denominator)
plt.plot((w/np.pi)*f_sample/2, np.absolute(h), linewidth=2)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain')
plt.title('Frequency Response')

# plt.show()

#' # The resulting transfer function:
#' \begin{equation}
#' H(z) = \frac{ 14.37z^{18}  - 147.1z^{17}  + 798.4z^{16}  - 2952z^{15}  + 8229z^{14}  - (1.821e+04)z^{13} + (3.301e+04)z^{12}  - (4.991e+04)z^{11}  + (6.37e+04)z^{10}  - (6.905e+04)z^9 + (6.37e+04)z^8 - (4.991e+04)z^7 + (3.301e+04)z^6 - (1.821e+04)z^5 + 8229z^4 - 2952z^3 + 798.4z^2 - 147.1z + 14.37}{(177.7 + -2.842e-14j)z^{18} + (-1331 + 2.274e-13j)z^{17}  + 5215z^{16} - (1.394e+04)z^{15}  + (2.822e+04 + 1.819e-12j)z^{14} + (-4.571e+04 + 1.819e-12j)z^{13}  + (6.11e+04 + -3.638e-12j)z^{12} + (-6.869e+04 + 1.091e-11j)z^{11}  + (6.569e+04 + -1.273e-11j)z^{10} + (-5.376e+04 + 1.364e-11j)z^9 + (3.772e+04 + -1.296e-11j)z^8 + (-2.263e+04 + 7.901e-12j)z^7 + (1.153e+04 + -3.865e-12j)z^6 + (-4932 + 1.478e-12j)z^5 + (1737 + -5.116e-13j)z^4 + (-488 + 1.279e-13j)z^3 + (103.8 + -1.421e-14j)z^2 + (-15.09 + 2.665e-15j)z + (1.163 + -2.776e-16j)}
#' \end{equation}