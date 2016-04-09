#' % Bandpass IIR filter design
#' % Khushhall Chandra Mahajan
#' % 7th April 2016
#' M = 141

import numpy as np
import scipy.signal as sg
import matplotlib.pyplot as plt

#' # Given specifications
M = 141
delta_p = 0.15
delta_s = 0.15
f_sample = 100000

# frequencies in digital domain
omega_s1 = 18.2*1000
omega_p1 = 20.2*1000
omega_p2 = 30.2*1000
omega_s2 = 32.2*1000

#' # Conversion of the given frequency
#' 1. Normalize the given digital frequency
#' 2. Convert digital frequency to analog frequency
#' 3. Bandpass to Lowpass conversion
omega = np.array([omega_s1,omega_p1,omega_p2,omega_s2],dtype='f')

# Normalized digital frequencies
omega_normalized=(omega/f_sample)*2*np.pi
#+ term=True
print omega_normalized
#+ term=False
# Conversion to analog bandpass
omega_analog = np.tan(omega_normalized/2)
#+ term=True
print omega_analog
#+ term=False
# Frequency Transformation to low pass domain
omega_o2= omega_analog[1] * omega_analog[2]
B = omega_analog[2] - omega_analog[1]
#+ term=True
print omega_o2
print B
#+ term=False
omega_l = (( omega_analog**2)-omega_o2)/( B*omega_analog)
omega_lpf = min(-1*omega_l[0],omega_l[3])
W_s = omega_lpf
W_p = omega_l[2]
#+ term=True
#Lowpass analog filter specifications
print W_p
print W_s
#+ term=False

#' # Designing the low pass filter - Chebyshev filter
D1 = (1/(1-delta_p)**2)-1
D2=(1/delta_s**2)-1
epsilon = np.sqrt(D1)
N = np.arccosh(np.sqrt(D2)/np.sqrt(D1))/np.arccosh(W_s/W_p)
N = np.ceil(N)
#+ term=True
print N
print epsilon
#+ term=False

# Finding the filter poles
# Total 2N poles
poles = np.zeros([2*N],dtype='complex64')
B_k = np.arcsinh(1/epsilon)/N
for k in range(2*int(N)):
	A_k = (2*k+1)*np.pi/(2*N)
	poles[k] = np.sin(A_k)*np.sinh(B_k) + (1.j)*np.cos(A_k)*np.cosh(B_k)

left_poles = np.zeros([0],dtype='complex64')
for s in poles:
  if(s.real<= 0):
    left_poles=np.append(left_poles,s)
#' Following are the poles of analog low pass filter:
#+ term=True
print left_poles
#+ term=False
# Plotting poles of analog low pass filter
plt.figure(1)
plt.grid(True)	
plt.scatter(left_poles.real,left_poles.imag,s=50,marker='x')
plt.title('Pole zero plot of Low Pass Filter')
plt.ylabel('Imaginary')
plt.xlabel('Real')
#gain calculation
gain = 1 + 0.j
for k in range(4):
	gain =  gain * W_p
gain = gain / (epsilon * 2**(N-1))
gain = gain.real
#+ term=True
print gain
#+ term=False

#' # Convert back to analog bandpass domain
#-------------------------------------------------------------
numerator = gain
denominator = 1+0.j
for p in left_poles:
    numerator = np.poly1d([B,0],r=0)*numerator 
    denominator = np.poly1d([1,-p*B, omega_o2],r=0)*denominator
z,p,k = sg.tf2zpk(numerator,denominator)

plt.figure(2)
plt.grid(True)
plt.scatter(p.real,p.imag,s=50,c='b',marker='x')
plt.scatter(z.real,z.imag,s=50,c='b',marker='o')
plt.title('Pole Zero plot of Analog Bandpass filter')
plt.ylabel('Imaginary')
plt.xlabel('Real')
#-------------------------------------------------------------
#' # Convert back to digital domain
#-------------------------------------------------------------
numerator = gain
denominator = 1+0.j
for p in left_poles:
    numerator = np.poly1d([B, 0, -B],r=0)*numerator 
    denominator = np.poly1d([(1 - (p*B) + omega_o2 ),(-2 + (2*omega_o2)),(omega_o2 + (B*p) + 1)],r=0)*denominator
z,p,k=sg.tf2zpk(numerator,denominator)

plt.figure(3)
plt.grid(True)
plt.scatter(p.real,p.imag,s=50,marker='x')
plt.scatter(z.real,z.imag,s=50,marker='o')
plt.title('Pole Zero plot of Digital Bandpass filter')
plt.ylabel('Imaginary')
plt.xlabel('Real')

# Plot Frequency response
plt.figure(4)
plt.clf()
plt.grid(True)
w,h= sg.freqz(numerator,denominator)
plt.plot((w/np.pi)*f_sample/2, np.absolute(h), linewidth=2)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain')
plt.title('Frequency Response')
plt.show()
#-------------------------------------------------------------
#' # The resulting transfer function:
#' \begin{equation}
#' H(z) = \frac{ 0.03793z^{8}  - 0.1517z^{6}  + 0.2276z^{4}  - 0.1517z^{2} + 0.03793}{24.26z^{8} + 2.283z^{7}  + 75.83z^{6} - 5.559z^{5} + 96.85z^{4} + 4.891z^{3}  + 58.96z^{2} + 1.538z^{1}  + 14.39}
#' \end{equation}
