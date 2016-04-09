#' % Bandstop FIR filter design
#' % Khushhall Chandra Mahajan
#' % 7th April 2016
#' M = 141

import numpy as np
import scipy.signal as sg
import matplotlib.pyplot as plt
import pylab

#' # Given specifications
M = 141
del_1 = 0.15
del_2 = 0.15
f_sampling = 100000

# frequencies in digital domain
omega_p1 = 9.4*1000.0
omega_s1 = 11.4*1000.0
omega_s2 = 21.4*1000.0
omega_p2 = 23.4*1000.0
transition = 2000.0 *2/f_sampling
f_digital = np.array([omega_p1,omega_s1,omega_s2,omega_p2],dtype='f')

#normalized digital frequencies
f_normalized = (f_digital/f_sampling)*2*np.pi
#+ term=True
print "Normalized:\n",f_normalized
#+ term=False

#' # Kaiser window coefficient 
A = -20*np.log10(del_1)
N, beta = sg.kaiserord(A,transition)
N = (N-1)/2
h_kaiser = sg.kaiser(2*N+1,beta)

#' # h[n] of Ideal BSF
w_cut_1 = (f_normalized[1]+f_normalized[0])*0.5
w_cut_2 = (f_normalized[3]+f_normalized[2])*0.5
iterable = ((np.sin(w_cut_1*k)-np.sin(w_cut_2*k))/(np.pi*k) for k in range(int(-N),int(N+1)))
h_ideal = np.fromiter(iterable,float)
h_ideal[N] = ((w_cut_1-w_cut_2)/np.pi)+1

h_fir = h_ideal*h_kaiser
#+ term=True
print "FIR Filter Coefficients:\n",h_fir
#+ term=False

#' # Plot the FIR filter coefficients
plt.figure(1)
plt.plot(h_fir, 'bs-', linewidth=3)
plt.title('Filter Coefficients (%d taps)' % (2*N+1))
plt.grid(True)
#' # Plot frequency response
plt.figure(2)
plt.clf()
plt.grid(True)
w,h= sg.freqz(h_fir)
plt.plot((w/np.pi)*f_sampling*0.5, np.absolute(h), linewidth=2)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain')
plt.title('Frequency Response')
plt.show()
#' # Plot phase response
plt.figure(3)
plt.grid(True)
# Take tan inverse to find the phase plot
hPhase = pylab.unwrap(np.arctan2(np.imag(h),np.real(h)))
plt.plot(w/max(w),hPhase, linewidth=2)
plt.ylabel('Phase (radians)')
plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
plt.title(r'Phase response')
#http://www.labbookpages.co.uk/audio/firWindowing.html
