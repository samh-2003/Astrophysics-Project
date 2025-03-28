# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 15:13:39 2025

@author: samhe
"""

from scipy.io.wavfile import write
import numpy as np
from math import floor
import matplotlib.pyplot as plt

import numpy as np
from math import floor, log10

def envelope(env, t, T, dB_adjust=0):
    """Calculates the envelope with an optional additive dB adjustment.
    
    Args:
        env: Envelope parameters [attack_ms, decay_ms, sustain_level, release_ms].
        t: Current time in samples.
        T: Total length of the note in samples.
        dB_adjust: Decibels to add (can be positive or negative).
    
    Returns:
        Amplitude after applying the envelope and dB adjustment.
    """
    T0 = 0
    T1 = floor(env[0] * samplerate / 1000.0)

    if env[1] == 0:
        T2 = T1
        T3 = T1
    else:
        T2 = floor((env[1] + env[0]) * samplerate / 1000.0)
        
    if env[2] == 0:
        T3 = T2
        V = 1
    else:
        T3 = floor((T - 3 * env[3]) * samplerate / 1000.0)
        V = env[2]
        
    if t < T1:  # Attack
        amplitude = (t * 1000.0) / (env[0] * samplerate)
    
    elif env[1] > 0 and t < T2:  # Decay
        m = (env[2] - 1.0) / env[1]
        b = 1.0 - m * env[0]
        amplitude = m * t * 1000.0 / samplerate + b
    
    elif env[2] > 0 and t < T3:  # Sustain
        amplitude = env[2]
    
    elif t >= T3:  # Release
        amplitude = V * np.exp(-1 * (t - T3) / (env[3] * samplerate / 1000.0))
    
    else:  # Shouldn't normally reach here
        amplitude = 0
    
    # Apply dB adjustment (additive in log space)
    if amplitude > 0:
        amplitude_dB = 20 * log10(amplitude) + dB_adjust
        return 10 ** (amplitude_dB / 20)
    else:
        return 0  # Avoid log(0)

def plot_envelope(note,instrument) :
    global BPM, samplerate, Amplitude_max
    print(note)
    duration = int(floor(samplerate*note[1]*60/BPM))
    transpose=instrument[2]
    print(duration)
    print(instrument[1])
    print(instrument[1][0]*samplerate/1000, instrument[1][1]*samplerate/1000,instrument[1][2], instrument[1][3]*samplerate/1000)
    n=[0 for i in range(duration)]
    T=duration*1000/samplerate
    print(T)
    for t in range(int(duration)):
        n[t]=(envelope(instrument[1],t,T))
    plt.figure()   
    plt.plot(n,label='envelope')
    plt.show()

    return 0

def create_note(note,intrument) :
    """ Takes an instrument definition and a note and creates an appropriate wav snippet from it. 
        uses global variables BPM, samplerate and Amplitude_max transpose multiplies frequency, 
        e.g. 0.25 transposes down two octaves, 2 will transpose up one octave"""
    
    global BPM, samplerate, Amplitude_max
    print(note)
    duration = (floor(samplerate*note[1]*60/BPM))
    transpose=instrument[2]
    print(duration)
    
    n=[0 for i in range(int(duration))]
    
    for f,a,ph in instrument[0] :
        amplitude=10**(a/10)*Amplitude_max # convert from dB back to a linear scale
        ph*=np.pi/180
        f=f*transpose
        
        for t in range(int(duration)):
            n[t]+=floor(amplitude*envelope(instrument[1],t,duration*1000/samplerate)*note[2]*np.sin(2.*np.pi*f*note[0]*t/samplerate)+ph)
        
    return n

# instrument[0] : list containing frequencies of harmonics (as multiples of the fundamental), intensities (in dB) and phase in degrees
# instrument[1] : contains the envelope: attack time (in ms), decay time (in ms, give 0 if no decay),  sustain (as a fraction of the total volume, 0 if no sustain), decay constant (in ms)
# e.g. inst[0]=[(1,-6,0), (2,-8,0), (3, -9,0), (4, -10,0)]  to give the fundamental, first second and third harmonic at -6, -8, -9 aand -10 dB.
# inst[1]=[20,5,0.8,3]
# inst[2] transpose factor: transpose multiplies frequency, e.g. 0.25 transposes down two octaves, 2 will transpose up one octave

string=[[(0.5,-10,0), (1,-12,90), (2,-14,45), (3.5, -14,10), (4, -30,50), (5,-20,120), (6,-18,20)],
           [20,10,0.8,30],
            1]

drum=[[(0.2413,-8,0), (1,-10,0), (2.328,-24,90), (4.11,-25,45), (6.30, -16,10),
        (1.73, -8,50), (3.91,-10,120), (6.71,-15,20), (10.07,-18,90),
        (7.34, -15,50), (11.4,-19,120)],
           [10,0,0,50],
          0.125]

organ=[[(0.5,-10,0), (1,-12,90), (2,-14,45), (3.5, -14,10), (4, -30,50), (5,-20,120), (6,-18,20)],
           [200,0,0.8,100],
           0.5]

steelDrum = [[(1.0, -31.6, 0), (1.1019, -43.2, 0), (1.3283, -51.5, 0), 
              (1.4717, -54.5, 0), (1.7774, -57.9, 0), (2.0415, -35.4, 0),
              (2.4189, -60.6, 0), (2.5472, -59.0, 0), (2.6453, -52.8, 0), (2.7811, -64.1, 0), 
              (2.8340, -63.1, 0), (2.9623, -32.6, 0), (3.1396, -54.1, 0), (3.5509, -64.2, 0),  
              (3.9396, -53.6, 0), (4.0906, -61.4, 0), (4.1962, -64.3, 0), (4.4113, -66.0, 0), 
              (4.6000, -63.6, 0), (5.0038, -60.5, 0), (5.1245, -50.2, 0), (5.2113, -52.2, 0), 
              (5.3887, -46.9, 0), (5.4717, -51.7, 0), (5.5358, -58.6, 0), (5.6038, -54.0, 0), 
              (5.7132, -51.4, 0), (5.9472, -29.7, 0), (6.0981, -42.6, 0), (6.2340, -57.6, 0), 
              (6.3622, -61.1, 0), (6.5849, -54.2, 0), (6.7132, -61.6, 0), (6.9509, -40.0, 0), 
              (7.1774, -53.2, 0), (7.3132, -63.4, 0), (7.4264, -68.2, 0), (7.5887, -59.4, 0), 
              (7.6415, -64.8, 0), (7.7472, -63.9, 0), (7.8377, -58.2, 0), (7.9434, -59.0, 0), 
              (8.1170, -54.3, 0), (8.2189, -62.4, 0), (8.3623, -53.8, 0), 
              (8.4906, -61.3, 0), (8.5585, -63.2, 0), (8.6415, -48.2, 0), (8.7057, -62.1, 0), 
              (8.8491, -50.1, 0), 
              (9.1774, -52.9, 0), (9.2943, -60.9, 0), (9.3623, -61.5, 0), (9.4302, -63.9, 0), 
              (9.5019, -65.5, 0), (9.6491, -68.6, 0), (9.8604, -59.2, 0), 
              (10.0868, -53.3, 0), (10.2906, -69.4, 0), (10.4226, -69.3, 0), 
              (10.5547, -64.4, 0), (10.6906, -76.9, 0), (10.8113, -73.4, 0), (10.9434, -53.6, 0), 
              (11.1887, -62.0, 0), (11.3396, -62.9, 0), (11.4340, -74.1, 0), (11.5509, -60.2, 0), 
              (11.6264, -56.1, 0), (11.7811, -60.5, 0), (11.8943, -43.0, 0), (12.0491, -60.9, 0), 
              (12.1962, -71.3, 0), (12.2868, -74.8, 0), (12.4642, -69.6, 0), 
              (12.6000, -72.3, 0), (12.7962, -64.6, 0), 
              (12.9736, -48.1, 0), (13.2453, -56.9, 0), (13.4981, -78.7, 0), 
              (13.8151, -81.1, 0), (13.9245, -69.1, 0), (14.0491, -67.0, 0), (14.1736, -52.6, 0), 
              (14.2905, -74.1, 0), (14.4189, -66.9, 0), (14.5434, -71.1, 0), 
              (14.7434, -82.7, 0), (14.8868, -56.1, 0), (15.1019, -65.3, 0), 
              (15.3396, -84.1, 0), (15.4302, -85.4, 0), (15.5962, -79.3, 0), 
              (15.7811, -77.0, 0), (15.8717, -77.5, 0), (15.9962, -76.7, 0), (16.1811, -81.5, 0), 
              (16.2642, -82.7, 0), (16.6264, -86.5, 0), (16.7962, -78.6, 0), 
              (17.0302, -79.3, 0), (17.1660, -78.6, 0), 
              (17.8717, -56.6, 0), (18.4453, -88.2, 0), (18.9245, -75.6, 0), 
              (19.0189, -70.1, 0),  
              (20.1887, -73.9, 0), (20.8491, -82.5, 0), (23.8301, -77.9, 0)], [10, 0, 0, 50], 1, 30]

# Global parameters
BPM = 150
samplerate = 44100
Amplitude_max = 32768

song=[]

# Twinkle, twinkle little star  (frequency in Hz, duration in beats, volume in fraction of max)
Notes = [(440,1,0.8),(440,1,0.8),(659,1,0.8),(659,1,0.8),(740,1,0.8),(740,1,0.8),(659,2,0.8),
         (587,1,0.8),(587,1,0.8),(554,1,0.8),(554,1,0.8),(494,1,0.7),(494,1,0.6),(440,2,0.5)]

instrument=steelDrum

if True : # Change to True if you want to see your envelope plotted 
    plot_envelope((440,1,0.8),instrument)

for note in Notes :
    print("working on ",note)
    note_wav=create_note(note,instrument)  # Get the waveform for a single note
    song=song + note_wav                # and append it into the output file

# Output in the correct format

wav=np.array(song,dtype=np.int16)
write("./example.wav", samplerate, wav)
print("Done")