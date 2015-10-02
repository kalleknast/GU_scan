# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 06:33:20 2013

Gaba/Glutamate Uncaging (GU) experiment script.
GU functions are located in GU_functions.

This script will scan an region of interest at increasing resolution in 3 steps

@author: Hjalmar K Turesson
"""
from pylab import *
from os import path
from GU_functions import get_ROI_4_sides, get_scanpos_given_4_sides
from GU_functions import get_scanpos_given_circles, get_input, scan
from GU_functions import get_scanpos_manual
import tkFileDialog
from Tkinter import Tk
from time import strftime, localtime, sleep

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Set defaults
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set the default path for log files and slice photos here
default_path = '/some/path/to/a/great/place'
# Set default UV stimulation duration in milliseconds
default_stm_dur = 50
# Set default diameter of UV stimulation spot in microns for 1st scan
default_stm_d1 = 500.
# Set default diameter of UV stimulation spot in microns for 2nd scan
default_stm_d2 = 250.
# Set default diameter of UV stimulation spot in microns for 3rd scan
default_stm_d3 = 133.
# Set the minumum number of stimulations per positon for 1st scan
min_n_per_pos1 = 3
# Set the minumum number of stimulations per positon for 2nd scan
min_n_per_pos2 = 3
# Set the minumum number of stimulations per positon for 3rd scan
min_n_per_pos3 = 3
# Maximum distance the microscope objective is allowed to pass to the right 
# of the pipette tip. To avoid the objective touching the pipette.
# Positive values for a pipette entering from right, negative for pipette
# entering from left. In microns
max_pipette_xovershoot = 300.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Set up and open log file 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if not default_path: default_path = path.expanduser('~')
# Path name separator. Probably path.sep on Linux and path.altsep on win
if not path.altsep: psep = path.sep
else: psep = path.altsep
# Get filename of log file. Extension will be .log
# The filename (minus extension) should be the same as the main part of the 
# corresponding data filename (.abf) minus extension and 
# the automatic numbering done by Clampex (000,001,...)
Tk().withdraw()
log_fn = tkFileDialog.asksaveasfilename(defaultextension='log',
                                        initialdir=default_path,
                                        title='Enter name common to both log and data file')
f = open(log_fn, 'w')
f.write('Photo stimulation, GU log\n')
f.write(strftime("%Y-%m-%d, %H:%M:%S\n", localtime()))
f.write('Data_filename: ' + log_fn.split(psep)[-1].strip('.log') + '_xxxx.abf\n')

# Stimulation duration (milliseconds)
stm_dur = get_input('Enter duration of UV stimulation (ms)',int,default_stm_dur)
print ' stm_dur set to:',stm_dur,'ms'
f.write('stm_dur:'+str(stm_dur)+' ms\n')

# Restion membrane potential (mV)
V_rest = get_input('Enter resting membrane potential (mV)', int, check=True)
print ' V_rest set to:',V_rest
f.write('V_rest:'+str(V_rest)+' mV\n')

# Spontaneously active cell
spont_active = get_input('Is the cell spontaneously active', bool, False)
print ' spont_active set to:', spont_active
f.write('spont_active:'+str(int(spont_active))+'\n')

# Subject age (days)
subj_age = get_input('Enter subject (rat) age (days)', int, check=True)
print ' subj_age set to:', subj_age
f.write('subj_age:'+str(subj_age)+' days\n')

# Subject ID
subj_ID = get_input('Enter subject ID (string)', str, check=True)
print ' subj_ID set to:', subj_ID
f.write('subj_ID:'+str(subj_ID)+' \n\n')

# Column headers
f.write('POS_NUM,POS_X,POS_Y,PSP,DIRECT,STM_DIAM,DELTA_X,DELTA_Y,'+\
        'TOTAL_X,TOTAL_Y,TIME,STM_NUM\n')
                     
                     
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Preparation of 1st scan
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ask about illumination diameter for 1st scan UV stimulation spot
stm_d1 = get_input('1st scan --- Enter diameter of the UV stimulation spot (micron)',
                   float,default_stm_d1)
print ' 1st scan -- UV stim diameter set to:',stm_d1,'microns'

# Select slice photo file and show it. The slice photo has to be a png
ok = False
fig = figure()
while not ok:
    fig.clf()
    Tk().withdraw()
    im_fn = tkFileDialog.askopenfilename(filetypes=[("Slice photos", ".png")],
                                         initialdir=default_path,
                                         title='Select slice photo to guide scan')
    ax = fig.add_subplot(111)
    ax.imshow(imread(im_fn),origin='upper')
    draw();sleep(1e-3)
    fig.show()
    ok = get_input('Is the correct slice photo selected',bool,True)

# Get and plot scan positions for the 1st scan.
ok = False
while not ok:
    # Clear the axes and redraw slice photo.
    # Meaningless for the 1st iteration of the while loop, but necessary for 
    # the 2nd (and onwards) iteration 
    # (when the experimenter choses to re-select stimulus positions).
    ax.cla()
    ax.imshow(imread(im_fn),origin='upper')
    # Get the 4 corners of region of interest
    corners = get_ROI_4_sides(ax)
    # Get the stimulus postions
    XY1, Circs1 = get_scanpos_given_4_sides(stm_d1/2., corners, ax,
                                            max_pipette_xovershoot = max_pipette_xovershoot)
    draw()
    sleep(1e-3)
    ok = get_input('1st scan -- Are the stimulation positions ok', bool, True)

# Dummy values to be changed later
XY_roi1, XY_roi2, XY_roi3 = [],[],[]


if len(XY1):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Start 1st scan
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    XY_roi1, Circ_roi1, stm_num = scan(XY1, Circs1, f, stm_d1,
                                       ax, min_n_per_pos1)

                                      
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Preparation for 2nd scan
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
scan2 = get_input('Do a 2nd scan of selected regions at higher resolution', 
                  bool, True)
if scan2:
    # Ask about illumination diameter for 2nd scan
    stm_d2 = get_input('2nd scan --- Enter diameter of the UV stimulation spot (micron)',
                       int, default_stm_d2)
    print ' 2nd scan -- UV stim diameter set to:',stm_d2,'microns'
    
    # Manually add stimulation positions
    XY2, Circs2 = get_scanpos_manual(stm_d2/2., ax=ax)
    
    if len(XY_roi1):
        XY2, Circs2 = get_scanpos_given_circles(XY_roi1,
                                                XY2,
                                                Circs2,
                                                stm_d1/2., 
                                                stm_d2/2., 
                                                ax=ax,
                                                max_pipette_xovershoot=max_pipette_xovershoot)
        
    if len(XY2):
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ## Start 2nd scan
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
        XY_roi2, Circ_roi2, stm_num = scan(XY2, Circs2, f, stm_d2, ax,
                                           min_n_per_pos2, stm_num=stm_num)                                           


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Preparation for 3rd scan
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scan3 = get_input('Do a 3rd scan of selected regions at even higher resolution',
                      bool, True)
    if scan3:
        # Ask about illumination diameter for 3rd scan
        stm_d3 = get_input('3rd scan --- Enter diameter of the UV stimulation spot (micron)',
                           int, default_stm_d3)
        print ' 3rd scan -- UV stim diameter set to:', stm_d3, 'microns'
        
        # Manually add stimulation positions
        XY3, Circs3 = get_scanpos_manual(stm_d3/2., ax=ax)
        
        if len(XY_roi2):             
            XY3, Circs3 = get_scanpos_given_circles(XY_roi2,
                                                    XY3,
                                                    Circs3,
                                                    stm_d2/2., 
                                                    stm_d3/2., 
                                                    ax=ax,
                                                    max_pipette_xovershoot=max_pipette_xovershoot)
         
           
        if len(XY3):
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## Start 3rd scan
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
            XY_roi3, Circ_roi3, stm_num = scan(XY3, Circs3, f, stm_d3, ax,
                                               min_n_per_pos3, stm_num=stm_num)                                        


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Finalize
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
f.close()
fig.set_size_inches([16,11.35])
fig.savefig(im_fn.strip('.png') +'_scanguide.png')
fig.savefig(im_fn.strip('.png') +'_scanguide.eps')
