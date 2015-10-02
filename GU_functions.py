# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 06:33:20 2013

GU functions for the Gaba/Glutamate Uncaging (GU) experiment script.

@aut
"""
import sys
import pyplot as plt
import numpy as np
from scipy.ndimage import label
from matplotlib.patches import Circle
from time import time, sleep
from pylibftdi import BitBangDevice
from multiprocessing import Process
from tmca_control import init_connection, experiment_default_settings, move

prestm_col = [.7,.7,.7]


def test_stm_locs():
    
    X_lores, Y_lores, stm_r_lores, ax = prep_stm_locs_lores()
    
    stm_r_hires = 133/2.0
    X_lores, Y_lores = np.array(X_lores), np.array(Y_lores)
    X_lores, Y_lores = X_lores[[1,2,10]], Y_lores[[1,2,10]]
    prep_stm_locs_hires( ax, X_lores, Y_lores, stm_r_lores, stm_r_hires )


def scan(XY, Circs, f, stm_d, ax = None, min_n_per_pos = 3, stm_num = 0):
    """
    Parameters
    ----------
    XY      : Stimulation positions, recaray with "X" and "Y"
    Circs   : A list of "Circle" objects returned from "get_scanpos_lores" or
              "get_scanpos_hires"
    stm_d   : Diameter of UV-stimulation spot
    ax      : Plotting. Axes to redraw/change Circs on
    min_n_per_pos
    stm_num : Counter for stimuli. Number where to start 
    
    Returns
    -------
    
    """
          
    t0 = time()
    ser = init_connection(port = 'COM5', verbose=False)
    experiment_default_settings(ser, verbose=False)

    # minimum number of light stimuli per stimulated position
    # inter stimulus interval, in seconds
    isi = 8.0
       
    # y: + -> door
    # x: + -> right
    tot_x, tot_y = 0, 0
    pos_x, pos_y = 0, 0

    X_roi, Y_roi = [],[]
    Circs_roi = []
    STM_NUMS = []
    
    # Counter for the number of moves
    s, n = '', 0
    N = XY.shape[0]
    RESPONSE = ['xxx']*N
    while (n < N) and (s.lower() != 'q'):
        
        t = time() - t0
                
        if n > 0: 
            dx, dy = XY['X'][n] - XY['X'][n-1], XY['Y'][n] - XY['Y'][n-1]
        else: 
            dx, dy = XY['X'][n], XY['Y'][n]
        
        # Keep moving or quit
        s = 'x'
        while s and (s.isalpha() and not (s.lower() in ['m','q',''])): 
            s = raw_input(' Move or quit ([m],q)?\n')
        if (not s) or (s.lower() == 'm'):
            move(ser, dx, dy)
            
            # Update x,y postion
            pos_x, pos_y = XY['X'][n],XY['Y'][n]
            # Update total distance moved.
            tot_x += abs(dx)
            tot_y += abs(dy)
            Circs[n].set_ec('w')
            draw(); sleep(0.5)         
            print '\n Pos:', n, '(' + str(n+1) + ' of ' + str(N) + ')'
            print ' Position (micron) x:', int(round(pos_x)),', y:', int(round(pos_y))
            print ' Dist moved (micron) dx:', int(round(dx)), ', dy:', int(round(dy)),'\n'

            # Stimulation loop
            n_stm = min_n_per_pos
            STM_NUMS = []
            while n_stm != 0:
                sys.stdout.write(' Trace: ')
                while n_stm > 0:
                    sys.stdout.write(str(n_stm)+', ')
                    STM_NUMS.append(stm_num)
                    send_ttl(10) # Stimulate, send 10 ms ttl
                    stm_num += 1
                    n_stm -= 1
                    sleep(isi)
                s2 = raw_input('\n Num add stim ([0],1,2 ...)?\n')
                while not (s2.isdigit() or s2 == ''):
                    print ' Not valid input'
                    s2 = raw_input(' Num add stim ([0],1,2 ...)?\n')
                if s2:
                    n_stm = int(s2)
                    print ' Num add stim: ', n_stm
                else: print ' No add stim'

            # Register if PSPs or APs were evoked at the stimulation location
            s3 = 'x'
            while not ((s3.isalpha() and (s3.lower() in ['a','p','n'])) or (not s3)):
                s3 = raw_input(' AP, PSP, or Nothing evoked by UV-stimulation (a/p/[n])?\n')
            
            # Update plot
            if s3 == 'p':
                RESPONSE[n] = 'psp'
                circ_col = 'r'
                Circs_roi.append(Circs[n])
                X_roi.append( pos_x )
                Y_roi.append( pos_y )
            elif s3 == 'a':
                RESPONSE[n] = 'dct'
                circ_col = 'y'
                Circs_roi.append(Circs[n])
                X_roi.append(pos_x)
                Y_roi.append(pos_y)
            else: 
                RESPONSE[n] = 'non'
                circ_col = 'b'
            Circs[n].set_ec(circ_col)
            Circs[n].set_ls('solid')
            draw()
            sleep(1e-3)
            
            print ' Evoked response:', RESPONSE[n]
            
            # Write to log file
            for snum in STM_NUMS:
                f.write('%(num)d,%(p_x)d,%(p_y)d,%(psp)d,%(dct)d,%(stm_d)f,%(dx)d,%(dy)d,%(tot_x)d,%(tot_y)d,%(t)f,%(snum)d\n' % \
                        {"num":n,"p_x":pos_x,"p_y":pos_y,"psp":int(RESPONSE[n]=='psp'),
                         "dct":int(RESPONSE[n]=='dct'),"stm_d":stm_d,"dx":dx,
                         "dy":dy,"tot_x":tot_x,"tot_y":tot_y,"t":t,"snum":snum})

            # Count stimulation position                         
            n += 1               
    
    print '\n-----------------------------------'
    print ' Stimulation diameter:',int(round(stm_d)),'microns'
    print ' Total number of moves:', n
    print ' Last stimulus postion (micron) x:', pos_x, ', y:', pos_y
    print ' Total distance moved (micron) x:', tot_x, ', y:', tot_y,'\n\n'
    
    XY_roi = np.recarray(len(X_roi),dtype=[('X',float),('Y',float)])
    XY_roi['X'],XY_roi['Y'] = X_roi,Y_roi
    # Return to zero
    move( ser, -pos_x, -pos_y )
    ser.close()
    
    if len(STM_NUMS): last_stm_num = STM_NUMS[-1]+1
    else: last_stm_num = 0
    
    return XY_roi, Circs_roi, last_stm_num

 
def get_ROI_4_sides(ax = None):
    """
    Returns four sides of an Region Of Interest
    """
    # Assumes scale bar length 500 microns
    scalebar_len_um = 500. # microns
    #prestm_col = [.3,.7,.3]

    # If no axes was set, plot on the currently active.
    if ax is None:
        ax = plt.gca()
        
    print 'Click on the top and the bottom of the scalebar.\n'
    sclbar_xy_px = np.array(plt.ginput(n=2))
    
    print 'Click on pipette tip/cell.\n'
    cell_xy_px = np.array(plt.ginput(n=1, timeout=60))[0]
    
    sclbar_len_px = abs(diff(sclbar_xy_px[:,1]))[0]
    um_to_px = sclbar_len_px/scalebar_len_um
    
    xlim_px, ylim_px = ax.get_xlim(), ax.get_ylim()
    xlim_um = ((xlim_px[0]-cell_xy_px[0])/um_to_px,
               (xlim_px[1]-cell_xy_px[0])/um_to_px)
    ylim_um = ((ylim_px[0]-cell_xy_px[1])/um_to_px,
               (ylim_px[1]-cell_xy_px[1])/um_to_px)
    IM = ax.get_images()[0]
    IM.set_extent((xlim_um[0], xlim_um[1], ylim_um[0], ylim_um[1]))
    xlabel('micron')
    ylabel('micron')
    plt.draw()
    sleep(1e-3)
    
    # Get and plot the outer borders of the area to stimulate
    print 'Click on 4 corners making up the border of region to stimulate.\n'
    region_border_xy = np.array(ginput(n=4))
    lol_bix = np.ones(4, dtype=bool)
    up_ix = np.argsort(region_border_xy[:,1])[2:]
    r_ix = np.argsort(region_border_xy[:,0])[2:]
    lol_bix[up_ix] = False
    lor_bix = lol_bix.copy()
    lol_bix[r_ix] = False
    lor_bix[lol_bix] = False  
    upl_bix = ~(lor_bix | lol_bix)
    upl_bix[r_ix] = False
    upr_bix = ~(lor_bix | lol_bix | upl_bix)
    lo_left = region_border_xy[lol_bix, :].flatten()
    up_right = region_border_xy[upr_bix, :].flatten()
    lo_right = region_border_xy[lor_bix, :].flatten()
    up_left = region_border_xy[upl_bix, :].flatten()
        
    plt.plot([lo_left[0], up_left[0], up_right[0], lo_right[0], lo_left[0]], 
             [lo_left[1], up_left[1], up_right[1], lo_right[1], lo_left[1]],
             c=prestm_col)

    ax.set_xlim(xlim_um)
    ax.set_ylim(ylim_um)
    corners = {'lo_left':lo_left,
               'up_left':up_left,
               'up_right':up_right,
               'lo_right':lo_right}          

    return corners
    

def get_scanpos_manual(stm_r, ax=None):
    """
    """
    # If no axes was set, plot on the currently active.
    tmp_col = [249/255., 77/255., 0/255.]
    if ax is None:
        ax = plt.gca()    

    XY = np.recarray(0, dtype=[('X',float), ('Y',float)])
    Circs = []
    man_add = get_input('Manually add stimulation position', bool, False)
    while man_add:
        retain = False
        while not retain:
            print 'Click on positions to stimulate\n'
            tmp_xy = np.array(ginput(n=1, timeout=60))[0]
            c = Circle(tmp_xy, stm_r, ec=tmp_col, lw=4, fill=False)
            tmp_circ = ax.add_patch(c)
            draw()
            sleep(1e-3)
            retain = get_input('Position ok', bool, True)

            if not retain:
                tmp_circ.remove()
                draw()
                sleep(1e-3)

        XY.resize(XY.shape[0]+1, refcheck=False)
        XY[-1] = tmp_xy[:]
        tmp_circ.set_lw(1)
        tmp_circ.set_ec(prestm_col)
        Circs.append(tmp_circ)
        draw()
        sleep(1e-3)
        man_add = get_input('Add more stimulation positios', bool, False)

    return XY, Circs


def get_scanpos_given_4_sides(stm_r,
                              corners,
                              ax=None,
                              max_pipette_xovershoot=300.):
    """
    """
    # If no axes was set, plot on the currently active.
    if ax is None:
        ax = plt.gca()
        
    step_len = np.floor(np.sqrt(2*stm_r**2))
    flk = step_len/2.0

    lo_l, up_l = corners["lo_left"], corners["up_left"]
    lo_r, up_r = corners["lo_right"], corners["up_right"]
    
    # Clean out x and ys outside outer border of ROI.
    lbord_m = (lo_l[1] - up_l[1]) / (lo_l[0] - up_l[0])
    lbord_b = lo_l[1] - lbord_m*lo_l[0]
    rbord_m = (lo_r[1] - up_r[1]) / (lo_r[0] - up_r[0])
    rbord_b = lo_r[1] - rbord_m*lo_r[0]
    upbord_m = (up_l[1] - up_r[1]) / (up_l[0] - up_r[0])
    upbord_b = up_l[1] - upbord_m*up_l[0]
    lobord_m = (lo_l[1] - lo_r[1]) / (lo_l[0] - lo_r[0])
    lobord_b = lo_l[1] - lobord_m*lo_l[0]
        
    lbord = lambda x, y: x > (y - lbord_b) / lbord_m - flk
    rbord = lambda x, y: x < (y - rbord_b) / rbord_m + flk
    lobord = lambda x, y: y > lobord_m*x + lobord_b - flk
    upbord = lambda x, y: y < upbord_m*x + upbord_b + flk
    
    X_tmp = np.arange(min(lo_l[0],up_l[0]),
                      max(lo_r[0]+100,up_r[0]+100),
                      step_len)
    Y_tmp = np.arange(min(lo_l[1],lo_r[1]),
                      max(up_l[1]+100,up_r[1]+100),
                      step_len)
    # Remove positions outside of max_pipette_xovershoot
    if max_pipette_xovershoot > 0:
        X_tmp = X_tmp[X_tmp <= max_pipette_xovershoot]
    elif max_pipette_xovershoot < 0:
        X_tmp = X_tmp[X_tmp >= max_pipette_xovershoot]

    X,Y = [],[]
    Circs = [];pos_num = 0
    for x in X_tmp:
        for y in Y_tmp:
            if lbord(x,y) and rbord(x,y) and lobord(x,y) and upbord(x,y) :
                c = Circle((x, y), stm_r,ec=prestm_col,
                           ls='dotted', fill=False)
                Circs.append(ax.add_patch(c))
                X.append(x)
                Y.append(y)
                plt.text(x, y, str(pos_num))
                pos_num += 1
    nstmloc = len(X)
    title(r'Stm pos: $N_{'+str(int(round(stm_r*2)))+'}='+str(nstmloc)+'$')
    plt.show()
    plt.draw()
    sleep(1e-3)
    XY = np.recarray(len(X), dtype=[('X',float),('Y',float)])
    XY['X'], XY['Y'] = X, Y
    
    return XY, Circs
    
    
def get_scanpos_given_circles(XY_roi, XY, Circs, stm_r_roi, stm_r,
                              ax=None, max_pipette_xovershoot=300.):

    """
    Parameters
    ----------
    XY_roi                 - center positions of circular regions of interest
    XY                     - 
    Circs                  - 
    stm_r_roi              - radius of circular regions of interest, in microns
    stm_r                  - radius of light stimulus, in microns
    ax                     - matplotlib axes object on which IM is plotted. 
    max_pipette_xovershoot - In microns. Max x distance to move to the right 
                             of pipette tip. Go to far and the objective will 
                             touch the pipette and destroy the patch. 
                             Assumes that the pipette enters from the right.
                             In microns.
    Returns
    -------
                         
    """
    
    # If no axes was set, plot on the currently active.
    if ax is None:
        ax = plt.gca()
    step_len = np.floor(sqrt(2*stm_r**2))
    xlim, ylim = np.sort(ax.get_xlim()), np.sort(ax.get_ylim())
    # Draw a grid with spacing 'step_len' covering the entire image
    y, x = np.ogrid[ylim[0]:ylim[1]:step_len,xlim[0]:xlim[1]:step_len]
    mask = np.zeros((y.shape[0],x.shape[1]))
    for xcntr,ycntr in XY_roi:
        mask = logical_or(mask,
                          ((y - ycntr)**2 + (x - xcntr)**2) <= stm_r_roi**2) 
    labeled_mask, n_labels = label( mask )
    mask_labels = np.unique( labeled_mask )[1:]
    ix0 = XY.shape[0]
    XY.resize(ix0+mask.sum(), refcheck=False)
    for ix, ml in enumerate(mask_labels):
        y_ix, x_ix = nonzero(labeled_mask == ml)
        ix1 = x_ix.shape[0] + ix0
        XY[ix0:ix1] = zip(x[0,x_ix],y[y_ix,0])
        ix0 = ix1
    # Remove positions outside of max_pipette_xovershoot
    if max_pipette_xovershoot > 0:
        XY = XY[XY['X'] <= max_pipette_xovershoot]
    elif max_pipette_xovershoot < 0:
        XY = XY[XY['X'] >= max_pipette_xovershoot]

    for x,y in XY:
        c = Circle((x,y), stm_r, ec=prestm_col, ls='dotted', fill=False)
        Circs.append(ax.add_patch(c)) 
    plt.title(ax.get_title()+r', $N_{'+str(int(round(stm_r*2)))+'}='+str(len(XY))+'$')
    plt.show()
    plt.draw()
    sleep(1e-3)
              
    return XY, Circs
    
    
def get_input(q_str, ans_type, default=None, check=False):
    
    """   
    
    Parameters
    ----------
    q_str       -- Question/input string, e.g.:
                   ' Are the stimulus positions OK'
                   ' 1st scan --- Enter diameter of the UV stimulation spot (micron)'
    ans_type    -- Type of the returned answer, e.g. int, float, str, bool
    default     -- Default value for output, if given it should have the same 
                   type as ans_type. For not defalut value enter None.
    check       -- Whether or not to ask if the input value is ok?
    
    Return
    ------
    ans         -- Reply to asked question, same type as ans_type

    Example
    -------
    ans = get_input(' Enter duration of UV stimulation (ms)', int, 50, True)
            
    """
  
    q_str = ' ' + q_str
    
    if type(default) in [float,int,str]:
        q_str += ': [' + str(default) + ']\n'
    elif type(default) is bool:
        if default:
            q_str += ' ([y]/n)?\n'
        else:
            q_str += ' (y/[n])?\n'
    elif default is None:
        q_str += ':\n'
        
    OK = False
    while not OK:
        
        s = raw_input(q_str)
        
        if not s:
            ans = default
            OK = True
        else:
            if type(default) is bool:
                if s.lower() == 'n': ans = False; OK = True
                elif s.lower() == 'y': ans = True; OK = True
                else: OK = False
            else:                
                try: ans = ans_type(s); OK = True
                except: print ' Invalid input'; OK = False
            
        if check:
            ok = 'x'
            while not ok.lower() in ['y','n','']:
                ok = raw_input(str(ans) + ' OK ([y]/n)?\n')
            if ok.lower() == 'n': OK = False
            
    return ans


def send_ttl(durs=100, ittli=0, blocking=False):
    """
    Sends a ttl over the usb port.
    
    NOTE
    TO DO
    
    Parameters
    ----------
    durs     - 
    ittli    -
    blocking -
    
    Returns
    -------
    Nothing
    
    Either sequences separated by ittli ms, or one ttl pulse.
    Sends a sequence of ttls durations over the usb port.
    Pauses ittli (inter-ttl-interval) between the ttls.
    
    When blocking=False:    
    Uses multiprocessing.Process to not block the process/thread from where it
    is called. It will spawn a new process that can run on another cpu.
    
    Reference:
        https://docs.python.org/3.4/library/multiprocessing.html
        http://hackaday.com/2009/09/22/introduction-to-ftdi-bitbang-mode/

    Hjalmar K. Turesson
    2014-10-13
    modified to handle sequences of ttls on 2015-01-27
    """
    # NOTE TODO:
    # write proper help & check state of bitbang before writing to it.
        
    try:
        ittli = float(ittli)
    except TypeError:
        print('TypeError: ittli needs to be a number.')
        return
    
    def ttl_out(durs, ittli):
        """
        """
        with BitBangDevice() as bb:
            for dur in durs:        
                bb.write('0x02')
                sleep(dur/1000.0)
                bb.write('0x08')
                sleep(ittli/1000.0)
                
    if not hasattr(durs, '__iter__'):
        durs = [durs]
    
    if blocking:
        ttl_out(durs, ittli)
    elif not blocking:
        p = Process(target=ttl_out, args=(durs,ittli))
        p.start()    
    