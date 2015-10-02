# import the Stimfit input/output module:
import stfio
import numpy as np
import matplotlib.pyplot as plt
# For patches/fields in figures
from matplotlib.patches import Rectangle, Circle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# For saving the list of dicts produced by plt_cGlu_pos
from pickle import dump, load
import csv
from scipy.interpolate import UnivariateSpline
from scipy.stats import linregress
from os import listdir
from os.path import expanduser

from spectral import rmlines, clip_lines

from scipy.ndimage.interpolation import rotate, shift, zoom
from numpy.linalg import norm


def read_cGlu_log(filename):
    """
    Reads log files from the caged glutamate experiment.
    For example 20111205_180824cGlu.txt
    File contents are returned in a dict with keys corresponding to column
    headers in the log file.

    Hjalmar Turesson, 11-12-12, mod 12-02-09
    """

    f = open(filename , 'r')

    data = {}
    data['notes'] = f.readline().strip()
    data['date'] = f.readline().strip()

    f.readline()

    data['data_filename'] = f.readline().replace("Data_filename:","").strip()
    data['photo_stim_duration'] = int(f.readline().replace("Photo_stm_duration:","").replace("ms","").strip())

    f.readline()

    headers = f.readline().strip().split(',')
    n_heads = len(headers)

    for head in headers:
        data[head] = []

    for line in f:
        fields = line.strip().split(',', n_heads - 1)
        for i, number in enumerate( fields ):
            data[headers[i]].append(eval(number))

    return data


def __flip_x_slice(slice_photo, transforms):

    # NOTE: find a faster way of flipping
    transforms.append({'type': 'flip_x', 'x': 0.0,
                       'y': 0.0, 'angle': 0.0})

    return slice_photo[:,-1::-1,:], transforms


def __shift_slice(slice_photo, transforms):

    print(' Click 1st on a point in the slice and then on\n ' +
          'the corresponding point in the template.\n')
    shift_xy = np.array(plt.ginput(n=2))
    shift_xy = np.diff(shift_xy, axis=0)
    shift_x = shift_xy[0][0]
    shift_y = shift_xy[0][1]

    slice_photo = shift(slice_photo, [shift_y, shift_x, 0.0],
                        mode='constant', cval=0.5)

    transforms.append({'type': 'shift', 'x': shift_x,
                       'y': shift_y, 'angle': 0.0})

    return slice_photo, transforms


def __zoom_slice(slice_photo, axis, transforms):

    zoom_x, zoom_y = 1.0, 1.0

    if axis.upper() == 'X' or axis.upper() == 'B':
        print(' Click 1st on 2 points along x in the slice and then on' +
              ' \n corresponding points in the template.\n')
        zoom_x = np.array(plt.ginput(n=4))
        zoom_x = np.diff(zoom_x[:,0])
        zoom_x = zoom_x[-1]/zoom_x[0]

    if axis.upper() == 'Y' or axis.upper() == 'B':
        print(' Click 1st on 2 points along y in the slice and then on' +
              ' \n corresponding points in the template.\n')
        zoom_y = np.array(plt.ginput(n=4))
        zoom_y = np.diff(zoom_y[:,1])
        zoom_y = zoom_y[-1]/zoom_y[0]

    slice_photo = zoom(slice_photo, [zoom_y, zoom_x, 1.0],
                       mode='constant', cval = 0.5)

    transforms.append({'type': 'zoom',
                       'x': zoom_x,
                       'y': zoom_y,
                       'angle': 0.0 })

    return slice_photo, transforms


def __reshape_slice_photo(old_sz, slice_photo):
    """
    """

    old_sz_y, old_sz_x, old_sz_z = old_sz
    new_sz_y, new_sz_x = slice_photo.shape[:2]
    dsz_x = new_sz_x - old_sz_x
    dsz_y = new_sz_y - old_sz_y
    offset_x = abs(dsz_x / 2)
    offset_y = abs(dsz_y / 2)

    if (dsz_x <= 0) and (dsz_y <= 0):
        reshaped_sp = np.ones((old_sz_y, old_sz_x, old_sz_z)) * 0.5
        reshaped_sp[offset_y : offset_y + new_sz_y,
                    offset_x : offset_x + new_sz_x, :] = slice_photo
    elif (dsz_x > 0) and (dsz_y > 0):
        reshaped_sp = slice_photo[offset_y : offset_y + old_sz_y,
                                  offset_x : offset_x + old_sz_x, :]
    elif (dsz_x <= 0) and (dsz_y > 0):
        reshaped_sp = np.ones((old_sz_y, old_sz_x, old_sz_z)) * 0.5
        reshaped_sp[:, offset_x : offset_x + new_sz_x, :] = \
                slice_photo[offset_y : offset_y + old_sz_y, :, :]
    elif (dsz_x > 0) and (dsz_y <= 0):
        reshaped_sp = np.ones((old_sz_y, old_sz_x, old_sz_z)) * 0.5
        reshaped_sp[offset_y : offset_y + new_sz_y, :, :] = \
                slice_photo[:, offset_x : offset_x + old_sz_x, :]

    return reshaped_sp


def __rotate_slice(slice_photo, transforms):

    s = raw_input(' Degrees to rotate (+- 180, or 0 - 360)?\n')
    alpha = float(s)
    slice_photo = rotate(slice_photo, -alpha, reshape=True,
                         mode='constant', cval=0.5)

    transforms.append({'type': 'rotation', 'x': 0.0,
                       'y': 0.0, 'angle': alpha})

    return slice_photo, transforms
    

def transform_xy(xy, transforms):
    """
    Parameters
    ----------
    xy          - A numpy array/ndarray of xy postions w. x in the 1st column
                  and y in the 2nd.
    transforms  - A recarray of transforms w. following fields:
                      'angle', 'type', 'x' & 'y'
    Returns
    -------
    xy          - A numpy array/ndarray of transformed xy postions 
                  w. x in the 1st column and y in the 2nd.
    """
    
    if (type(xy) is tuple) or (type(xy) is list):
        xy = np.array(xy)
        
    for trf in transforms:
        if trf[1]:
            xy = __update_xy(xy, trf)
            if xy[0] == 0 and xy[1] == 0: import pdb;pdb.set_trace()

    return xy    


def __update_xy(previous_xy, transform):
    """
    """

    if transform['type'] == 'flip_x':
        curr_x = - previous_xy[0]
        curr_y = previous_xy[1]
    elif transform['type'] == 'shift':
        curr_x = previous_xy[0] + transform['x']
        curr_y = previous_xy[1] + transform['y']
    elif transform['type'] == 'zoom':
        curr_x = previous_xy[0] * transform['x']
        curr_y = previous_xy[1] * transform['y']
    elif transform['type'] == 'rotation':
        curr_alpha = transform['angle']
        d = norm(previous_xy)
        previous_alpha = np.angle(previous_xy[0] + 1j * previous_xy[1])
        curr_x = d * np.cos(np.deg2rad(curr_alpha) + previous_alpha)
        curr_y = d * np.sin(np.deg2rad(curr_alpha) + previous_alpha)
    else:
        return 0

    return np.array([curr_x, curr_y])


def __update_xy_list(start_xy, transforms):

    previous_xy = start_xy
    if transforms:
        for transform in transforms:
            curr_xy = __update_xy(previous_xy, transform)
            previous_xy = curr_xy
    else:
        curr_xy = start_xy

    return curr_xy


def __undo_last_trans(slice_photo, transforms):

    utrans = transforms[-1]
    if utrans['type'] == 'flip_x':
        slice_photo = slice_photo[:,-1::-1,:]
    elif utrans['type'] == 'shift':
        shift_xyz = [ -utrans['y'],-utrans['x'], 0.0 ]
        slice_photo = shift(slice_photo, shift_xyz, mode='constant', cval=0.5)
    elif utrans['type'] == 'zoom':
        zoom_yxz = [ 1.0 / utrans['y'], 1.0 / utrans['x'], 1.0 ]
        slice_photo = zoom(slice_photo, zoom_yxz, mode='constant', cval=0.5)
    elif utrans['type'] == 'rotation':
        slice_photo = rotate(slice_photo, utrans['angle'],
                             reshape=False, mode='constant', cval=0.5)

    if len(transforms) == 1:
        transforms = []
    else:
        transforms.pop()

    return slice_photo, transforms


def __draw_template(ax, xlim, ylim, template_fn, template_obj=None):
    """
    """

    # Show template fig on top of slice photo
    if template_obj:
        selections = ''
        for ix, tmp in enumerate(template_fn):
            selections += (' ' + tmp.split('/')[-1] +
                           ':\t' + str(ix) + '\n')
        s = raw_input(' Select template file:\n' + selections)
        if s.isdigit() and (int(s) < len(template_fn)):
            # Load selected template fig
            selected_template_fn = template_fn[int(s)]
            template = plt.imread(selected_template_fn)
            template = template[-1::-1,:,:] # Flip y-axis
            template_obj.set_data(template)
        else:
            # NOTE Try recursive
            # It should let the user try to select again, if he/she messed up.
            __draw_template(ax, xlim, ylim, template_fn, template_obj)
    else:
        # Default template fig is 1st in list
        selected_template_fn = template_fn[0]
        template = plt.imread(selected_template_fn)
        template = template[-1::-1,:,:] # Flip y-axis
        template_obj = ax.imshow(template, origin = 'lower')

    template_obj.set_extent((xlim[0], xlim[1], ylim[0], ylim[1]))
    plt.draw()

    return template_obj, selected_template_fn


def match_slice_to_template(slice_photo_fn, log_fn=None, fig=None):

    # Whats left:
    #    TEST!!!

    debug = False

    tmpl_dir = expanduser('~') + '/Pictures/BST_templates'
    template_fn = [tmpl_dir + '/BST_Fig35_Pare.png',
                    tmpl_dir + '/BST_Fig33_Pare.png',
                    tmpl_dir + '/BST_Fig34_Pare.png',
                    tmpl_dir + '/BST_Fig32_Pare.png']

    # Load slice photo
    slice_photo = plt.imread(slice_photo_fn)
    # Flip y-axis
    # arrays are counted from top to bottom, but in the fig we want to count
    # from bottom to top
    slice_photo = slice_photo[-1::-1,:,:]
    if not fig:
        # Draw slice figure, and get its original size
        fig = plt.figure()

    ax = fig.add_subplot(111)
    sp_obj = ax.imshow(slice_photo, origin='lower')
    original_slice_sz = slice_photo.shape
    # Center figure. Find center pixel and place origo there
    xlim = ax.get_xlim()
    xlim -=  np.diff(xlim) / 2.0
    ylim = ax.get_ylim()
    ylim -=  abs(np.diff(ylim)) / 2.0
    sp_obj.set_extent((xlim[0], xlim[1], ylim[0], ylim[1]))

    # Get conversion from physical distance to pixels
    # from scalebar in slice_photos
    ax.axis([-345,-305+100,-50,80])
    plt.draw()
    print 'Click on the top and the bottom of the scalebar.\n'
    sclbar_xy = plt.ginput(n=2)
    ax.axis((xlim[0], xlim[1], ylim[0], ylim[1]))
    sclbar_x_px = (sclbar_xy[0][0] + sclbar_xy[1][0])/2.0
    sclbar_y_px = np.array([sclbar_xy[0][1], sclbar_xy[1][1]])
    sclbar_len_px = abs(np.diff(sclbar_y_px))[0]
    # Assumes scale bar length 500 microns
    um_to_px = sclbar_len_px / 500.0
    # Plot scalebar
    ax.plot([sclbar_x_px-5, sclbar_x_px+5],
            [sclbar_y_px.min(), sclbar_y_px.min()], '-r')
    ax.plot([sclbar_x_px-5, sclbar_x_px+5],
            [sclbar_y_px.min() + sclbar_len_px,
             sclbar_y_px.min() + sclbar_len_px], '-r')
    ax.plot([sclbar_x_px, sclbar_x_px],
            [sclbar_y_px.min(), sclbar_y_px.min() + sclbar_len_px], '-r', lw=2)
    ax.set_title(slice_photo_fn.split('/')[-1].split('.')[0])

    # Get pippette tip/cell location, stimulation sites will be drawn relative
    # to this
    print ' Click on pipette tip/cell.\n'
    cell_xy = np.array(plt.ginput(n=1, timeout=60))[0]
    print cell_xy
    cell_xy_raw = cell_xy.copy()
    if debug:
        ax.plot([cell_xy[0]-5, cell_xy[0]+5],
                [cell_xy[1]-5,cell_xy[1]+5], '-w')
        ax.plot([cell_xy[0]-5, cell_xy[0]+5],
                [cell_xy[1]+5,cell_xy[1]-5], '-w')
        ax.plot([cell_xy[0], 0], [cell_xy[1], 0],'-g')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    # Draw template_fn[0] as a default
    template_obj, selected_template_fn = __draw_template(ax, xlim, ylim, template_fn)

    # Initialize list to record transforms
    transforms = []

    s = ''
    s = raw_input( ' Select transformation\n' +
                   ' ---------------------\n' +
                   '  Flip slice along x-axis:\tf\n' +
                   '  Shift along x and y axes:\ts\n' +
                   '  Zoom x-axis:\t\t\tzx\n' +
                   '  Zoom y-axis:\t\t\tzy\n' +
                   '  Rotate around center:\t\tr\n' +
                   '  Change template:\t\tt\n' +
                   '  Undo last action:\t\tu\n' +
                   '  Quit:\t\t\t\tq\n')

    while s.lower() != 'q':

        trans = True
        if s.lower() == 'f':
            # Flip slice photo to right hemisphere.
            slice_photo, transforms = __flip_x_slice(slice_photo, transforms)
        elif s.lower() == 's':
            # Shift slice photo in x and y to be approximately aligned on top of
            # slice photo.
            slice_photo, transforms = __shift_slice(slice_photo, transforms)
        elif s.lower() == 'zx':
            # Scale/zoom template fig in x to be of approximately of the same
            # size as the BNTS in slice photo.
            slice_photo, transforms = __zoom_slice(slice_photo,
                                                   'x', transforms)
        elif s.lower() == 'zy':
            # Scale/zoom template fig in y to be of approximately of the same
            # size as the BNTS in slice photo.
            slice_photo, transforms = __zoom_slice(slice_photo,
                                                   'y', transforms)
        elif s.lower() == 'r':
            # Rotate template fig to fig slice photo.
            slice_photo, transfroms = __rotate_slice(slice_photo, transforms)
        elif s.lower() == 't':
            # Change template
            template_obj, selected_template_fn = __draw_template(ax,
                                                                 xlim,
                                                                 ylim,
                                                                 template_fn,
                                                                 template_obj)
            trans = False
        elif s.lower() == 'u':
            # Undo last action/transform
            slice_photo, transforms = __undo_last_trans(slice_photo,
                                                        transforms)
        else:
            trans = False

        if trans:
            # Refresh slice fig, but show only the original number of pixels
            reshaped_sp = __reshape_slice_photo(original_slice_sz, slice_photo)
            sp_obj.set_data(reshaped_sp)
            plt.draw()


        s = raw_input('\n flip [f], shift [s], zoom_x [zx],' +
                      ' zoom_y [zy], rotate [r], template [t],' +
                      ' undo [u], quit [q]\n')

    if log_fn:
        stm_xy = __draw_scan_grid(ax, log_fn, transforms, um_to_px, cell_xy)[0]

    cell_xy = __update_xy_list(cell_xy, transforms)

    # Mark the location of recorded cell
    ax.plot([cell_xy[0] - 5, cell_xy[0]+5],
            [cell_xy[1]-5, cell_xy[1]+5], '-w', lw=2)
    ax.plot([cell_xy[0] - 5, cell_xy[0]+5],
            [cell_xy[1]+5, cell_xy[1]-5], '-w', lw=2)
    if debug:
        ax.plot([0, cell_xy[0]], [0, cell_xy[1]], '-g')

    ax.axis((xlim[0], xlim[1], ylim[0], ylim[1]))

    plt.draw()

    if not transforms:
        reshaped_sp = slice_photo


    if log_fn:
        return stm_xy, cell_xy, cell_xy_raw, transforms, \
                selected_template_fn, um_to_px, fig, reshaped_sp
    else:
        return cell_xy, cell_xy_raw, transforms, \
                selected_template_fn, um_to_px, fig, reshaped_sp


def __draw_scan_grid(ax,
                     filename,
                     transforms,
                     um_to_px,
                     cell_xy,
                     response_xy=None,
                     num_on=True,
                     path_on=True,
                     circ_on=False):

    #NOTE: Fix scaling of illum circles, so that they match template

    # radius of illuminated spot in micrometer scaled to px
    illum_radius = 100.0 * um_to_px

    data = read_cGlu_log(filename)

    # convert from microns to pixel dist.
    stm_x = np.array(data['position_x']) * um_to_px
    # The minus sign below is needed because "up" in the slice picture, and
    # thus up in the figure, is towards the wall in the recording rig, and
    # towards the wall is negative y-steps when controllig the stepper motors
    # from tmca_experiment.py.
    stm_y = - np.array(data['position_y']) * um_to_px
    # Places the scan grids zero pos at cell_xy
    stm_x += cell_xy[0]; stm_y += cell_xy[1]

    #  Transform stimulation coordinates to match template
    stm_xy = np.c_[stm_x, stm_y]
    for ix, xy in enumerate(stm_xy):
        stm_xy[ix] = __update_xy_list(xy, transforms)

    # Plot stimulation positions
    path_n = len(stm_xy) - 1
    for ix, xy in enumerate(stm_xy):
        if path_on and ix < path_n:
            ax.plot([xy[0], stm_xy[ix+1,0]], [xy[1], stm_xy[ix+1,1]], c='k')
        if circ_on:
            # Create and draw a circle showing photo stimulation area
            circ = Circle((xy[0], xy[1]), illum_radius, ec='r', fill=False)
            ax.add_patch(circ)
        if num_on:
            ax.text(xy[0], xy[1], str(ix+1))

    # If response coordinates are supplied
    # Scale, shift, transform and then plot them
    if response_xy:
        # convert from microns to pixel dist.
        response_xy[:,0] *= um_to_px
        response_xy[:,1] *= -um_to_px # See above for the magic minus
        # Places the scan grids zero pos at cell_xy
        response_xy[:,0] += cell_xy[0]
        response_xy[:1] += cell_xy[1]

        for ix, xy in enumerate(response_xy):
            xy = __update_xy_list(xy, transforms)
            response_xy[ix] = xy
            # Create and draw a circle showing photo stimulation area
            circ = Circle(xy[0], xy[1], illum_radius, ec='none', fc=[1,0,0])
            ax.add_patch(circ)

    else: response_xy = np.empty(0)

    return stm_xy, response_xy


def fill_cell_position(cell_xy, ar, xlim, ylim, um_to_px,
                       fill_value=1, radius=50.0):

    dxy = [np.diff(xlim)/ar.shape[1], np.diff(ylim)/ar.shape[0]]
    x_index = np.arange(xlim[0], xlim[1], dxy[0], dtype = int)
    y_index = np.arange(ylim[0], ylim[1], dxy[1], dtype = int)
    spot_x, spot_y = __get_template_spot(radius, dxy, um_to_px)

    sp_x = spot_x + int(round(cell_xy[0]))
    sp_y = spot_y + int(round(cell_xy[1]))
    sp_bool_x = np.logical_and(sp_x >= min(x_index), sp_x <= max(x_index))
    sp_bool_y = np.logical_and(sp_y >= min(y_index), sp_y <= max(y_index))

    sp_x = sp_x[np.logical_and(sp_bool_x, sp_bool_y)]
    sp_y = sp_y[np.logical_and(sp_bool_x, sp_bool_y)]

    x = x_index[sp_x]
    y = y_index[sp_y]
    ar[y,x] = fill_value

    return ar


def fill_positions_fancy(pos_xy, ar, xlim, ylim, radius=150.0):
    
    smooth_factor = 15.0 #
    min_radius, max_radius = radius - smooth_factor, radius + smooth_factor
    
    radii = np.linspace(min_radius, max_radius, 15, endpoint=False)
    fill_values = np.cos(np.linspace(0,np.pi,15, endpoint=False)) + 1.0
    fill_values /= fill_values.sum()
    
    for ix, r in enumerate(radii):
        fill_positions(pos_xy, ar, xlim, ylim,
                       fill_value=fill_values[ix], radius=r)
    

def fill_positions(pos_xy, ar, xlim, ylim, fill_value=1.0, radius=150.0):
    """
    """
    nrow, ncol = ar.shape
    if pos_xy.ndim == 1: pos_xy = pos_xy.reshape(1, 2)
    n = pos_xy.shape[0]
    y, x = np.ogrid[-radius: radius, -radius: radius]
    index = x**2 + y**2 <= radius**2
    
    k_x = np.float(ncol) / abs(np.diff(xlim))
    k_y = np.float(nrow) / abs(np.diff(ylim))
    m_x, m_y = - xlim[0], - ylim[1]
    posxy = pos_xy.copy()  
    posxy[:,0] = (pos_xy[:,0] + m_x) * k_x
    posxy[:,1] = - (pos_xy[:,1] + m_y) * k_y

    X0 = (posxy[:,0] - radius).round()
    X1 = (posxy[:,0] + radius).round()
    Y0 = (posxy[:,1] - radius).round()
    Y1 = (posxy[:,1] + radius).round()
    
    rm_bix_x = np.logical_and((X1 > 0), (X0 < ncol))
    rm_bix_y = np.logical_and(((Y1 > 0)), (Y0 < nrow))
    rm_bix = np.logical_and(rm_bix_x,rm_bix_y)
    
    X0, X1 = X0[rm_bix], X1[rm_bix]
    Y0, Y1 = Y0[rm_bix], Y1[rm_bix]
    n = X0.shape[0]
    for ix in range(n):
        y_ix = np.arange(Y0[ix],Y1[ix])
        x_ix = np.arange(X0[ix],X1[ix])
        y_bix = (y_ix >= 0) & (y_ix < nrow)
        x_bix = (x_ix >= 0) & (x_ix < ncol)
        x0, x1 = x_ix[x_bix][0],x_ix[x_bix][-1]+1
        y0, y1 = y_ix[y_bix][0],y_ix[y_bix][-1]+1
        ar[y0:y1,x0:x1][index[y_bix,:][:,x_bix]] += fill_value

    
def __get_template_spot(radius, dxy, um_to_px, npts=256):

    #    illum_radius = 150.0

    circ_x = np.sin(np.linspace(-np.pi,np.pi, npts)) * radius * um_to_px
    circ_y = np.cos(np.linspace(-np.pi,np.pi, npts)) * radius * um_to_px

    x_index = np.arange(min(circ_x)-1.0, max(circ_x)+1.0, dxy[0], dtype=int)
    spot_x, spot_y = np.array([], dtype=int), np.array([], dtype=int)
    # first left side
    for ix, circy in enumerate(circ_y[0:npts/2]):
        tmp_x = x_index[np.logical_and(x_index >= circ_x[ix], x_index < 0.0)]
        tmp_y = np.zeros(tmp_x.shape, dtype=int) + int(circy)
        spot_x = np.r_[spot_x, tmp_x]
        spot_y = np.r_[spot_y, tmp_y]
    spot_x = np.r_[spot_x, -spot_x]
    spot_y = np.r_[spot_y, spot_y]
    # the two right quadrants
    # midlines
    spot_x = np.r_[spot_x, np.zeros(len(np.unique(spot_y)), dtype=int)]
    spot_y = np.r_[spot_y, np.unique(spot_y)]
    spot_x = np.r_[spot_x, np.unique(spot_x)]
    spot_y = np.r_[spot_y, np.zeros(len(np.unique(spot_x)), dtype=int)]

    return spot_x, spot_y


def plot_cGlu_scan_grid(filename, response_pos=None, grid_on=True, num_on=True,
                        path_on=False, sb_on=True, box_on=True, naked=False):
    """
    Plots a grid representing the scanned area, with circles at the stimulated
    positions.

    Parameters
    ----------
    filename     : filename and path to a cGlu log file.
                   E.g. '/home/hjalmar/Data/111205_1/20111205_150456cGlu.txt'
    response_pos : A list or tuple of lists, arrays or tuples of x & y
                   positions where photo stimulation resulted in responses
                   (PSPs). Optional.
    grid_on      : Whether or not to draw the scan grid.
    num_on       : Add position numbers to the circles, so as to show in which
                   order the positions where stimulated.
    path_on      : Whether or not to draw the scan path.
    sb_on        : Scalebar on/off. Plots 500 micrometer horizontal and vertical
                   scalebars.
    box_on       : Whether or not to plot x and y axes and a frame around the grid.

    Hjalmar Turesson, 11-12-12
    """

    # radius of illuminated spot in micrometer
    illum_radius = 75
    # length in micrometer of plotted scalebar
    sblen = 500

    data = read_cGlu_log(filename)

    if naked:
        grid_on = False
        path_on = False
        box_on = True
        sb_on = True

    xpos = data['position_x']
    ypos = data['position_y']
    x_uniq = np.unique(xpos)
    y_uniq = - np.unique(ypos)
    x_step = abs(x_uniq[x_uniq != 0]).min()
    y_step = abs(y_uniq[y_uniq != 0]).min()

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')

    if grid_on:
        xlim = [x_uniq.min(), x_uniq.max()]
        ylim = [y_uniq.min(), y_uniq.max()]
        X = np.c_[x_uniq,x_uniq].T
        Y = np.c_[y_uniq,y_uniq].T
        XLIM = np.tile(xlim,(Y.shape[1],1)).T
        YLIM = np.tile(ylim,(X.shape[1],1)).T
        ax.plot(X, YLIM, color='k')
        ax.plot(XLIM, Y, color='k')

    if path_on:
        for i in range(len(xpos) - 1):
            ax.plot([xpos[i], xpos[i+1]], [-ypos[i], -ypos[i+1]], color='k')

    if len(xpos) == len(ypos):
        for ix in range(len(xpos)):
            # Create and draw a circle showing photo stimulation area
            circ = Circle((xpos[ix],-ypos[ix]), illum_radius,
                           ec=[1, 0, 0], fc='none')
            ax.add_patch(circ)

            if num_on:
                ax.text(xpos[ix], - ypos[ix], str(ix+1))
    else:
        raise ValueError('Something funny with the log file')

    if response_pos:
        xrespos = response_pos[0]
        yrespos = response_pos[1]
        for ix in range(len(xrespos)):
            # Create and draw a circle showing photo stimulation area
            circ = Circle((xrespos[ix],-yrespos[ix]),illum_radius,
                          edgecolor = None, facecolor = [1, 0, 0])
            ax.add_patch(circ)

    if sb_on:
        # draw scalebar
        sb_pos_x = x_uniq.min() - x_step * 0.95
        sb_pos_y = y_uniq.min() - y_step * 0.95
        ax.plot([sb_pos_x, sb_pos_x], [sb_pos_y, sb_pos_y+sblen],
                'b', linewidth=2)
        ax.plot([sb_pos_x, sb_pos_x + sblen], [sb_pos_y, sb_pos_y],
                'b', linewidth=2)
        ax.text(sb_pos_x + x_step * 0.1, sb_pos_y + y_step * 0.1,
                str(sblen) + '$\mu m$')

    if not box_on:
        # remove box around plot
        ax.set_frame_on(True)
        # drop axes
        for loc, spine in ax.spines.iteritems():
            if loc in ['left','bottom']:
                spine.set_position(('outward',10)) # outward by 10 points
            elif loc in ['right','top']:
                spine.set_color('none') # don't draw spine
            else:
                raise ValueError('unknown spine location: %s'%loc)

        # turn off ticks where there is no spine
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        # Maybe remove some more stuff to really clean up the figure.

    ax.set_xticks(x_uniq)
    ax.set_yticks(y_uniq)
    ax.set_yticklabels(- y_uniq)
    ax.set_xlabel(r'x distance ($\mu m$)')
    ax.set_ylabel(r'y distance ($\mu m$)')

    if naked:
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')


    plt.show()

    return ax


def __get_data(fname, trace_ix=None, start_ix=0,
               end_ix=None, get_ts=False, clean_signal=0):
    """
    Reads in some traces.
    trace_ix has to be a list, numpy array or tuple
    """

    rec = stfio.read(fname)
    Fs = 1000.0/rec.dt # dt is in ms

    if not trace_ix:
        trace_ix = range(len(rec[0]))

    if not np.iterable(trace_ix):
        trace_ix = [trace_ix]

    ntraces = len(trace_ix)

    trace_list = [0]*ntraces
    n = np.zeros(ntraces)

    for ix, tr_ix in enumerate(trace_ix):
        trace_list[ix] = rec[0][tr_ix].asarray()[start_ix:end_ix]
        n[ix] = len(trace_list[ix])

    end_ix = n.min()
    # If traces are of unequal lenght, shortent all to the length of the
    # shortest.
    if np.any(end_ix != n ):
        for ix in range(ntraces):
            trace_list[ix] = trace_list[ix][0:end_ix]


    if clean_signal == 1:
        lines = np.array([100.708])
        for ix, trace in enumerate(trace_list):
            trace_list[ix] = rmlines( trace,
                                      Fs,
                                      lines = 100.708,
                                      pval = 0.01,
                                      pad = 2 )[0]
    elif clean_signal == 2:
        lines = np.array([50.38, 100.708])
        for ix, trace in enumerate(trace_list):
            trace_list[ix] = clip_lines( trace,
                                         Fs,
                                         lines = lines )[0]

    # Create the time axis (times of sample points) only if requested.
    if get_ts:
        ts = np.arange(0, end_ix*rec.dt, rec.dt)
        return np.array(trace_list).T, Fs, ts
    else:
        return np.array(trace_list).T, Fs


def plt_cGlu_pos( fname_log, fname_data,
                  save_dir = '/home/hjalmar/Data/tmp_cGlu/' ):
    """
    Hjalmar K Turesson, 14-02-12
    """


    start_t = 750.0
    end_t = 1050.0
    stm_on_t = 800.0

    # alt for early recordings.
    # 111115_1 & 111115_2
    #    start_t = 400
    #    end_t = 750
    #    stm_on_t = 500

    #start_t = 150
    #end_t = 450
    #stm_on_t = 200

    if save_dir[-1] != '/': save_dir += '/'

    xlim = [start_t, end_t]

    log = read_cGlu_log( fname_log )

    if log['data_filename'] != fname_data.split("/")[-1]:
        print 'Data file name is not consistent ' + \
                ' with log file record.'
        print 'log file record:', log['data_filename']
        print 'data file name:', fname_data.split("/")[-1]
        return 0

    stm_dur = log['photo_stim_duration']
    npos = len(log['stimulus_codes'])
    global_events = ['']*(npos+3)
    global_events[0] = {'data_filename': log['data_filename'],
                        'log_filename': fname_log.split("/")[-1],
                        'cell_name': log['data_filename'][:8],
                        'n_PSP': 0,
                        'PSP_index': [],
                        'n_Direct': 0,
                        'Direct_index':[],
                        'n_None': 0,
                        'None_index': [],
                        'PSP_x':[],
                        'PSP_y':[],
                        'PSP_reversal':[]}


    dt = 0.1
    Fs = 1/dt
    start_ix = np.int_( xlim[0] * Fs )
    end_ix = np.int_( xlim[1] * Fs )
    stm_on_ix = int( round(stm_on_t * Fs - start_ix)  )

    fig_size = (8.5,11.0)
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot( 1, 1, 1 )

    for pos_ix, trace_ix_raw in enumerate(log['stimulus_codes']):
        # NOTE: very ugly hack. Best would be to change the log files to be
        # zero based. Have fun
        for ix in range(len(trace_ix_raw)): # to zeros base
            trace_ix_raw[ix] -= 1

        traces_raw, _, ts = __get_data( fname_data,
                                        trace_ix = trace_ix_raw,
                                        start_ix = start_ix,
                                        end_ix = end_ix,
                                        get_ts = True,
                                        clean_signal = 0)
 
        ntraces = len( trace_ix_raw )
        events_state = np.empty(ntraces,dtype=np.bool)
        events_state[:] = True

        traces_clean = rmean( traces_raw, 5 )

        ymin = np.floor(traces_clean.min() / 5.0) * 5.0 - 3
        # Not above - 35 mV, will cut spikes
        ymax = min(np.ceil(traces_clean.max() / 5.0) * 5.0 + 3, 50.0)
        ylim = [ymin, ymax]

        # Create and draw a rectangle showing photo stimulation period
        rect = Rectangle((stm_on_t, ylim[0]), stm_dur, np.diff(ylim), facecolor="#aaaaaa",
                          edgecolor='none')
        ax.add_patch( rect )
        lines = ax.plot( ts, traces_clean )
        ax.set_xlabel('Time (ms)')
        # store original line colors so as to allow reset.
        lines_color = [0]*len(lines)
        for ix, l in enumerate(lines): lines_color[ix] = l.get_color()

        ax.legend( trace_ix_raw, loc=2, labelspacing=0.1, numpoints = 2 )
        ax.set_ylim(ylim)

        events = __detect_events(traces_clean,stm_on_ix,start_ix,
                                 2.0,Fs,events_state,trace_ix_raw)

        arrows, amp_lines = __draw_events(ax, lines, events )

        txt_pos_x = np.diff(xlim)*0.17 + xlim[0]
        txt_pos_y = ylim[1] - np.diff(ylim)*0.01
        txt = global_events[0]['cell_name'] + '\n' + \
              r'$n = ' + str(log['Number'][pos_ix]) + '\, of \,' + str(npos) + \
              '$\n$x =' + str(log['position_x'][pos_ix]) + '\mu m' + \
              '$\n$y =' + str(log['position_y'][pos_ix]) + '\mu m$'
        ax.text(txt_pos_x,txt_pos_y,txt,verticalalignment='top',family='sans-serif')
        ax.set_xlim(xlim)
        plt.draw()

        # Remove bad traces
        s = '' # dummy, to be replaced by command line input
        while s != 'n':
            s = raw_input( ' Traces to exclude?\n Comma separated' +
                           ' traces indeci, "n" for none or "a" for all. \n [n], a, '
                           + str(trace_ix_raw).strip('[]') + ':\n ')
            if not s: s = 'n' # default option

            if s[0] == 'a':
                s = str(trace_ix_raw).strip('[]')
            if s[0].strip(",").isdigit():
                __rm_traces_update_plot(ax,lines,lines_color,
                                        trace_ix_raw,s,events_state)
                events = __detect_events(traces_clean,
                                         stm_on_ix,start_ix,2.0,Fs,
                                         events_state,trace_ix_raw)
                __move_events( events, arrows, amp_lines, events_state )
                plt.draw()


        # Accept/reject/change events
        # 2 clicks per trace, onset and peak/offset
        # No timeout
        # Right clicking cancels last input
        s = 'x' # x is a dummy value
        while s:
            s = raw_input( ' Events\n Accept: [Enter], reject: ' +
                          '"r,tr_num1, tr_num2,..."' +
                          ' or change: "tr_num"?\n ')
            if s:
                if s[0] is 'r': # Reject traces
                    if len(s) == 1:
                        events_state[:] = False
                    elif s.strip("r,")[0].isdigit():
                        trace_rm_evt_ix = __intstring_2_intlist(s[1:])
                        line_rm_evt_ix = [0]*len(trace_rm_evt_ix)
                        for ix, tr_rm_ix in enumerate(trace_rm_evt_ix):
                            # Check that entered index is in list
                            if trace_ix_raw.count(tr_rm_ix):
                                line_rm_evt_ix[ix] = trace_ix_raw.index(tr_rm_ix)
                                events_state[line_rm_evt_ix[ix]] = False

    #                    __rm_events( events, events_state, trace_ix_raw)
                    __move_events( events, arrows, amp_lines, events_state )
                    plt.draw()
                elif s[0].strip(',').isdigit(): # Change on and offsets of events.
                    trace_ch_evt_ix = __intstring_2_intlist(s)
                    line_ch_evt_ix = [0]*len(trace_ch_evt_ix)
                    for ix, tr_ch_ix in enumerate(trace_ch_evt_ix):
                        if trace_ix_raw.count(tr_ch_ix):
                            line_ch_evt_ix[ix] = trace_ix_raw.index(tr_ch_ix)
                            events_state[line_ch_evt_ix[ix]] = True
                    for l_ix in line_ch_evt_ix:
                        # Make the selected trace fat in order to highlight it
                        lines[l_ix].set_linewidth(3)
                        # Keep zorder for reset
                        line_zorder = lines[l_ix].get_zorder()
                        # bring line to front for better visibility
                        lines[l_ix].set_zorder(100)
                        print (' Click in figure on the fat trace to' +
                               ' select event onset (1st) and event' +
                               ' offset/peak (2nd).')
                        while True:
                            try:
                                evt_on, evt_off = plt.ginput(n=2, timeout=-1,
                                                             mouse_pop=None,
                                                             mouse_stop=None)
                                break
                            except:
                                print ' Input failed. Try again.'
                        # Update events, if trace is in evets['trace_ix']
    #                        if events['trace_ix'].count(trace_ix_raw[l_ix]):
    #                        evt_ix = events['trace_ix'].index(trace_ix_raw[l_ix])

                        events['on_time'][l_ix], on_ts_ix = \
                                findnearest( lines[l_ix].get_data()[0], evt_on[0])
                        events['on_amplitude'][l_ix] = \
                                np.median(lines[l_ix].get_data()[1]
                                          [on_ts_ix-5:on_ts_ix+5])
                        events['off_time'][l_ix], off_ts_ix = \
                                findnearest( lines[l_ix].get_data()[0],
                                             evt_off[0])
                        events['off_amplitude'][l_ix] = \
                                np.median(lines[l_ix].get_data()[1]
                                          [off_ts_ix-5:off_ts_ix+5])
                        events['on_latency'][l_ix] = \
                                events['on_time'][l_ix] - stm_on_t
                        events['off_latency'][l_ix] = \
                                events['off_time'][l_ix] - stm_on_t
                        events['delta_amplitude'][l_ix] = \
                                events['off_amplitude'][l_ix] - \
                                events['on_amplitude'][l_ix]
                        events['duration'][l_ix] = \
                                events['off_time'][l_ix] - events['on_time'][l_ix]

                        # Update the plot
                        __move_events( events, arrows, amp_lines, events_state )
                        # Return line to normal width once change is done.
                        lines[l_ix].set_linewidth(1)
                        plt.draw()

        events['x_position'] = log['position_x'][pos_ix]
        events['y_position'] = log['position_y'][pos_ix]
        events['number'] = log['Number'][pos_ix]

        s = 666 # dummy value
        while (int(s) != 0) and (int(s) != 1) and (int(s) != 2):
            s = raw_input(' Which event type?\n' +
                          ' None: [0], PSP: 1, Direct: 2\n ')
            if not s: s = '0'
            if s.isdigit():
                if int(s) is 0:
                    events['type'],events['type_num'] = 'None', 0
                elif int(s) is 1:
                    events['type'], events['type_num'] = 'PSP', 1
                elif int(s) is 2:
                    events['type'], events['type_num'] = 'Direct', 2
            else:
                s = 666
                print ' Input should be 0, 1 or 2.\n Try again.'

        # Make a global count of event types
        if events['type_num'] == 0:
            global_events[0]['n_None'] += 1
            global_events[0]['None_index'].append(pos_ix)

        elif events['type_num'] == 1:
            global_events[0]['n_PSP'] += 1
            global_events[0]['PSP_index'].append(pos_ix)
            global_events[0]['PSP_x'].append(events['x_position'])
            global_events[0]['PSP_y'].append(events['y_position'])
            V_rev,r,p = reversal_potential(events['delta_amplitude'],
                                           events['on_amplitude'],
                                           plot=False)
            global_events[0]['PSP_reversal'].append(V_rev)
        elif events['type_num'] == 2:
            global_events[0]['n_Direct'] += 1
            global_events[0]['Direct_index'].append(pos_ix)

        # Store the data for all illuminaiton positions in global_events
        global_events[pos_ix+1] = __clean_events_dict(events,events_state)

        fig_name = save_dir + log['data_filename'][:-4] + '_' + \
                '{:0>3}'.format(str(events['number'])) + '_' + \
                events['type'].upper().replace('E','').replace('I','').replace('R','')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        plt.draw()
        plt.savefig( fig_name + '.svg', format='svg')
        plt.savefig( fig_name + '.png', format='png')

        ax.cla()

    # Save global_events
    out_fn = save_dir + global_events[0]['data_filename'].strip('.abf') + '.pik'
    output_f = open(out_fn,'wb')
    dump(global_events,output_f,-1)

    return global_events


def __clean_events_dict(events, events_state):
    """
    Clean out all non events
    """
    for item in events.iteritems():
        if (type(item[1]) == list):
            tmp_array = np.array(item[1])
            events[item[0]] = list(tmp_array[events_state])
        elif (type(item[1]) == np.ndarray):
            events[item[0]] = item[1][events_state]

    return events


def __intstring_2_intlist(int_str):
    """
    Takes a string of integers separated by commas
    and returns a list of integers
    """
    int_str = int_str.strip(' [],').replace(' ',',')
    int_list = int_str.replace(',,',',').split(',')
    for ix, int_val in enumerate(int_list):
        int_list[ix] = int(int_val)
    return int_list


def __rm_events( events, events_state, traces_ix_raw):
    """
    Removes events given by events_state from dict events.
    """

    events['delta_amplitude'][np.invert(events_state)] = np.nan
    events['duration'][np.invert(events_state)] = np.nan
    events['off_amplitude'][np.invert(events_state)] = np.nan
    events['off_latency'][np.invert(events_state)] = np.nan
    events['off_time'][np.invert(events_state)] = np.nan
    events['on_amplitude'][np.invert(events_state)] = np.nan
    events['on_latency'][np.invert(events_state)] = np.nan
    events['on_time'][np.invert(events_state)] = np.nan
    events['trace_ix'] = list(np.array(traces_ix_raw)[events_state])


def __draw_events( ax, lines, events ):

    events_on_t = events['on_time']
    events_off_t = events['off_time']
    events_on_amp = events['on_amplitude']
    events_off_amp = events['off_amplitude']
    events_delta_amp = events['delta_amplitude']

    nlines =  len(lines)

        # Draw arrows at the onset of detected events
    arrow_length = 0.07*abs(np.diff(ax.get_ylim()))[0]

    arrows = {'on':[0]*nlines, 'off':[0]*nlines}
    amp_lines = {'on':[0]*nlines, 'off':[0]*nlines, 'delta':[0]*nlines}
    for ix in range(nlines):
        arrows['on'][ix] = ax.arrow( events_on_t[ix],
                                 arrow_length + events_on_amp[ix],
                                 0, -arrow_length,
                                 color = [0.3,0.3,0.3],
                                 head_width = 5,
                                 shape = 'left',
                                 length_includes_head = True,
                                 width = 1,
                                 head_length = 0.4 )

        amp_lines['on'][ix] = ax.plot( [events_on_t[ix]-10,
                                     events_off_t[ix]],
                                    [events_on_amp[ix],
                                     events_on_amp[ix]], 'k')[0]

        arrows['off'][ix] = ax.arrow( events_off_t[ix],
                                  arrow_length + events_off_amp[ix],
                                  0, -arrow_length,
                                  color = [0.3,0.3,0.3],
                                  head_width = 5,
                                  shape = 'right',
                                  length_includes_head = True,
                                  width = 1,
                                  head_length = 0.4 )

        amp_lines['off'][ix] = ax.plot( [events_off_t[ix]-10,
                                      events_off_t[ix]+10],
                                     [events_off_amp[ix],
                                      events_off_amp[ix]], 'k')[0]

        amp_lines['delta'][ix] = ax.plot( [events_off_t[ix],
                                          events_off_t[ix]],
                                          [events_on_amp[ix],
                                           events_on_amp[ix] +
                                           events_delta_amp[ix]],':k')[0]

    return arrows, amp_lines


def __move_events( events, arrows, amp_lines, events_state ):
    """
    Moves arrows and amp_lines to positions given in events,
    and hides arrows and amp_lines associated with rejected traces.
    """

    for ix, e_state in enumerate(events_state):
        if e_state:
            arrows['on'][ix].set_visible(e_state)
            arrows['off'][ix].set_visible(e_state)
            amp_lines['on'][ix].set_visible(e_state)
            amp_lines['off'][ix].set_visible(e_state)
            amp_lines['delta'][ix].set_visible(e_state)
            new_on_t = events['on_time'][ix]
            new_on_amp = events['on_amplitude'][ix]

            new_off_t = events['off_time'][ix]
            new_off_amp = events['off_amplitude'][ix]
            new_delta = events['delta_amplitude'][ix]

            __move_arrow(arrows['on'][ix], new_on_t, new_on_amp)
            __move_arrow(arrows['off'][ix], new_off_t, new_off_amp)
            __move_amp_off_line(amp_lines['off'][ix], new_off_t, new_off_amp )
            __move_amp_delta_line(amp_lines['delta'][ix],new_off_t,new_off_amp,new_delta)
            __move_amp_on_line(amp_lines['on'][ix], new_on_t, new_off_t, new_on_amp)
        else: # Set arrows and lines to invisible
            arrows['on'][ix].set_visible(e_state)
            arrows['off'][ix].set_visible(e_state)
            amp_lines['on'][ix].set_visible(e_state)
            amp_lines['off'][ix].set_visible(e_state)
            amp_lines['delta'][ix].set_visible(e_state)


def __move_amp_off_line( amp_line, new_x, new_y ):
    """
    Moves an amp_off_line from its old position to a new.
    Only for off lines, but not for on and delta lines.
    """
    x_len = np.diff(amp_line.get_data()[0])[0]
    new_x = [new_x-x_len/2.0, new_x+x_len/2.0]
    new_y = [new_y, new_y]
    amp_line.set_data((new_x,new_y))


def __move_amp_on_line( amp_line, new_on_x, new_off_x, new_y ):
    """
    Moves an amp_on_line from its old position to a new.
    Only for on lines, but not for off and delta lines.
    """
    x_len = new_off_x - new_on_x
    new_x = [new_on_x-10, new_on_x+x_len]# NOTE: hack make better if I have time.
    new_y = [new_y, new_y]
    amp_line.set_data((new_x,new_y))


def __move_amp_delta_line( amp_delta_line, new_x, new_y, new_delta):
    """
    Moves an amp_delta_line.
    """
    new_x = [new_x, new_x]
    new_y = [new_y - new_delta, new_y]
    amp_delta_line.set_data((new_x,new_y))


def __move_arrow(ar, new_x, new_y):
    """
    Moves an arrow from its old position to a new.
    """
    offset_x = new_x - ar.get_xy()[0,0]
    offset_y = new_y - min(ar.get_xy()[:,1])
    new_x = ar.get_xy()[:,0] + offset_x
    new_y = ar.get_xy()[:,1] + offset_y
    ar.set_xy(np.c_[new_x,new_y])


def __detect_events(traces_clean,stm_on_ix,
                    start_ix,thresh_nstd,Fs,events_state,trace_ix_raw):
    """
    Only to be called by plt_cGlu_pos
    """
    cross_th_nt = 10 # ms
    cross_th_nix = int(round(cross_th_nt * Fs)) # Fs is in samples per ms, not Hz
    # ms, don't detect events earlier than start_offset_nt after stm onset.
    start_skip_nt = 10
    start_skip_nix = int(round(start_skip_nt * Fs))
    ntraces, lentraces = min(traces_clean.shape), max(traces_clean.shape)
    # Detect events as significant deviations from baseline
    events_on_ix, events_off_ix = np.zeros(ntraces, dtype = int ), np.zeros(ntraces, dtype = int )
    events_on_amp, events_off_amp = np.zeros(ntraces), np.zeros(ntraces)
    len_detect_segment = lentraces - stm_on_ix + 1 - cross_th_nix
    for tr_ix in range(ntraces):
        trace = traces_clean[:,tr_ix]

        if events_state[tr_ix]:
            threshold = trace[:stm_on_ix].std() * thresh_nstd
            trace_norm = abs(trace[stm_on_ix:] - trace[:stm_on_ix].mean())

            for D_ix in np.arange(start_skip_nix,len_detect_segment):
                if np.all(trace_norm[D_ix:D_ix+cross_th_nix] > threshold):
                    events_on_ix[tr_ix] = D_ix
                    events_on_amp[tr_ix] = trace[D_ix+stm_on_ix]
                    events_state[tr_ix] = True
                    evt_trace = smooth( trace_norm[D_ix:],20)
                    for ix in np.arange(len(evt_trace[:-51])):
                        # Find the peak of the event as the first point where
                        # the following 50 points are lower. 50 points should
                        # be 5 ms.
                        if np.all(evt_trace[ix+1:ix+51] < evt_trace[ix]):
                            events_off_ix[tr_ix] = ix + D_ix
                            events_off_amp[tr_ix] = \
                            trace[events_off_ix[tr_ix]+stm_on_ix]
                            break
                    break

    offset_t = int(round((stm_on_ix + start_ix) / Fs))
    events_on_t = (events_on_ix / Fs) + offset_t
    events_off_t = (events_off_ix / Fs) + offset_t
    events_on_t[np.invert(events_state)] = np.nan
    events_off_t[np.invert(events_state)] = np.nan

    events_dict = {'on_time': events_on_t,
                   'on_latency': events_on_t-offset_t,
                   'on_amplitude':events_on_amp,
                   'off_time': events_off_t,
                   'off_latency': events_off_t-offset_t,
                   'off_amplitude':events_off_amp,
                   'delta_amplitude':events_off_amp-events_on_amp,
                   'duration': events_off_t - events_on_t,
                   'trace_ix': trace_ix_raw}
    return events_dict


def __rm_traces_update_plot(ax,lines,lines_color,traces_ix_raw,s,events_state):
    """
    Only to be called by plt_cGlu_pos
    """
    traces_ix_rm = __intstring_2_intlist(s)
    nlines = len(traces_ix_raw)
    for ix, tr_ix in enumerate(traces_ix_raw):
        if traces_ix_rm.count(tr_ix):
            events_state[ix] = False
            lines[ix].set_color([0.8,0.8,0.8])
        else:
            # reset original color
            lines[ix].set_color(lines_color[ix])
            events_state[ix] = True

    ax.legend_ = None
    ax.legend( traces_ix_raw, loc=2, labelspacing=0.1, numpoints = 2 )
    plt.draw()


def rmean( traces, binwidth ):
    """
    Computes a running average of traces. 
    
    Parameters
    ----------
    traces      - a numpy array of traces with time in the first dimension and 
                  different traces in the 2nd dim. No checking is done.
    binwidth    - The number of sample over which to average. In points/samples.

    """
    
    # creates a destination numpy array for smoothed data
    straces = np.empty( traces.shape )

    npts = traces.shape[0]
    # running mean algorithm
    for i in np.arange( npts ):

        if ( npts - i ) > binwidth:
            # append to list the running mean of `binwidth` values
            # np.mean(sweep) calculates the mean of list
            straces[i, :] = traces[i:(binwidth+i), :].mean( axis = 0)

        else:
            # use all remaining points for the average:
            straces[i, :] = traces[i:, :].mean( axis = 0 )

    return straces


def reversal_potential(dV,V,plot=True):
    """
    Linear regression of PSP amplitude and holding potential.
    Reversal potential at the intercept.

    Returns
    -------
    V_rev   : Reversal potential
    r       : correlation coefficient
    pval    : two-sided p-value for a hypothesis test whose H0 is
              that the slope is zero.
    """

    k, m, r, pval = linregress(V,dV)[:4]
    V_rev = -m/k

    if plot:
        plt.figure()
        ax = plt.subplot(1,1,1)
        ax.plot(V,dV,'or')
        ax.plot(V,k*V+m,'k')
        ax.set_xlabel('Membrane potential (mV)')
        ax.set_ylabel('PSP amplitude (mV)')
        txt = r'$V_{reversal} = ' + str(round(V_rev)) + 'mV$\n' \
                '$r_{corr} = ' + str(round(r,2)) + '$\n' \
                '$p = ' + str(round(pval,4)) + '$'

        ylim = ax.get_ylim()
        xlim = ax.get_xlim()
        x_pos = xlim[0] + np.diff(xlim)*0.1
        y_pos = ylim[0] + np.diff(ylim)*0.1
        ax.text(x_pos,y_pos,txt)

    return V_rev, r, pval


def compile_synopsis_dict( data_dir, header ):
    """
    """

    if not header:
        header = ['cell_name',
                 'n_PSP',
                 'PSP_reversal',
                 'PSP_x',
                 'PSP_y',
                 'PSP_index']

    out_dict = {header[0]:[],
                header[1]:[],
                header[2]:[],
                header[3]:[],
                header[4]:[],
                header[5]:[]}

    data_fnames = listdir(data_dir)
    data_fnames.sort()

    for dat_fn in data_fnames:
            #(dat_fn != "120106_2_cGlu_BSTov_0003.pik") and \
            #(dat_fn != "120106_4_cGlu_BSTov_0003.pik") and \
           #(dat_fn != "120107_2_cGlu_BSTov_0003.pik"):
        f = open(data_dir + dat_fn, 'rb')
        events = load(f)
        f.close()
        for ix, head in enumerate(header):
            if head == 'PSP_reversal':
                if len(events[0][head]) > 1:
                    tmp = [0]*len(events[0][head])
                    for tmp_ix, val in enumerate(events[0][head]):
                        tmp[tmp_ix] = round(val,1)
                else:
                    if events[0][head]:
                        tmp = round(events[0][head][0],1)
                    else:
                        tmp = events[0][head]
                out_dict[head].append(tmp)
            else:
                out_dict[head].append(events[0][head])

    return out_dict


def write_synopsis_to_csv( fname, data_dir, save_dir =
                          '/home/hjalmar/Data/tmp_cGlu/' ):

    out_f = open(save_dir + fname + '_synopsis.csv', 'wb')

    header = ['cell_name',
              'n_PSP',
              'PSP_reversal',
              'PSP_x',
              'PSP_y',
              'PSP_index']

    row = [0]*len(header)
    synopse_dict = compile_synopsis_dict( data_dir, header )
    csv_writer = csv.DictWriter(out_f, header)
    csv_writer.writeheader()
    for row_ix, cname in enumerate(synopse_dict['cell_name']):
        for head_ix, head in enumerate(header):
            row[head_ix] = str(synopse_dict[head][row_ix]).strip('[]')
        row_dict = dict(zip(header,row))
        csv_writer.writerow(row_dict)

    out_f.close()


def write_events_to_csv( filename, events,
                        save_dir = '/home/hjalmar/Data/clean_BST_cGlu/' ):
    """
    Write output from plt_cGlu_pos, global_events, to a spreadsheet.
    """


    # Openoffice readable csv dialect
    csv.Dialect.delimiter = ','
    csv.Dialect.doublequote = False
    csv.Dialect.lineterminator = '\r\n'
    csv.Dialect.quotechar = '"'
    csv.Dialect.quoting = 0
    csv.Dialect.quotechar = '"'
    csv.Dialect.skipinitialspace = False

    sheet1_f = open(save_dir + filename + '_summary.csv', 'wb')

    # Write the 1st sheet, with summary values
    sheet1_header = ['cell_name',
                     'n_PSP',
                     'PSP_x',
                     'PSP_y',
                     'PSP_reversal',
                     'data_filename',
                     'log_filename',
                     'PSP_index']

    sheet1_values = [events[0][sheet1_header[0]],
                     events[0][sheet1_header[1]],
                     events[0][sheet1_header[2]],
                     events[0][sheet1_header[3]],
                     events[0][sheet1_header[4]],
                     events[0][sheet1_header[6]],
                     events[0][sheet1_header[7]]]

    sheet1_dict = dict(zip(sheet1_header,sheet1_values))

    csv_writer = csv.DictWriter(sheet1_f, sheet1_header)
    csv_writer.writeheader()
    csv_writer.writerow(sheet1_dict)

    sheet1_f.close()


def plot_connections(fn = '/home/hjalmar/Data/BNST/cGlu_BNST_data.pik'):

    f = open(fn,'rb')
    rec_dict = load( f )
    f.close()

    n_pos = 6

    tmp = {'tot':np.zeros((40,n_pos)),
           'non':np.zeros((40,n_pos)),
           'und':np.zeros((40,n_pos)),
           'inh':np.zeros((40,n_pos)),
           'exc':np.zeros((40,n_pos))}

    tmp['tot'] = rec_dict['stm'] > 0
    tmp['non'] = rec_dict['stm'] == 1
    tmp['und'] = rec_dict['stm'] == 2
    tmp['inh'] = rec_dict['stm'] == 3
    tmp['exc'] = rec_dict['stm'] == 4

    rec_dict['psp_type'] = tmp

    #    psp_types = ['non','inh','exc']
    psp_types = ['inh','exc']
    #    psp_types = ['non']

    fig_size = (7.0,11.5)
    fig = plt.figure(figsize=fig_size)
    
    rec_pos = 1
    ax = fig.add_subplot( 321 , aspect = 'equal' )
    psp_counts = count_psp( rec_dict, n_pos, rec_pos )
    __plot_psp_map_color( ax, psp_types, psp_counts, rec_pos, colorbar_on = False,
                    total_on = True )
    ax.set_title('Anterior Lateral')

    rec_pos = 2
    ax = fig.add_subplot( 322, aspect = 'equal' )
    psp_counts = count_psp( rec_dict, n_pos, rec_pos )
    __plot_psp_map_color( ax, psp_types, psp_counts, rec_pos, colorbar_on = True,
                    total_on = True )
    ax.set_title('Anterior Medial')

    rec_pos = 3
    ax = fig.add_subplot( 323, aspect = 'equal' )
    psp_counts = count_psp( rec_dict, n_pos, rec_pos )
    __plot_psp_map_color( ax, psp_types, psp_counts, rec_pos, colorbar_on = False,
                    total_on = True )
    ax.set_title('Anterior Lateral, close to AC')

    rec_pos = 4
    ax = fig.add_subplot( 324, aspect = 'equal' )
    psp_counts = count_psp( rec_dict, n_pos, rec_pos )
    __plot_psp_map_color( ax, psp_types, psp_counts, rec_pos, colorbar_on = True,
                    total_on = True )
    ax.set_title('Anterior Medial, close to AC')

    rec_pos = 5
    ax = fig.add_subplot( 325, aspect = 'equal' )
    psp_counts = count_psp( rec_dict, n_pos, rec_pos )
    __plot_psp_map_color( ax, psp_types, psp_counts, rec_pos, colorbar_on = False,
                    total_on = True )
    ax.set_title('Ventral, lateral part')

    rec_pos = 6
    ax = fig.add_subplot( 326, aspect = 'equal' )
    psp_counts = count_psp( rec_dict, n_pos, rec_pos )
    __plot_psp_map_color( ax, psp_types, psp_counts, rec_pos, colorbar_on = True,
                    total_on = True )
    ax.set_title('Ventral, medial part')

    plt.show()


    fig_size = (7.0,11.5)
    fig = plt.figure(figsize=fig_size)

    rec_pos = 1
    ax = fig.add_subplot( 321 , aspect = 'equal' )
    psp_counts = count_psp( rec_dict, n_pos, rec_pos )
    __plot_psp_map_bw( ax, psp_types, psp_counts, rec_pos )
    ax.set_title('Anterior Lateral')

    rec_pos = 2
    ax = fig.add_subplot( 322, aspect = 'equal' )
    psp_counts = count_psp( rec_dict, n_pos, rec_pos )
    __plot_psp_map_bw( ax, psp_types, psp_counts, rec_pos )
    ax.set_title('Anterior Medial')

    rec_pos = 3
    ax = fig.add_subplot( 323, aspect = 'equal' )
    psp_counts = count_psp( rec_dict, n_pos, rec_pos )
    __plot_psp_map_bw( ax, psp_types, psp_counts, rec_pos )
    ax.set_title('Anterior Lateral, close to AC')

    rec_pos = 4
    ax = fig.add_subplot( 324, aspect = 'equal' )
    psp_counts = count_psp( rec_dict, n_pos, rec_pos )
    __plot_psp_map_bw( ax, psp_types, psp_counts, rec_pos )
    ax.set_title('Anterior Medial, close to AC')

    rec_pos = 5
    ax = fig.add_subplot( 325, aspect = 'equal' )
    psp_counts = count_psp( rec_dict, n_pos, rec_pos )
    __plot_psp_map_bw( ax, psp_types, psp_counts, rec_pos )
    ax.set_title('Ventral, lateral part')

    rec_pos = 6
    ax = fig.add_subplot( 326, aspect = 'equal' )
    psp_counts = count_psp( rec_dict, n_pos, rec_pos )
    __plot_psp_map_bw( ax, psp_types, psp_counts, rec_pos )
    ax.set_title('Ventral, medial part')

    plt.show()


def __plot_psp_map_bw(ax, psp_types, psp_counts,
                      rec_pos, arrow_on=True, total_on=True):

    if type(psp_types) is not list:
        psp_types = [psp_types]

    n_types = len(psp_types)
    # Plot positions
    xstep = 1.5
    ystep = xstep * n_types
    x = np.array([0,xstep,0,xstep,0,xstep])
    y = np.array([0,0,-ystep,-ystep,-2*ystep,-2*ystep])
    xmin = x.min() - xstep
    ymax = y.max() + ystep


    for p_ix, ptype in enumerate(psp_types):
        fracs, total = count2fractions( psp_counts, ptype )
        # Size of markers, correspond to the fraction of cells that responded
        # to stimulation in a given sector of BNST.
        print ptype
        print 'fracs: ',fracs.ravel()
        print 'total: ',total

        if arrow_on:
            ax.arrow( x[rec_pos-1]-xstep*0.8,
                      y[rec_pos-1]+xstep*0.8,
                      0.5*xstep,
                      -0.5*xstep,
                      length_includes_head = True,
                      head_width = xstep*0.3,
                      head_starts_at_zero = True,
                      head_length = xstep*0.3,
                      color = [0.6, 0.6, 0.6],
                      width = xstep/50.0)

        if total_on:
            for pos_ix, tot in enumerate(total.ravel()):
                ax.text(x[pos_ix],y[pos_ix],int(tot),
                        fontsize='medium',
                        color=[0.9,0,0],
                        weight='semibold')

        if n_types > 1 and p_ix+1 < n_types:
            ax.plot([x.max()+xstep,x.max()+xstep],
                    [y.min()-0.5*ystep,y.max()+0.5*ystep],
                    '--k')
        if n_types == 3:
            fontsize = 'x-small'
        elif n_types == 2:
            fontsize = 'small'
        else:
            fontsize = 'medium'
        txt_x = x.min()+xstep*0.5
        txt_y = y.max()+ystep*0.65
        if ptype is 'non':
            ax.text(txt_x,txt_y,'No PSPs',size=fontsize,ha='center')
        elif ptype is 'inh':
            ax.text(txt_x,txt_y,'Inhib PSPs',size=fontsize,ha='center')
        elif ptype is 'exc':
            ax.text(txt_x,txt_y,'Excit PSPs',size=fontsize,ha='center')

        x += xstep*3

    xmax = x.max() - 2*xstep
    ymin = y.min() - ystep
    ax.set_ylim([ymin, ymax])
    ax.set_xlim([xmin, xmax])
    ax.set_xticks([])
    ax.set_yticks([])


def __plot_psp_map_color(ax,
                         psp_types,
                         psp_counts,
                         rec_pos,
                         colorbar_on=False,
                         arrow_on=True,
                         total_on=False ):

    if type(psp_types) is not list:
        psp_types = [psp_types]

    n_types = len(psp_types)
    # Plot positions
    xstep = 1.5
    ystep = xstep * n_types
    x = np.array([0,xstep,0,xstep,0,xstep,0,0])         # Last two values are just
    y = np.array([0,0,-ystep,-ystep,-2*ystep,-2*ystep,0,0])     # a hack for the colorscale
    xmin = x.min() - xstep
    ymax = y.max() + ystep


    for p_ix, ptype in enumerate(psp_types):
        fracs, total = count2fractions( psp_counts, ptype )
        # Size of markers, correspond to the number of samples in a given
        # sector of BNST.
        area = np.append((np.pi*2*total.ravel())**2,[0,0])/(2*n_types)
        # Color of markers, correspond to the fraction of cells that responded
        # to stimulation in a given sector of BNST
        c =np.append(fracs.ravel(),[0,0.75]) # Append 0, 1 in order to fix colorscale

        s = ax.scatter(x,y,s=area,marker='o',c = c)
        if arrow_on:
            ax.arrow( x[rec_pos-1]-xstep*0.8,
                      y[rec_pos-1]+xstep*0.8,
                      0.5*xstep,
                      -0.5*xstep,
                      length_includes_head = True,
                      head_width = xstep*0.3,
                      head_starts_at_zero = True,
                      head_length = xstep*0.3,
                      color = 'k',
                      width = xstep/50.0)

        if total_on:
            for pos_ix, tot in enumerate(total.ravel()):
                ax.text(x[pos_ix],y[pos_ix],int(tot),
                        fontsize='medium',
                        color=[0.7,0,0],
                        weight='semibold')

        if n_types > 1 and p_ix+1 < n_types:
            ax.plot([x.max()+xstep,x.max()+xstep],
                    [y.min()-0.5*ystep,y.max()+0.5*ystep],
                    '--k')
        if n_types == 3:
            fontsize = 'x-small'
        elif n_types == 2:
            fontsize = 'small'
        else:
            fontsize = 'medium'
        txt_x = x.min()+xstep*0.5
        txt_y = y.max()+ystep*0.65
        if ptype is 'non':
            ax.text(txt_x,txt_y,'No PSPs',size=fontsize,ha='center')
        elif ptype is 'inh':
            ax.text(txt_x,txt_y,'Inhib PSPs',size=fontsize,ha='center')
        elif ptype is 'exc':
            ax.text(txt_x,txt_y,'Excit PSPs',size=fontsize,ha='center')

        x += xstep*3

    xmax = x.max() - 2*xstep
    ymin = y.min() - ystep
    ax.set_ylim([ymin, ymax])
    ax.set_xlim([xmin, xmax])
    ax.set_xticks([])
    ax.set_yticks([])


    if colorbar_on:
        axins = inset_axes(ax,
                           width="5%", # width = 10% of parent_bbox width
                           height="100%", # height : 100%
                            loc=3,
                           bbox_to_anchor=(1.05, 0., 1, 1),
                           bbox_transform=ax.transAxes,
                           borderpad=0 )
        cb = plt.colorbar(s, cax=axins, ticks = [0,0.25,0.5,0.75])
        cb.set_label('Fraction')


def count_psp(rec_dict, n_pos, curr_pos):
    """
    Helper for plot_connetions
    """
    # dict to keep track of PSP types
    psp_counts = {'tot':np.zeros(n_pos),
                  'non':np.zeros(n_pos),
                  'und':np.zeros(n_pos),
                  'inh':np.zeros(n_pos),
                  'exc':np.zeros(n_pos)}
    pos_ix = rec_dict['pos'] == curr_pos
    for psp_type in psp_counts:
        for pos in range(n_pos):
            psp_counts[psp_type][pos] = \
                np.sum(rec_dict['psp_type'][psp_type][pos][pos_ix])

    return psp_counts


def count2fractions(psp_counts, psp_type):
    """
    Helper for __plot_psp_map
    """
    
    totals = np.zeros((3,2))
    totals[0,0] = psp_counts['tot'][0]
    totals[0,1] = psp_counts['tot'][1]
    totals[1,0] = psp_counts['tot'][2]
    totals[1,1] = psp_counts['tot'][3]
    totals[2,0] = psp_counts['tot'][4]
    totals[2,1] = psp_counts['tot'][5]
    
    fractions = np.zeros((3,2))
    fractions[0,0] = psp_counts[psp_type][0]
    fractions[0,1] = psp_counts[psp_type][1]
    fractions[1,0] = psp_counts[psp_type][2]
    fractions[1,1] = psp_counts[psp_type][3]
    fractions[2,0] = psp_counts[psp_type][4]
    fractions[2,1] = psp_counts[psp_type][5]
        
    return fractions/totals, totals


def bst_region(xy, region, dv_lim = -55.0):
    """
    region: 'medial', 'lateral', 'ventral'
    """
    
    k_in, m_in = 2.65, -40.0  
    
    if region is 'medial':
        # Medial external border
        X = np.array([- 150, -130, -110,-65,-10,30,50])
        Y = np.array([dv_lim, 0, 50,100,150,200, 250])
        s = UnivariateSpline(X, Y, s=1)
        y_ex, x_in = s(xy[0]), (xy[1]-m_in)/k_in
        return (xy[1] > dv_lim) and (x_in > xy[0]) and ( y_ex>xy[1] )
    elif region is 'lateral':
        # Lateral external border
        X = [125, 150]
        Y = [250,dv_lim]
        k_ex =(Y[1]-Y[0])/(X[1]-X[0])
        m_ex = 1800
        # inner border
        x_ex, x_in = (xy[1]-m_ex)/k_ex, (xy[1]-m_in)/k_in
        return (xy[1] > dv_lim) and (x_in <= xy[0]) and (xy[0] < x_ex)
    elif region is 'ventral':
        # Ventral external border
        X = [- 150, -140, -50, 50, 120,150]
        Y = [-55, -100, -220, -220, -150, -55]
        s = UnivariateSpline(X, Y, s=1)
        y_ex = s(xy[0])
        return (xy[1] <= dv_lim) and (xy[1] > y_ex)


def smooth(x, window_len=10, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    import numpy as np    
    t = np.linspace(-2,2,0.1)
    x = np.sin(t)+np.random.randn(len(t))*0.1
    y = smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   

    From: cookb_signalsmooth.py, http://scipy.org/Cookbook/SignalSmooth

    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = getattr(np, window)(window_len)
    y = np.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]

    
def findnearest(array, value):
    ix = (np.abs(array-value)).argmin()

    return array[ix], ix
        