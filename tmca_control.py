from numpy import int_
from platform import system
import serial, time, glob
import sys


def disp_ports():

    print "Found ports:"
    for name in list_serial_ports(): print name


def list_serial_ports():
    if system() == 'Windows':
        # Scan for avilable ports
        available = []
        for ix in range(256):
            try:
                s = serial.Serial(ix)
                available.append((ix, s.portstr))
                s.close()
            except serial.SerialException:
                pass
        return available
    elif system() == 'Darwin': # Mac
        return glob.glob('/dev/tty*') + glob.glob('dev/cu*')
    elif system() == 'Linux':
        return glob.glob('/dev/ttyS*') + glob.glob('/dev.ttyUSB*')
            
    

def reply_bytes_to_data(rpl_bytes):

    rpl_data = 256**3 * rpl_bytes[5] + 256**2 * \
            rpl_bytes[4] + 256 * rpl_bytes[3] + rpl_bytes[2]

    if rpl_bytes[5] > 127:
        rpl_data = rpl_data - 256**4 # Handles negative data

    return(rpl_data )


def send( inst, ser ):
   # send instruction
   # inst must be a list of 6 bytes (no error checking)
   for i in range (6):
       ser.write(chr(inst[i]))
   return


def receive(ser):
   # return 6 bytes from the receive buffer
   # there must be 6 bytes to receive (no error checking)
    nbytes = 6
    r = [0]*nbytes
    for i in range(nbytes):
        r[i] = ord(ser.read(1))
    return r


def motor_specs(model = 'NM17C12S-MC6'):

    m_specs = {}

    if model == 'NM17C12S-MC6':
        m_specs['microstep_sz'] = 0.028125     # Default resolution, degrees
        m_specs['max_speed'] = 7910            # Degree/s
        m_specs['min_speed'] = 0.1318          # Degree/s
        m_specs['speed_resolution'] = 0.1318   # Degree/s
        m_specs['motor_steps_per_rev'] = 200
    else:
        print 'Data for the model is lacking'

    return m_specs


def stage_specs():
    """
    Specifications for individual xy-stages
    """
    # TODO: make better

    s_specs = {'um_per_turn': 2000.0/3}

    return s_specs


def microns_to_microsteps(um, model='NM17C12S-MC6'):
    """
    Parameters
    ----------
    um - microns
    model - model of motor connected to T-MCA

    Returns
    -------
    n_microsteps - number of microsteps to move um microns
    """

    # Full turn: 2000.0/3 um
    # step: 100 um -> 100 / (2000.0/3) turns (approx 0.15)
    # step in degrees: 0.15 * 360
    # in microsteps: 0.15 * 360 / m_specs['microstep_sz']

    s_specs = stage_specs()
    m_specs = motor_specs(model)
    n_microsteps = int_(round((um/s_specs['um_per_turn']*360) /
                        m_specs['microstep_sz']))

    return n_microsteps


def experiment_default_settings(ser, model='NM17C12S-MC6', verbose=True):

    m_specs = motor_specs(model)

    target_speed = 2.0                          # Degree/s
    target_speed_bytes = cmd_data_to_bytes(
        int_(round(target_speed * 1.6384 / m_specs['microstep_sz'])))

    if verbose:
        print 'Target speed is set to: ', target_speed, 'degree/s.'
        print 'The command bytes are: ', target_speed_bytes

    instruction = [1, 42]
    instruction.extend(target_speed_bytes)
    if verbose:
        print instruction
    send(instruction, ser)
    instruction[0] = 2
    if verbose:
        print instruction
    send(instruction, ser)

    acceleration = 100                          # Degree/s
    acc_bytes = cmd_data_to_bytes(
        int_(round(acceleration/(10000/1.6384*m_specs['microstep_sz']))))
    if verbose:
        print 'Acceleration is set to: ', acceleration, 'degreee/s^2.'
        print 'The command bytes are: ', acc_bytes

    instruction = [1, 43]
    instruction.extend(acc_bytes)
    if verbose:
        print instruction
    send(instruction, ser)
    instruction[0] = 2
    if verbose:
        print instruction
    send(instruction, ser)

    # TODO: set max relative move, cmd 46


def cmd_data_to_bytes(cmd_data):

    cmd_bytes = [0, 0, 0, 0]

    if cmd_data < 0:                # Handles negative data
        cmd_data = 256**4 + cmd_data

    cmd_bytes[3] = cmd_data / 256**3
    cmd_data = cmd_data - 256**3 * cmd_bytes[3]
    cmd_bytes[2] = cmd_data / 256**2
    cmd_data = cmd_data - 256**2 * cmd_bytes[2]
    cmd_bytes[1] = cmd_data / 256
    cmd_data = cmd_data - 256 * cmd_bytes[1]
    cmd_bytes[0] = cmd_data

    return cmd_bytes 


def init_connection(port='COM5', verbose=True):

    # Open serial port
    try:
        ser = serial.Serial(port, 9600, 8, 'N', 1, timeout=5)
    except:
        print 'Error opening serial port. Quitting.'
        return sys.exit(0)

    if verbose:
        print "Opening " + ser.portstr

    # Renumber devices
    instruction = [0, 2, 0, 0, 0, 0] # All devices renumber (must be done if a
                                     # new device was added or removed from the
                                     # daisy chain)
    if verbose:
        print "Sending renumbering instruction: ", instruction
    send(instruction, ser)

    time.sleep(1)

    return ser 


def move(ser, x_dist=0, y_dist=0, verbose=False):
    """

    Move the stage a relative distance specified by x_dist and y_dist.

    ser     : ABCMeta, open comm port
    x_dist  : distance to move in microns
    y_dist  : distance to move in microns
    verbose : print communication with devices
    TODO:
    set max relative move (cmd 46)

    Parameters
    ----------

    Returns
    -------

    """

    x_step = 0
    y_step = 0

    # If x/y_dist is positive then max_x/y_step will be positive, if x/y_dist
    # is negative then max_x/y_step will be negative.
    if x_dist != 0:
        max_x_step = 500 * x_dist/abs(x_dist)
    if y_dist !=0:
        max_y_step = 500 * y_dist/abs(y_dist)

    while abs(x_dist) > 0:

        if (abs(max_x_step) - abs(x_dist)) < 0:
            x_step = max_x_step
            x_dist = x_dist - max_x_step
        else:
            x_step = x_dist
            x_dist = 0

        instruction = [1, 21]
        # From microns to microsteps
        instruction.extend(cmd_data_to_bytes(microns_to_microsteps(x_step)))
        send(instruction, ser)
        time.sleep(10.0 * x_step/max_x_step)
        if verbose:
            try:
                reply = receive(ser)
                print 'Instruction sent: ', instruction
                print 'Device number: ', reply[0]
                print 'Command number: ', reply[1]
                print 'Supply voltage: ', reply[2]/10.0, 'V'
                print 'Data: ', reply_bytes_to_data( reply )
            except:
                print "No reply was received."

    while abs(y_dist) > 0:

        if (abs(max_y_step) - abs(y_dist)) < 0:
            y_step = max_y_step
            y_dist = y_dist - max_y_step
        else:
            y_step = y_dist
            y_dist = 0

        instruction = [2, 21]
        # From microns to microsteps
        instruction.extend(cmd_data_to_bytes(microns_to_microsteps(y_step)))
        send(instruction, ser)
        time.sleep(10.0 * y_step/max_y_step)
        if verbose:
            try:
                ### command number 255: ERROR!
                reply = receive(ser)
                print 'Instruction sent: ', instruction
                print 'Device number: ', reply[0]
                print 'Command number: ', reply[1]
                print 'Supply voltage: ', reply[2]/10.0, 'V'
                print 'Data: ', reply_bytes_to_data(reply)
            except:
                print "No reply was received."

def solitary_move(x_dist=0, y_dist=0, verbose=False, port='COM5'):
    """

    Initiate port move the stage a distance specified by x_dist and y_dist and
    close port.

    x_dist  : distance to move in microns
    y_dist  : distance to move in microns
    verbose : print communication with devices
    port : { string, defalut: COM1 (linux:/dev/ttyUSB0) }
    TODO:
    set max relative move (cmd 46)

    Parameters
    ----------

    Returns
    -------

    """

    # open serial port
    # replace "/dev/ttyUSB0" with "COM1", "COM2", etc in Windows
    try:
        ser = serial.Serial(port, 9600, 8, 'N', 1, timeout=5)
    except:
        print "Error opening serial port. Quitting."
        return sys.exit(0)

    print "Opened " + ser.portstr

    move(ser, x_dist, y_dist)
    print "Closing " + ser.portstr
    ser.close()

def test(port='COM5', instruction=[1, 52, 0, 0, 0, 0]):
    """
    Tests communication with T-MCA Stepper Motor Controllers

    Parameters
    ----------

    port : { string, defalut: /dev/ttyUSB0 }
    """
    # open serial port
    # replace "/dev/ttyUSB0" with "COM1", "COM2", etc in Windows
    try:
        ser = serial.Serial(port, 9600, 8, 'N', 1, timeout=5)
    except:
        print "Error opening serial port. Quitting."
        return sys.exit(0)

    print "Opening " + ser.portstr

    # Byte 1 - Device number
    # Byte 2 - Command number
    # Byte 3 - Data - Least Significant Byte (LSB)
    # Byte 4 - Data
    # Byte 5 - Data
    # Byte 6 - Data - Most Significant Byte (MSB)

    #instruction = [1,21,20,0,0,0]              # move 1st motor 20 steps
    #instruction = [1,52,0,0,0,0]               # return power supply voltage
    #instruction = [0,2,0,0,0,0]                # All devices renumber (must
                                                # be done if a new device was
                                                # added or removed from the
                                                # daisy chain)
    print "Sending instruction", instruction
    send(instruction, ser)

    time.sleep(1)                                   # wait for 1 second

    try:
        reply = receive(ser)
        print "Receiving reply", reply
        print "Device number:", reply[0]
        print "Command number:", reply[1]
        print "Supply voltage:", reply[2]/10.0, "V"
        print 'Data: ', reply_bytes_to_data(reply)
    except:
        print "No reply was received."

    print "Closing " + ser.portstr
    ser.close()
