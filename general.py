import numpy as np

def ztov(zc, zo):
    """
    Transform the redshift measurement to velocity.
    
    Parameters:
    ------
    zc: the cosmological redhsift of the source
    zo: the observed (absorber or emitter) redshift 
    
    Return:
    ------
    The returned value is the velocity 'v' of the outflow in units of
    km/s, here the 'v' is assumed to be positive for outflow.
    
    Formula:
    ------
    The formula used to do the calculatios is: (1+zo) = (1+za)*(1+zc),
    where 'za' is the intrinsic absorber redshift in the source 
    restframe and 1+za = sqrt((1-b)/(1+b)) and b = v/c shoud be positive
    for outflow.
    
    reference: Tombesi et al. 2011 "EVIDENCE FOR ULTRA-FAST OUTFLOWS IN RADIO-QUIET ACTIVE GALACTIC NUCLEI. II. ..."
    """
    
    za = (1+zo)/(1+zc)-1.0
    b = (1-(1+za)**2)/(1+(1+za)**2)
    c = 3.0e5 # units: km/s
    v = b*c
    
    return v

def z2d(zc):
    """
    To obtain the distance of object when cosmological redshift is known.
    
    Parameters:
    -----------
    zc: float
        Cosmological redshift
    
    Return:
    -------
    Distance: float
        Distance in unit Mpc
        
    NOTE: This funcation only works for low redshit for now
    
    Reference:
    [1] Hubble's law: https://en.wikipedia.org/wiki/Hubble%27s_law
    [2] Cosmological redshift: https://astronomy.swin.edu.au/cosmos/c/cosmological+redshift
    """
    c=3.0e5 #km/s
    H0=70 #km/s/Mpc, The Hubble constant
    
    v=zc*c #km/s, the recession velocity, here is the approximation that only valid for low z.
    
    D=v/H0 # Hubble's law
    
    return D


def LEDD(mass:float):
    """
    To calculate the Eddington luminosity of a star with a given mass.
    
    Parameter:
    ------
    The mass should be given in the form of solar mass: M_solar = 2.0*10^33g
    
    Return:
    ------
    The luminosity is given in unit of erg/sec.
    
    Reference: http://www-ppl.s.chiba-u.jp/lecture/radiation/node2.html
    """
    
    return mass*1.2e38

def luminosity(flux:float, distance:float):
    """
    Compute the intrinsic luminosity.
    
    Parameters:
    ------
    flux: the observed flux in a certain energy range in unit of erg/cm^2/sec.
    distance: distance to the source in unit of kpc.
    
    Return:
    ------
    The source luminosity in unit of erg/sec.
    """
    # kpc to cm
    dis = distance*3.086e21
    return flux*4*3.1416*dis**2

def Edd_ratio(flux:float, distance:float, mass:float):
    """
    Compute the Eddington fraction.
    
    Parameters:
    ------
    flux: the observed flux in unit: erg/cm^2/sec
    distance: the distance in unit: kpc
    mass: the mass in solar mass
    
    Return:
    ------
    The Eddington fraction.
    """
    return luminosity(flux, distance)/LEDD(mass)

def Rg(mass):
    """
    Compute the gravitational radius with black hole mass.
    
    Parameter
    ------
    mass: mass of the object in the unit of solar mass
    
    Return:
    ------
    The returned raidus is in unit of m.
    """
    
    M_sun = 2.0e30 # kg
    G = 6.67e-11 # m^3 kg^-1 s^-2
    c = 3.0e8 # m/s
    
    return G*(mass*M_sun)/c**2

def divided_uncertainty(x, xerr, y, yerr, version=2):
    """
    Calculate the uncertainty of u = x/y given the value and uncertainties of x and y.
    
    Parameters:
    ------
    x: 1D array
    xerr: 1D array with the same size of x
    y: 1D array with the same size of x
    yerr: 1D array with the same size of x
    version: default 2
        when 1: if u = x/y, du/u = dx/x + dy/y
        when 2: if u = x/y, (du/u)**2 = (dx/x)**2 + (dy/y)**2
    
    Return:
    ------
    u: 1D array (x/y)
    uerr: 1D array, uncertainty of u
    
    Reference:
    ------
    [1] https://www.cnblogs.com/heaventian/archive/2012/11/24/2786241.html
    [2] https://phas.ubc.ca/~oser/p509/Lec_10.pdf
    [3] https://sciencing.com/how-to-calculate-uncertainty-13710219.html
    
    """
    u = x/y
    
    if (version == 1):
        uerr = u * (xerr/x + yerr/y)
    elif (version == 2):
        uerr = u * np.sqrt((xerr/x)**2 + (yerr/y)**2)
    else:
        raise ValueError("\'version\' parameter can only be 1 or 2. Check if there is any mistake!")
    
    return u, uerr

def add_uncertainty(x, xerr, y, yerr):
    
    u = x + y
    uerr = np.sqrt(xerr**2 + yerr**2)
    
    return u, uerr
    


def keV2Angs(ene):
    """
    Convert keV to Angstrom
    
    keV = 12.398521 / Angstrom
    """
    
    return 12.398521/ene

def yyyymmdd2mjd(string):
    
    """
    Covert datetime in yyyymmdd to MJD
    
    Parameters:
    -----------
    string: str or list
    
    Return:
    -------
    MJD: str or list
    
    """
    
    import datetime
    from astropy.time import Time
    
    if type(string)==str:
    
        time_obj = datetime.datetime.strptime(string, '%Y%m%d')
        tmp = Time(time_obj, format='datetime')
        t = str(tmp.mjd)
        
    elif type(string)==list:
        t = []
        
        for ele in string:
            time_obj = datetime.datetime.strptime(ele, '%Y%m%d')
            tmp = Time(time_obj, format='datetime')
            t.append(str(tmp.mjd))
    
    return t

def yyyy_mm_dd2mjd(string):
    
    """
    Covert datetime in yyyymmdd to MJD
    
    Parameters:
    -----------
    string: str or list
    
    Return:
    -------
    MJD: str or list
    
    """
    
    import datetime
    from astropy.time import Time
    
    if type(string)==str:
    
        time_obj = datetime.datetime.strptime(string, '%Y-%m-%d')
        tmp = Time(time_obj, format='datetime')
        t = str(tmp.mjd)
        
    elif type(string)==list:
        t = []
        
        for ele in string:
            time_obj = datetime.datetime.strptime(ele, '%Y-%m-%d')
            tmp = Time(time_obj, format='datetime')
            t.append(str(tmp.mjd))
    
    return t

def mjd2yyyymmdd(mjd_obj):
    
    from astropy.time import Time
    
    if type(mjd_obj)==str:
    
        tmp = Time(mjd_obj, format='mjd')
        m = tmp.iso

        t=''
        for ele in m.split(' ')[0].split('-'):
            t=t+ele
    elif type(mjd_obj)==list:
        t = []
        for obj in mjd_obj:
            tmp = Time(obj, format='mjd')
            m = tmp.iso
            
            t1=''
            for ele in m.split(' ')[0].split('-'):
                t1=t1+ele
                
            t.append(t1)
    return t

def newpar_list(x, y, z):
    for i in range(z):
        string = f"newpar {int(x+i)}=p{int(y+i)}"
        print(string)


def gamma_kT_2_tau(gamma, kTe):
    """
    Calculate the optical depth tau given gamma (photon index)
    and corona temperature kTe (in keV)
    
    See the equation in Sec. 3.1 in 2101.08043 (Nandi et al. 2020 on Ark 120)
    """
    
    electron_rest_ene = 511.0 # keV
    theta_e = kTe/electron_rest_ene
    
    tau = np.sqrt(9.0/4.0 + 3/(theta_e*(gamma+2)*(gamma-1))) - 1.5
    
    return tau

def delete_head(infile, skip_head=3):
    import subprocess as sp

    comm = f"sed -i -e \"1,{skip_head}d\" {infile}"
    sp.check_call(comm, shell=True)
    
    return True

def rewrite_ipl(infile:str):
    """
    Rewrite the output file of iplot in XSPEC. Each data group will be stored in
    the file 'temp-N.txt' (N starts from 0).

    Parameters:
    -----------
    infile: str
        The filename of the input file
    
    Output:
    -------
    Total number of files
    """
    j = 0
    fw = open('temp-0.txt','w+')
    with open(infile,'r') as f:
        for i in f:
            if (i.split()[0]=='NO'):
                fw.close()
                j = j + 1
                outfile=f"temp-{j}.txt"
                fw = open(outfile,'w+')
            else:
                fw.write(i)           
        fw.close()
    delete_head('temp-0.txt', skip_head=3)
    return j+1

def h2v(infile, outfile):
    """
    Rewrite the horizontal txt file to vertical file (columns).
    """
    f = open(infile)
    s = f.readlines()
    length = len(s[0].split(' '))
    width = len(s)
    fw = open(outfile, "w+")

    for i in range(length):
        string = ''
        for j in range(width):
            string = string + s[j].strip('\n').split(' ')[i] + '\t'
        string = string + '\n'
        #print(string)
        fw.write(string)

    #print(s)
    f.close()
    fw.close()

    return True

def list_from_txt(infile):
    f = open(infile)

    lines = f.readlines()
    f.close()

    lines = list(map(lambda s: s.strip(), lines))

    return lines

def list_to_txt(filename, string_list):
    f = open(filename, "w+")

    for line in string_list:
        f.write(line+'\n')

    f.close()

    return True


def hz2keV(hz):
    h = 6.626e-34 # J s
    keV = 1.6e-16 # J
    
    e=h*hz/keV
    
    return e

def keV2hz(e):
    h = 6.626e-34 # J s
    keV = 1.6e-16 # J
    
    hz=e*keV/h
    
    return hz
