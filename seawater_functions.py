def dens0(s, t):
    """
    Density of Sea Water at atmospheric pressure.

	Copied from the Seawater Package: https://pypi.org/project/seawater/

    Parameters
    ----------
    s(p=0) : array_like
             salinity [psu (PSS-78)]
    t(p=0) : array_like
             temperature [℃ (ITS-90)]

    Returns
    -------
    dens0(s, t) : array_like
                  density  [kg m :sup:`3`] of salt water with properties
                  (s, t, p=0) 0 db gauge pressure

    Examples
    --------
    >>> # Data from UNESCO Tech. Paper in Marine Sci. No. 44, p22
    >>> import seawater as sw
    >>> from seawater.library import T90conv
    >>> s = [0, 0, 0, 0, 35, 35, 35, 35]
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> sw.dens0(s, t)
    array([  999.842594  ,   999.842594  ,   995.65113374,   995.65113374,
            1028.10633141,  1028.10633141,  1021.72863949,  1021.72863949])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J. and  Poisson, A. International one-atmosphere
       equation of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
       doi:10.1016/0198-0149(81)90122-9

    """
    import numpy as np
    
    s, t = list(map(np.asanyarray, (s, t)))

    T68 = T68conv(t)

    # UNESCO 1983 Eqn.(13) p17.
    b = (8.24493e-1, -4.0899e-3, 7.6438e-5, -8.2467e-7, 5.3875e-9)
    c = (-5.72466e-3, 1.0227e-4, -1.6546e-6)
    d = 4.8314e-4
    return (smow(t) + (b[0] + (b[1] + (b[2] + (b[3] + b[4] * T68) * T68) *
            T68) * T68) * s + (c[0] + (c[1] + c[2] * T68) * T68) * s *
            s ** 0.5 + d * s ** 2)

def dens(s, t, p):
    """
    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.
    
    Copied from the Seawater Package: https://pypi.org/project/seawater/
    

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].

    Returns
    -------
    dens : array_like
           density  [kg m :sup:`3`]

    Examples
    --------
    >>> # Data from Unesco Tech. Paper in Marine Sci. No. 44, p22.
    >>> import seawater as sw
    >>> from seawater.library import T90conv
    >>> s = [0, 0, 0, 0, 35, 35, 35, 35]
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> p = [0, 10000, 0, 10000, 0, 10000, 0, 10000]
    >>> sw.dens(s, t, p)
    array([  999.842594  ,  1045.33710972,   995.65113374,  1036.03148891,
            1028.10633141,  1070.95838408,  1021.72863949,  1060.55058771])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J., Chen, C.T., Bradshaw, A., and Schleicher, K. A new
       high pressure equation of state for seawater. Deap-Sea Research., 1980,
       Vol27A, pp255-264. doi:10.1016/0198-0149(80)90016-3

    """
    import numpy as np
    
    s, t, p = list(map(np.asanyarray, (s, t, p)))

    # UNESCO 1983. Eqn..7  p.15.
    densP0 = dens0(s, t)
    K = seck(s, t, p)
    p = p / 10.  # Convert from db to atm pressure units.
    return densP0 / (1 - p / K)

def seck(s, t, p=0):
    """
    Secant Bulk Modulus (K) of Sea Water using Equation of state 1980.
    UNESCO polynomial implementation.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].

    Returns
    -------
    k : array_like
        secant bulk modulus [bars]

    Examples
    --------
    >>> # Data from Unesco Tech. Paper in Marine Sci. No. 44, p22.
    >>> import seawater as sw
    >>> s = [0, 0, 0, 0, 35, 35, 35, 35]
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> p = [0, 10000, 0, 10000, 0, 10000, 0, 10000]
    >>> sw.seck(s, t, p)
    array([ 19652.21      ,  22977.2115    ,  22336.0044572 ,  25656.8196222 ,
            21582.27006823,  24991.99729129,  23924.21823158,  27318.32472464])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J. and  Poisson, A. International one-atmosphere equation
       of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
       doi:10.1016/0198-0149(81)90122-9

    """
    import numpy as np
    
    s, t, p = list(map(np.asanyarray, (s, t, p)))

    # Compute compression terms.
    p = p / 10.0  # Convert from db to atmospheric pressure units.
    T68 = T68conv(t)

    # Pure water terms of the secant bulk modulus at atmos pressure.
    # UNESCO Eqn 19 p 18.
    # h0 = -0.1194975
    h = [3.239908, 1.43713e-3, 1.16092e-4, -5.77905e-7]
    AW = h[0] + (h[1] + (h[2] + h[3] * T68) * T68) * T68

    # k0 = 3.47718e-5
    k = [8.50935e-5, -6.12293e-6, 5.2787e-8]
    BW = k[0] + (k[1] + k[2] * T68) * T68

    # e0 = -1930.06
    e = [19652.21, 148.4206, -2.327105, 1.360477e-2, -5.155288e-5]
    KW = e[0] + (e[1] + (e[2] + (e[3] + e[4] * T68) * T68) * T68) * T68

    # Sea water terms of secant bulk modulus at atmos. pressure.
    j0 = 1.91075e-4
    i = [2.2838e-3, -1.0981e-5, -1.6078e-6]
    A = AW + (i[0] + (i[1] + i[2] * T68) * T68 + j0 * s ** 0.5) * s

    m = [-9.9348e-7, 2.0816e-8, 9.1697e-10]
    B = BW + (m[0] + (m[1] + m[2] * T68) * T68) * s  # Eqn 18.

    f = [54.6746, -0.603459, 1.09987e-2, -6.1670e-5]
    g = [7.944e-2, 1.6483e-2, -5.3009e-4]
    K0 = (KW + (f[0] + (f[1] + (f[2] + f[3] * T68) * T68) * T68 +
                (g[0] + (g[1] + g[2] * T68) * T68) * s ** 0.5) * s)  # Eqn 16.
    return K0 + (A + B * p) * p  # Eqn 15.

def T68conv(T90):
    """
    Convert ITS-90 temperature to IPTS-68

    :math:`T68  = T90 * 1.00024`

    Parameters
    ----------
    t : array_like
           temperature [℃ (ITS-90)]

    Returns
    -------
    t : array_like
           temperature [℃ (IPTS-68)]

    Notes
    -----
    The International Practical Temperature Scale of 1968 (IPTS-68) need to be
    correct to the ITS-90. This linear transformation is accurate within
    0.5 ℃ for conversion between IPTS-68 and ITS-90 over the
    oceanographic temperature range.

    Examples
    --------
    >>> import seawater as sw
    >>> T68conv(19.995201151723585)
    20.0

    References
    ----------
    .. [1] Saunders, P. M., 1991: The International Temperature Scale of 1990,
       ITS-90. WOCE Newsletter, No. 10, WOCE International Project Office,
       Southampton, United Kingdom, 10.

    """
    import numpy as np
    
    T90 = np.asanyarray(T90)
    return T90 * 1.00024

def smow(t):
    """
    Density of Standard Mean Ocean Water (Pure Water) using EOS 1980.

    Parameters
    ----------
    t : array_like
        temperature [℃ (ITS-90)]

    Returns
    -------
    dens(t) : array_like
              density  [kg m :sup:`3`]

    Examples
    --------
    >>> # Data from UNESCO Tech. Paper in Marine Sci. No. 44, p22.
    >>> import seawater as sw
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> sw.smow(t)
    array([ 999.842594  ,  999.842594  ,  995.65113374,  995.65113374,
            999.842594  ,  999.842594  ,  995.65113374,  995.65113374])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J. and  Poisson, A. International one-atmosphere equation
       of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
       doi:10.1016/0198-0149(81)90122-9

    """
    import numpy as np
    
    t = np.asanyarray(t)

    a = (999.842594, 6.793952e-2, -9.095290e-3, 1.001685e-4, -1.120083e-6,
         6.536332e-9)

    T68 = T68conv(t)
    return (a[0] + (a[1] + (a[2] + (a[3] + (a[4] + a[5] * T68) * T68) * T68) *
            T68) * T68)
