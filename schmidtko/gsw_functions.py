'''
This file contains functions from python-gsw package
SP_from_SA was adjusted for values around the coast of Antarctica
(reese@pik-potsdam.de)

Copyright Notice and Statement for the gsw project:

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''

import numpy as np


################## constants ##############################

SSO = 35.16504
"""
SSO is the Standard Ocean Reference Salinity (35.16504 g/kg.)
SSO is the best estimate of the Absolute Salinity of Standard Seawater
when the seawater sample has a Practical Salinity, SP, of 35
(Millero et al., 2008), and this number is a fundamental part of the
TEOS-10 definition of seawater.
References:
-----------
.. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
   seawater - 2010: Calculation and use of thermodynamic properties.
   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
   UNESCO (English), 196 pp. See appendices A.3, A.5 and Table D.4.
.. [2] Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008:
   The composition of Standard Seawater and the definition of the
   Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.
   See Table 4 and section 5.
"""

sfac = 0.0248826675584615
"""
sfac = 1 / (40. * (SSO / 35.))
"""

cp0 = 3991.86795711963
"""
The "specific heat" for use with Conservative Temperature. cp0 is the ratio
of potential enthalpy to Conservative Temperature.
See Eqn. (3.3.3) and Table D.5 from IOC et al. (2010).
"""

T0 = Kelvin = 273.15
"""
The Celsius zero point; 273.15 K.  That is T = t + T0 where T is the
Absolute Temperature (in degrees K) and t is temperature in degrees C.
"""

######################## functions #######################


def pt_from_CT(SA, CT):
    '''
    Copied from gsw.pt_from_CT()
    '''
    SA, CT, mask = strip_mask(SA, CT)
    SA = np.maximum(SA, 0)

    s1 = SA * 35. / SSO

    a0 = -1.446013646344788e-2
    a1 = -3.305308995852924e-3
    a2 = 1.062415929128982e-4
    a3 = 9.477566673794488e-1
    a4 = 2.166591947736613e-3
    a5 = 3.828842955039902e-3

    b0 = 1.000000000000000e+0
    b1 = 6.506097115635800e-4
    b2 = 3.830289486850898e-3
    b3 = 1.247811760368034e-6

    a5CT = a5 * CT
    b3CT = b3 * CT
    CT_factor = (a3 + a4 * s1 + a5CT)
    pt_num = a0 + s1 * (a1 + a2 * s1) + CT * CT_factor
    pt_den = b0 + b1 * s1 + CT * (b2 + b3CT)
    pt = pt_num / pt_den

    dCT_dpt = pt_den / (CT_factor + a5CT - (b2 + b3CT + b3CT) * pt)

    # 1.5 iterations through the modified Newton-Rapshon iterative method
    CT_diff = CT_from_pt(SA, pt) - CT
    pt_old = pt
    pt = pt_old - CT_diff / dCT_dpt  # 1/2-way through the 1st modified N-R.
    ptm = 0.5 * (pt + pt_old)

    # This routine calls gibbs_pt0_pt0(SA, pt0) to get the second derivative of
    # the Gibbs function with respect to temperature at zero sea pressure.

    dCT_dpt = -(ptm + Kelvin) * gibbs_pt0_pt0(SA, ptm) / cp0
    pt = pt_old - CT_diff / dCT_dpt  # End of 1st full modified N-R iteration.
    CT_diff = CT_from_pt(SA, pt) - CT
    pt_old = pt
    pt = pt_old - CT_diff / dCT_dpt  # 1.5 iterations of the modified N-R.
    # Abs max error of result is 1.42e-14 deg C.
    return np.ma.array(pt, mask=mask, copy=False)



def CT_from_pt(SA, pt):
    """
    Calculates Conservative Temperature of seawater from potential
    temperature (whose reference sea pressure is zero dbar).
    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]
    Returns
    -------
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (ITS-90)]
    Examples
    --------
    >>> import gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> pt = [28.7832, 28.4209, 22.7850, 10.2305, 6.8292, 4.3245]
    >>> gsw.CT_from_pt(SA, pt)
    array([ 28.80992302,  28.43914426,  22.78624661,  10.22616561,
             6.82718342,   4.32356518])
    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
       of seawater - 2010: Calculation and use of thermodynamic properties.
       Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
       UNESCO (English), 196 pp. See section 3.3.
    """

    SA, pt, mask = strip_mask(SA, pt)

    pot_enthalpy = pot_enthalpy_from_pt(SA, pt)

    CT = pot_enthalpy / cp0

    return np.ma.array(CT, mask=mask, copy=False)

def pot_enthalpy_from_pt(SA, pt):
    """
    Calculates the potential enthalpy of seawater from potential
    temperature (whose reference sea pressure is zero dbar).
    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]
    Returns
    -------
    pot_enthalpy : array_like
                   potential enthalpy [J kg :sup:`-1`]
    Notes
    -----
    TODO
    Examples
    --------
    >>> import gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> pt = [28.7832, 28.4209, 22.7850, 10.2305, 6.8292, 4.3245]
    >>> gsw.pot_enthalpy_from_pt(SA, pt)
    array([ 115005.40853458,  113525.30870246,   90959.68769935,
             40821.50280454,   27253.21472227,   17259.10131183])
    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
       of seawater - 2010: Calculation and use of thermodynamic properties.
       Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
       UNESCO (English), 196 pp. See section 3.2.
    """

    SA, pt, mask = strip_mask(SA, pt)

    SA = np.maximum(SA, 0)

    x2 = sfac * SA
    x = np.sqrt(x2)
    y = pt * 0.025  # Normalize for F03 and F08

    pot_enthalpy = (61.01362420681071 + y * (168776.46138048015 +
    y * (-2735.2785605119625 + y * (2574.2164453821433 +
    y * (-1536.6644434977543 + y * (545.7340497931629 +
    (-50.91091728474331 - 18.30489878927802 * y) * y))))) +
    x2 * (268.5520265845071 + y * (-12019.028203559312 +
    y * (3734.858026725145 + y * (-2046.7671145057618 +
    y * (465.28655623826234 + (-0.6370820302376359 -
    10.650848542359153 * y) * y)))) +
    x * (937.2099110620707 + y * (588.1802812170108 +
    y * (248.39476522971285 + (-3.871557904936333 -
    2.6268019854268356 * y) * y)) +
    x * (-1687.914374187449 + x * (246.9598888781377 +
    x * (123.59576582457964 - 48.5891069025409 * x)) +
    y * (936.3206544460336 +
    y * (-942.7827304544439 + y * (369.4389437509002 +
    (-33.83664947895248 - 9.987880382780322 * y) * y)))))))

    # The above polynomial for pot_enthalpy is the full expression for
    # potential enthalpy in terms of SA and pt, obtained from the Gibbs function
    # as below.  It has simply collected like powers of x and y so that it is
    # computationally faster than calling the Gibbs function twice as is done in
    # the commented code below. When this code below is run, the results are
    # identical to calculating pot_enthalpy as above, to machine precision.

    # g000 = gibbs(n0, n0, n0, SA, pt, 0)
    # g010 = gibbs(n0, n1, n0, SA, pt, 0)
    # pot_enthalpy = g000 - (Kelvin + pt) * g010

    # This is the end of the alternative code
    # %timeit gsw.CT_from_pt(SA, pt)
    # 1000 loops, best of 3: 1.34 ms per loop <- calling gibbs
    # 1000 loops, best of 3: 254 us per loop <- standard

    return np.ma.array(pot_enthalpy, mask=mask, copy=False)


def gibbs_pt0_pt0(SA, pt0):
    """
    Calculates the second derivative of the specific Gibbs function with
    respect to temperature at zero sea pressure or _gibbs(0,2,0,SA,t,0).
    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt0 : array_like
          potential temperature relative to 0 dbar [:math:`^\circ` C (ITS-90)]
    Returns
    -------
    gibbs_pt0_pt0 : array_like
                    TODO: write the eq. for the second derivative of the
                    specific Gibbs function. FIXME: [units]
    Notes
    -----
    This library function is called by both "pt_from_CT(SA,CT)"
    and "pt0_from_t(SA,t,p)".
    """

    SA, pt0, mask = strip_mask(SA, pt0)
    x2 = sfac * SA
    x = np.sqrt(x2)
    y = pt0 * 0.025
    g03 = (-24715.571866078 +
           y * (4420.4472249096725 +
           y * (-1778.231237203896 +
           y * (1160.5182516851419 +
           y * (-569.531539542516 + y * 128.13429152494615)))))

    g08 = (x2 * (1760.062705994408 + x * (-86.1329351956084 +
           x * (-137.1145018408982 + y * (296.20061691375236 +
           y * (-205.67709290374563 + 49.9394019139016 * y))) +
           y * (-60.136422517125 + y * 10.50720794170734)) +
           y * (-1351.605895580406 + y * (1097.1125373015109 +
           y * (-433.20648175062206 + 63.905091254154904 * y)))))
    gibbs_pt0_pt0 = (g03 + g08) * 0.000625
    return np.ma.array(gibbs_pt0_pt0, mask=mask, copy=False)




def SP_from_SA_Antarctica(SA):
    '''
    based on gsw.SP_from_SA
    This function uses SAAR (the Salinity Anomaly Ratio) which is obtained from a data
    set/look-up table and not availabe close to the Antarctic continent. 
    I therefore compute the average saar value where availabe and used this as a proxy 
    '''
    
    saar_for_Antarctica = 0.00023004095217 # from calculate_reference_saar.py
    SP = (35.0 / 35.16504) * SA / (1.0 + saar_for_Antarctica)

    return SP


#### original frunction from gsw package
# def SP_from_SA(SA, p, lon, lat):
#     """
#     TODO: docstring
#     """
#     # maybe add input checking...

#     saar, in_ocean = SAAR(p, lon, lat)
#     SP = (35 / 35.16504) * SA / (1.0 + saar)

#     SP_baltic = SP_from_SA_Baltic(SA, lon, lat)
#     bmask = SP_baltic.mask
#     if bmask is not np.ma.nomask and not bmask.all():
#         inbaltic = ~bmask
#         SP[inbaltic] = SP_baltic[inbaltic]

# return SP, in_ocean


def p_from_z(z, lat, geo_strf_dyn_height=0):
    """
    Copied from gsw.p_from_z()

    Calculates sea pressure from height using computationally-efficient
    48-term expression for density, in terms of SA, CT and p (McDougall et al.,
    2011).  Dynamic height anomaly, geo_strf_dyn_height, if provided, must be
    computed with its pr=0 (the surface.)
    Parameters
    ----------
    z : array_like
        height [m]
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    geo_strf_dyn_height : float, optional
                          dynamic height anomaly [ m :sup:`2` s :sup:`-2` ]
                          The reference pressure (p_ref) of geo_strf_dyn_height
                          must be zero (0) dbar.
    Returns
    -------
    p : array_like
        pressure [dbar]
    Examples
    --------
    >>> import gsw
    >>> z = [-10., -50., -125., -250., -600., -1000.]
    >>> lat = 4.
    >>> gsw.p_from_z(z, lat)
    array([   10.05572704,    50.28354425,   125.73185732,   251.54028663,
             604.2099135 ,  1007.9900587 ])
    >>> -gsw.z_from_p(gsw.p_from_z(z, lat), lat)
    array([  10.,   50.,  125.,  250.,  600., 1000.])
    Notes
    -----
    Height (z) is NEGATIVE in the ocean. Depth is -z. Depth is not used in the
    gibbs library.
    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
       of seawater - 2010: Calculation and use of thermodynamic properties.
       Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
       UNESCO (English), 196 pp.
    .. [2] McDougall T.J., P.M. Barker, R. Feistel and D.R. Jackett, 2011: A
       computationally efficient 48-term expression for the density of seawater
       in terms of Conservative Temperature, and related properties of seawater.
    .. [3] Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74,
       128-133.
    .. [4] Saunders, P. M., 1981: Practical conversion of pressure to depth.
       Journal of Physical Oceanography, 11, 573-574.
    """

    X = np.sin(lat * DEG2RAD)
    sin2 = X ** 2
    gs = 9.780327 * (1.0 + (5.2792e-3 + (2.32e-5 * sin2)) * sin2)

    # get the first estimate of p from Saunders (1981)
    c1 = 5.25e-3 * sin2 + 5.92e-3
    p = -2 * z / ((1 - c1) + np.sqrt((1 - c1) * (1 - c1) + 8.84e-6 * z))

    df_dp = db2Pascal * specvol_SSO_0_p(p)  # Initial value for f derivative.

    f = (enthalpy_SSO_0_p(p) + gs *
         (z - 0.5 * gamma * (z ** 2)) - geo_strf_dyn_height)

    p_old = p
    p = p_old - f / df_dp
    p_mid = 0.5 * (p + p_old)
    df_dp = db2Pascal * specvol_SSO_0_p(p_mid)
    p = p_old - f / df_dp

    # After this one iteration through this modified Newton-Raphson iterative
    # procedure, the remaining error in p is at computer machine precision,
    # being no more than 1.6e-10 dbar.

    return p


def strip_mask(*args):
    """
    copied from python-gsw utilities

    Process the standard arguments for efficient calculation.
    Return unmasked arguments, plus a mask.
    The first argument, SA, is handled specially so that it can be
    This could be absorbed into a decorator, but it would
    require redefining functions to take the additional
    mask argument or kwarg.
    """

    mask = np.ma.getmaskarray(args[-1])
    SA = args[0]
    if SA.shape:
        SA = np.ma.asarray(SA)
        SA[SA < 0] = np.ma.masked
        for a in args[:-1]:
            mask = np.ma.mask_or(mask, np.ma.getmask(a))
        newargs = [SA.filled(0)]
    elif SA < 0:
        SA = 0
        for a in args[1:-1]:
            mask = np.ma.mask_or(mask, np.ma.getmask(a))
        newargs = [SA]
    newargs.extend([np.ma.filled(a, 0) for a in args[1:]])
    newargs.append(mask)
    return newargs