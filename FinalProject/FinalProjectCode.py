# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 13:45:10 2024

@author: krist
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from scipy.optimize import brentq, root, newton


class Spacetime:

    """
    Define an axisymmetric spacetime.

    For an axisymmetric, vacuum spacetime with Brill-Lindquist singularities
    the only parameters that matter is the locations of the singularities
    (i.e. their z-location) and their bare masses.

    Parameters
    ----------

    z_positions : list of float
        The location of the singularities on the z-axis.
    masses : list of float
        The bare masses of the singularities.
    reflection_symmetry : bool, optional
        Is the spacetime symmetric across the x-axis.

    See also
    --------

    TrappedSurface : class defining the trapped surfaces on a spacetime.

    Examples
    --------

    >>> schwarzschild = Spacetime([0.0], [1.0], True)

    This defines standard Schwarzschild spacetime with unit mass.

    >>> binary = Spacetime([-0.75, 0.75], [1.0, 1.1])

    This defines two black holes, with the locations mirrored but different
    masses.
    """

    def __init__(self, z_positions, masses, reflection_symmetric=False):
        """
        Initialize the spacetime given the location and masses of the
        singularities.
        """
        self.reflection_symmetric = reflection_symmetric
        self.z_positions = np.array(z_positions)
        self.masses = np.array(masses)
        self.N = len(z_positions)


class TrappedSurface:

    r"""
    Store any trapped surface, centred on a particular point.

    The trapped surface is defined in polar coordinates centred on a point
    on the z-axis; the z-axis is :math:`\theta` = 0 or :math:`\theta` =
    :math:`\pi`.

    Parameters
    ----------

    spacetime : Spacetime
        The spacetime on which the trapped surface lives.
    z_centre : float
        The z-coordinate about which the polar coordinate system describing
        the trapped surface is defined.

    See also
    --------

    Spacetime : class defining the spacetime.

    Notes
    -----

    With the restricted spacetime considered here, a trapped surface
    :math:`h(\theta)` satisfies a boundary value problem with the
    boundary conditions :math:`h'(\theta = 0) = 0 = h'(\theta = \pi)`.
    If the spacetime is reflection symmetric about the x-axis then the
    boundary condition :math:`h'(\theta = \pi / 2) = 0` can be used
    and the domain restricted to :math:`0 \le \theta \le \pi / 2`.

    The shooting method is used here. In the reflection symmetric case
    the algorithm needs a guess for the initial horizon radius,
    :math:`h(\theta = 0)`, and a single condition is enforced at
    :math:`\pi / 2` to match to the boundary condition there.

    In the general case we guess the horizon radius at two points,
    :math:`h(\theta = 0)` and :math:`h(\theta = \pi)` and continuity
    of both :math:`h` *and* :math:`h'` are enforced at the matching point
    :math:`\pi / 2`. The reason for this is a weak coordinate singularity
    on the axis at :math:`\theta = 0, \pi` which makes it difficult to
    integrate *to* these points, but possible to integrate *away* from them.

    Examples
    --------

    >>> schwarzschild = Spacetime([0.0], [1.0], True)
    >>> ts1 = TrappedSurface(schwarzschild)
    >>> ts1.find_r0([0.49, 0.51])
    >>> ts1.solve_given_r0()
    >>> print(round(ts1.r0[0], 9))
    0.5

    This example first constructs the Schwarzschild spacetime which, in this
    coordinate system, has the horizon with radius 0.5. The trapped surface
    is set up, the location of the trapped surface at :math:`\theta = 0` is
    found, which is (to the solver accuracy) at 0.5.
    """

    def __init__(self, spacetime, z_centre=0.0, reflection_symmetric=False):
        """
        Initialize a horizon centred on a particular point.
        """
        self.z_centre = z_centre
        self.spacetime = spacetime

    def expansion(self, theta, H):
        """
        Compute the expansion for the given spacetime at a fixed point.

        This function gives the differential equation defining the
        boundary value problem.

        Parameters
        ----------

        theta : float
            The angular location at this point.
        H : list of float
            A vector of :math:`(h, h')`.
        """

        h = H[0]
        dhdtheta = H[1]

        z_i = self.spacetime.z_positions
        m_i = self.spacetime.masses

        distance_i = np.zeros_like(z_i)
        z0_minus_zi = np.zeros_like(z_i)
        for i in range(len(z_i)):
            z0_minus_zi[i] = self.z_centre - z_i[i]
            distance_i[i] = np.sqrt(h ** 2 +
                                    2.0 * z0_minus_zi[i] * h * np.cos(theta)
                                    + z0_minus_zi[i] ** 2)

        C2 = 1.0 / (1.0 + (dhdtheta / h) ** 2)
        if (abs(theta) < 1e-16) or (abs(theta - np.pi) < 1e-16):
            cot_theta_dhdtheta_C2 = 0.0
        else:
            cot_theta_dhdtheta_C2 = dhdtheta / (np.tan(theta) * C2)

        psi = 1.0
        dpsi_dr = 0.0
        dpsi_dtheta = 0.0
        for i in range(len(m_i)):
            psi += 0.5 * m_i[i] / distance_i[i]
            dpsi_dr -= 0.5 * m_i[i] * (h + z0_minus_zi[i] * np.cos(theta)) /\
                distance_i[i] ** 3
            dpsi_dtheta += 0.5 * m_i[i] * h * z0_minus_zi[i] * np.sin(theta) /\
                distance_i[i] ** 3

        dHdtheta = np.zeros_like(H)
        dHdtheta[0] = dhdtheta
        dHdtheta[1] = 2.0 * h - cot_theta_dhdtheta_C2 + \
            4.0 * h ** 2 / (psi * C2) * \
            (dpsi_dr - dpsi_dtheta * dhdtheta / h ** 2) + \
            3.0 * dhdtheta ** 2 / h

        return dHdtheta

    # Define the shooting function if using matching (0 <= theta <= pi)
    def shooting_function_full(self, r0):
        r"""
        The function used in the shooting algorithm.

        This is the full algorithm from integrating over
        :math:`0 \le \theta \le \pi`. The difference between the
        solution and its derivative at the matching point is the
        error to be minimized.

        Parameters
        ----------

        r0 : list of float
            Initial guess for the horizon radius, as outlined above.

        Returns
        -------

        list of float
            The error at the matching point.
        """
        # First half of the horizon
        H0 = np.array([r0[0], 0.0])
        solver1 = ode(self.expansion)
        solver1.set_integrator("dopri5", atol=1.e-8, rtol=1.e-6)
        solver1.set_initial_value(H0, 0.0)
        solver1.integrate(np.pi / 2.0)
        # Second half of the horizon
        H0 = np.array([r0[1], 0.0])
        solver2 = ode(self.expansion)
        solver2.set_integrator("dopri5", atol=1.e-8, rtol=1.e-6)
        solver2.set_initial_value(H0, np.pi)
        solver2.integrate(np.pi / 2.0)

        return solver1.y - solver2.y

    # Define the shooting function if symmetric (0 <= theta <= pi/2)
    def shooting_function(self, r0):
        r"""
        The function used in the shooting algorithm.

        This is the symmetric algorithm from integrating over
        :math:`0 \le \theta \le \pi / 2`. The difference between the
        derivative at the end point and the boundary condition is the
        error to be minimized.

        Parameters
        ----------

        r0 : float
            Initial guess for the horizon radius, as outlined above.

        Returns
        -------

        float
            The error at the end point.
        """

        H0 = np.array([r0, 0.0])
        solver1 = ode(self.expansion)
        solver1.set_integrator("dopri5", atol=1.e-8, rtol=1.e-6)
        solver1.set_initial_value(H0, 0.0)
        solver1.integrate(np.pi / 2.0)

        return solver1.y[1]

    def find_r0(self, input_guess, full_horizon=False):
        r"""
        Given some initial guess, find the correct starting location
        for the trapped surface using shooting.

        This finds the horizon radius at :math:`\theta = 0` which,
        together with the differential equation, specifies the trapped
        surface location.

        Parameters
        ----------

        input_guess : list of float
            Two positive reals defining the guess for the initial radius.

            Note that the meaning is different depending on whether this
            is a "full" horizon or not. For a full horizon the numbers
            correspond to the guesses at :math:`\theta = 0, \pi`
            respectively. In the symmetric case where only one guess is
            needed the vector defines the interval within which a *unique*
            root must lie.

        full_horizon : bool, optional
            If the general algorithm is needed (ie, the domain should be
            :math:`0 \le \theta \le \pi` instead of
            :math:`0 \le \theta \le \pi / 2`).

            This parameter is independent of the symmetry of the spacetime.
            If the spacetime is not symmetric this parameter will be
            ignored and the general algorithm always used. If the spacetime
            is symmetric it may still be necessary to use the general
            algorithm: for example, for two singularities it is possible to
            find a trapped surface surrounding just one singularity.
        """

        # Now find the horizon given the input guess
        self.r0 = []
        if (full_horizon or
                not self.spacetime.reflection_symmetric or
                abs(self.z_centre) > 1.e-15):
            sol = root(self.shooting_function_full, input_guess, tol=1.e-12)
            self.r0 = sol.x
        else:
#            sol = brentq(self.shooting_function, input_guess[0],
#                         input_guess[1])
            sol = newton(self.shooting_function, input_guess[1])
            self.r0 = [sol]

    def Calch(self, dtheta, ai2, hi):
        k2 = dtheta * (ai2)
        
        return hi + k2
        
    def CalcParams(self, ai, hi, theta, dtheta, z_i, m_i):
        distance_i = np.zeros_like(z_i)
        z0_minus_zi = z_i
        #x = h[i-1,0] * np.cos(theta)
        #z = h[i-1,0] * np.sin(theta)
        for c in range(len(z_i)):
            distance_i[c] = np.sqrt(hi ** 2 +
                                   2.0 * z0_minus_zi[c] * hi * np.cos(theta)
                                   + z0_minus_zi[c] ** 2)
            #print("distance " + str(distance_i[c]))
           
        psi = 1.0
        dpsi_dr = 0.0
        dpsi_dtheta = 0.0
        for b in range(len(m_i)):
            psi += 0.5 * m_i[b] / distance_i[b]
            dpsi_dr -= 0.5 * m_i[b] * (hi + z0_minus_zi[b] * np.cos(theta)) /\
                distance_i[b] ** 3
            dpsi_dtheta += 0.5 * m_i[b] * hi * z0_minus_zi[b] * np.sin(theta) /\
               distance_i[b] ** 3
           
        C2 = 1.0 / (1.0 + (ai / hi) ** 2)
           
        if (abs(theta) <= np.pi/150) or (abs(theta - np.pi) <= np.pi/150) \
            or (abs(2*np.pi - theta) <= np.pi/150):
            cot_theta_dhdtheta_C2 = 0.0
        else:
            cot_theta_dhdtheta_C2 = ai / (np.tan(theta) * C2)
               
        return psi, dpsi_dr, dpsi_dtheta, C2, cot_theta_dhdtheta_C2

    def CalcA(self, ai, hi, theta, dtheta, z_i, m_i):
        psi, dpsi_dr, dpsi_dtheta, C2, cot_theta_dhdtheta_C2 = \
            self.CalcParams(ai, hi, theta, dtheta, z_i, m_i)
        
        k1 = dtheta * (2*hi - cot_theta_dhdtheta_C2 + \
                    (4.0 * hi ** 2 / (psi * C2)) * \
                    (dpsi_dr - dpsi_dtheta * ai / hi ** 2) + \
                    3.0 * ai ** 2 / hi)
        
        psi, dpsi_dr, dpsi_dtheta, C2, cot_theta_dhdtheta_C2 = \
            self.CalcParams((ai + k1/2), hi, (theta + (dtheta/2)), dtheta, z_i, m_i)
            
        k2 = dtheta * (2*hi - cot_theta_dhdtheta_C2 + \
                    (4.0 * hi ** 2 / (psi * C2)) * \
                    (dpsi_dr - dpsi_dtheta * (ai + (k1/2)) / hi ** 2) + \
                    3.0 * (ai + (k1/2)) ** 2 / hi)
        return ai + k2

    def CalcHorizon(self, dtheta, FullCircle):
        """
        Parameters
        ----------
        dtheta : float
            This is the size of the change in theta for each step
        FullCircle : float
            This is what theta will go to to complete a full circle = 2pi

        Returns
        -------
        None.

        """
        # define step size
        Steps = round(FullCircle/dtheta)
        
        a = np.zeros((Steps+1,1),'float')
        h = np.zeros((Steps+1,2),'float')
        
        # initial conditions
        a[0,0] = 0
        h[0,0] = self.r0[0]
        theta = 0 + dtheta
        
        #print(self.r0[0])
        
        z_i = self.spacetime.z_positions
        m_i = self.spacetime.masses
        
        for i in range(1,Steps+1):
            a2 = self.CalcA(a[i-1,0], h[i-1,0], theta, dtheta, z_i, m_i)
            
            h2 = self.Calch(dtheta, a2, h[i-1,0])
                    
            if (i == 249 or i == 498):
                a[i,0] = 0
                h[i,0] = self.r0[1]
                #print(h[i,0])
                #h[i,0] = h2
            else:
                a[i,0] = a2
                    
                #print(str(a[i,0]) + " a")
                h[i,0] = h2
                #print(str(h[i,0]) + " h " + str(i))
            h[i,1] = theta
                
            theta += dtheta
        
        self.H = h
        return None

    def convert_to_cartesian(self):
        """
        When the solution is known in r, theta coordinates, compute
        the locations in cartesian coordinates (2 and 3d).

        This function assumes that the trapped surface has been located and
        solved for.

        See also
        --------

        solve_given_r0 : find the trapped surface location in polar
                         coordinates.
        """

        self.x = self.H[:, 0] * np.sin(self.H[:,1])
        self.z = self.z_centre + self.H[:, 0] * np.cos(self.H[:,1])

        return None
        
def find_horizon_binary(z, mass1, mass2):
    r"""
    Utility function to find horizons for the general case.

    This returns the horizon for a spacetime with precisely two singularities
    of mass [mass1, mass2] located at :math:`\pm z`. That is, we work in the
    frame where the location of the horizons is symmetric.

    Notes
    -----

    The initial guess for the horizon location is based on fitting a cubic
    to the results constructed for :math:`0 \le z \le 0.75` for the unit
    mass case. The radius should scale with the mass. For larger separations
    we should not expect a common horizon.

    Parameters
    ----------

    z : float, optional
        The distance from the origin of the singularities (ie the two
        singularities are located at [-z, +z]).
    mass : float, optional
        The mass of the singularities.

    Returns
    -------

    ts : TrappedSurface
        Only returns the single surface found, expected to be the common
        horizon.
    """

    st = Spacetime([-z, z], [mass1, mass2])
    ts = TrappedSurface(st, 0.0)
    # An empirical formula for the required initial guess
    # (ie the value of r0, or h, at theta = 0)
    # This really is just a guess based on the symmetric case.
    zom = 2.0 * z / (mass1 + mass2)
    r0_empirical = (1.0 - 0.0383 * zom + 0.945 * zom ** 2 -
                    0.522 * zom ** 3) * \
        (mass1 + mass2) / 2.0
    r0_empirical = max(r0_empirical, z + 0.5 * max(mass1, mass2))
    initial_guess = [r0_empirical, r0_empirical]
    ts.find_r0(initial_guess, True)
    ts.CalcHorizon(np.pi/249, 2*np.pi)
    ts.convert_to_cartesian()
    return ts
"""
zs = [0.76, 0.66, 0.56, 0.46, 0.36, 0.26, 0.16, 0.06]
print(zs)

figs = []
for a in range(len(zs)):
    figs.append(find_horizon_binary(zs[a], 1, 1))
"""
ts = find_horizon_binary(0.78, 1, 1)

fig, host = plt.subplots(figsize=(6,8), layout='constrained')

host.set_title("Simulated Apparent Horizon")
host.set_xlabel("x-axis")
host.set_ylabel("z-axis")

host.plot(ts.x, ts.z, 'b-')
for z, m in zip(ts.spacetime.z_positions, ts.spacetime.masses):
    host.plot(0.0, z,
            'kx', markersize=10, markeredgewidth=1 + int(round(m)))
"""    
count = 0
for ts in figs:
    word = "d = " + str(zs[count]*2)    
    host.plot(ts.x, ts.z, label = word)
    count += 1
    for z, m in zip(ts.spacetime.z_positions, ts.spacetime.masses):
        host.plot(0.0, z,
                'kx', markersize=1, markeredgewidth=1 + int(round(m)))
     """  
plt.legend()
