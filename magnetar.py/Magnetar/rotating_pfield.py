
import numpy as np

def oldpsicalc(mass, radius, b):
    from scipy.integrate import quad
    MoverR = mass / radius
    a2 = (1 - 2 * MoverR) * MoverR * MoverR
    x = b / radius * (1 - 2 * MoverR)**0.5
    # print(x)
    x2 = x * x
    f = (lambda u: 1 / (a2 - (1 - 2 * u) * u * u * x2)**0.5)
    return x * quad(f, 0, MoverR)[0]


# for a given impact parameter b over what angle psi does the light travel [radians]
def psicalc(mass, radius, b):
    from scipy.integrate import quad
    MoverR = mass / radius
    x = b / mass
    # print(x)
    x2 = x * x
    f = (lambda u: (1 - (1 - 2 * u) * u * u * x2)**-0.5)
    return x * quad(f, 0, MoverR)[0]


# Eq 1 of https://arxiv.org/pdf/astro-ph/0201117.pdf
def psiapprox(mass, radius, b):
    MoverR = mass / radius
    rinf = radius / (1 - 2 * MoverR)**0.5
    sinalpha = b / rinf
    cosalpha = (1 - sinalpha**2)**0.5
    cospsi = 1 - (1 - cosalpha) / (1 - 2 * MoverR)
    return np.arccos(cospsi)


# how much is a photon delayed from the surface for a given impact parameter b compared to zero impact parameter (sub-observer point) [same units as mass]
def timecalc(mass, radius, b):
    from scipy.integrate import quad
    MoverR = mass / radius
    x = b / mass
    # print(x)
    x2 = x * x
    f = (
        lambda u: ((1 - (1 - 2 * u) * u * u * x2)**-0.5 - 1) / (1 - 2 * u) / u / u
    )
    return mass * quad(f, 0, MoverR)[0]


# Approximate results inspired by https://arxiv.org/pdf/astro-ph/0201117.pdf
# good to two percent for R>5M (much better if b<0.9 rinf)
def timeapprox(mass, radius, b):
    uu = 1 - np.cos(psiapprox(mass, radius, b))
    return radius * radius * uu / (radius - mass * uu / 3.0)


# calculate both together
def calc_psiandtime(mass, radius, b):
    res1 = psiapprox(mass, radius, b)
    uu = 1 - np.cos(res1)
    return np.degrees(res1), radius * radius * uu / (radius - mass * uu / 3.0)

from Magnetar.pfield import pfield

class rotating_pfield(pfield):
    def setrotationvector(self, rate, inclination, field_angle):
        field_angle = np.abs(field_angle)
        inclination = np.abs(inclination)
        assert(self.theta<=np.abs(field_angle+inclination) and self.theta>=np.abs(field_angle-inclination)),\
            'Magnetic field direction must lie between |inclination=-field_angle| and inclination+field_angle'

        self.rate = rate
        self.inclination = np.radians(inclination)
        self.field_angle = np.radians(field_angle)
        self.phase = 0
        self.phase_delay = []
        self.psi = []
        self.time = []
        self.energy = []
        self.spin_alt = []
        self.spin_colat = []
        self.alpha = []
        self.zeta = []
        self.spin_azim = []
        self.pt = 0
        self.pr = []
        self.pth = []
        self.pph = []
        self.ut = []
        self.ur = []
        self.uth = []
        self.uph = []
        self.delay_mag = []
        self.theta_delay = []
        self.eta_delay = []

        #=== SPIN FRAME ===#
        # calculate the value of psi and light-travel time  at each point
        self.psi, self.time = calc_psiandtime(self.mass, self.radius,
                                              self.data['b'])
        # calculate the phase using spherical trig at the time that light leaves the sub-observer point
        self.phase = np.arccos(
            (np.cos(np.radians(self.theta)) -
             np.cos(self.field_angle) * np.cos(self.inclination)) /
            (np.sin(self.field_angle) * np.sin(self.inclination)))
        # calculate the colatitude of each surface elemente relative to the spin axis
        self.eta = np.arccos(
            (np.cos(self.field_angle) -
             np.cos(self.inclination) * np.cos(np.radians(self.theta))) /
            (np.sin(self.inclination) * np.sin(np.radians(self.theta))))
        self.spin_colat = np.arccos(
            np.cos(self.inclination) * np.cos(np.radians(self.psi)) +
            np.sin(self.inclination) * np.sin(np.radians(self.psi)) * np.
            cos(np.pi - self.data['beta'] - self.eta))
        # zeta angle between n-b plane and n-los plane
        beta = np.where(self.data['beta'] > np.pi, self.data['beta'] - np.pi,
                        self.data['beta'])
        sin_zeta = np.sin(np.radians(self.theta)) * np.sin(
            np.pi - beta) / np.sin(np.radians(self.data['mag_colat']))
        zeta = np.where(sin_zeta > 1.0, np.pi / 2.0, np.arcsin(sin_zeta))
        self.zeta = np.where(
            np.cos(np.radians(self.theta)) <
            np.cos(np.radians(self.data['mag_colat'])) * np.cos(
                np.radians(self.psi)), np.pi - zeta, zeta)
        # alpha angle between n-s plane and n-los plane
        xi = np.where((self.data['beta'] < (np.pi - self.eta)),
                      np.pi - self.data['beta'] - self.eta,
                      self.eta + self.data['beta'] - np.pi)
        xi = np.where((self.data['beta'] > (2.0 * np.pi - self.eta)),
                      2.0 * np.pi + (np.pi - self.data['beta'] - self.eta), xi)
        sin_alpha = np.sin(xi) * np.sin(self.inclination) / np.sin(
            self.spin_colat)
        alpha = np.where(sin_alpha > 1.0, np.pi / 2.0, np.arcsin(sin_alpha))
        self.alpha = np.where(
            np.cos(self.inclination) <
            np.cos(self.spin_colat) * np.cos(np.radians(self.psi)),
            np.pi - alpha, alpha)
        # calculate the azimuth in the rotating frame
        #        aux_azim1 = np.pi - (
        #            self.alpha + np.radians(self.data['phi']) + self.zeta - np.pi)
        #        aux_azim2 = np.pi - (
        #            self.alpha + np.radians(self.data['phi']) - self.zeta)

        aux_azim1 = (
            self.alpha + np.radians(self.data['phi']) + self.zeta - np.pi)
        aux_azim2 = (self.alpha + np.radians(self.data['phi']) - self.zeta)
        self.spin_azim = np.where(self.data['mag_colat'] < 90.0, aux_azim1,
                                  aux_azim2)
        self.spin_alt = np.pi / 2.0 - np.radians(self.data['theta'])

        #=== DOPPLER SHIFT ===#
        # Diagonal components (+,-,-,-) of the metric tensor g_ii
        gtt = 1 - 2 * self.mass / self.radius
        grr = -1.0 / gtt
        gthth = -self.radius**2.0
        gphph = -(self.radius * np.sin(self.spin_colat))**2.0
        # calculate the four velocity u^i= (ut, ur, uth, uph).
        # self.rate  with 1/c factor to convert [hz] to [1/cm]
        c = 3.0e+10
        self.ut = 1.0 / np.sqrt(gtt - gphph * (self.rate / c)**2)
        self.ur = 0
        self.uth = 0
        self.uph = self.ut * (self.rate / c)
        # calculate the photon four momentum p_i=(pt,pr,pth,pph) in the spin frame using p_t=1
        self.pt = 1.0
        #self.pr   = not needed
        #self.pth  = not needed
        self.pph = np.sqrt((-self.pt**2 * gphph / gtt) /
                           ((np.tan(self.spin_alt))**2 *
                            (1 + 1 / (np.tan(self.spin_azim))**2) + 1 /
                            (np.tan(self.spin_azim))**2 + 1))
        # defines where pph is positive or negative on the sphere
        self.pph = np.where((self.data['beta'] < np.pi - self.eta) |
                            (self.data['beta'] > 2.0 * np.pi - self.eta),
                            -self.pph, self.pph)
        #self.pph  = np.where(self.spin_azim>np.pi/2.0,-self.pph,self.pph)
        # calculate the photon energy E = p_i * u^i at the star surface
        self.energy = self.pt * self.ut + self.pph * self.uph

        #=== TIME DELAY ===#
        # because of the time delay  and the star rotation, for every b>0 there is a slight change in the direction of the
        # mag axis,  which should be recomputed from phase_delay. One option is to assign self.phase to b=0, then for b>=0
        self.phase_delay = self.phase - self.rate * self.time / c
        # compute the theta_delay angle between the magnetic axis and LOS
        cos_theta_delay = np.cos(self.field_angle) * np.cos(
            self.inclination) + np.sin(self.field_angle) * np.sin(
                self.inclination) * np.cos(self.phase_delay)
        self.theta_delay = np.arccos(cos_theta_delay)
        # same as the lines before but for the eta angle
        cos_eta_delay = (
            np.cos(self.field_angle) -
            np.cos(self.inclination) * np.cos(self.theta_delay)) / (
                np.sin(self.inclination) * np.sin(self.theta_delay))
        self.eta_delay = np.arccos(cos_eta_delay)
        # cartesians components of the spin axis
        sx = np.sin(self.inclination) * np.cos(self.eta)
        sy = np.sin(self.inclination) * np.sin(self.eta)
        sz = np.cos(self.inclination)
        # the components of the spin axis and delayed magnetic axis satisfy the system
        #     (1) bx*sx + by*sy + bz*sz = cos(field_angle)
        #     (2) by**2 + by**2 = (sin(theta_delay))**2
        #     (3) bz = cos(theta_delay)
        # this yields an equation a1*bx**2 + a2*bx + c = 0, whose coeficients are
        a1 = sx**2.0 + sy**2.0
        a2 = -2.0 * sx * (
            np.cos(self.field_angle) - np.cos(self.theta_delay) * sz)
        a3 = (np.cos(self.field_angle) - np.cos(self.theta_delay) * sz
              )**2.0 - (np.sin(self.theta_delay) * sy)**2.0
        # calculate the solution of the polynomial equation and the components (bx,by,bz) of the delayed magnetic axis
        delta1 = a2**2.0 - 4.0 * a1 * a3
        self.bx = np.where(delta1 < 0, -a2 / (2.0 * a1),
                           (-a2 + np.sqrt(delta1)) / (2.0 * a1))
        delta2 = (np.sin(self.theta_delay))**2.0 - self.bx**2.0
        self.by = np.where(delta2 < 0, 0, -np.sqrt(delta2))
        self.bz = np.cos(self.theta_delay)
        # cartesians components of the surface element
        nx = np.sin(np.radians(self.psi)) * np.cos(np.pi - self.data['beta'])
        ny = np.sin(np.radians(self.psi)) * np.sin(np.pi - self.data['beta'])
        nz = np.cos(np.radians(self.psi))
        # calculate the magnetic colatitude using cos(mag_col) = b*n
        self.mag_delay = np.arccos(self.bx * nx + self.by * ny + self.bz * nz)
##############################################################################################################################
        #calculate  zeta delay. I need theta delay and mag_delay
        # zeta angle between n-b plane and n-los plane
#        beta = np.where(self.data['beta'] > np.pi, self.data['beta'] - np.pi,
#                        self.data['beta'])
#        beta = np.where(self.data['beta'] > np.pi, self.data['beta'],
#                        self.data['beta'])
#        sin_zeta_delay = np.sin(self.theta_delay) * np.sin(
#            np.pi - beta) / np.sin(self.mag_delay)
#        zeta = np.where(sin_zeta_delay > 1.0, np.pi / 2.0, np.arcsin(sin_zeta_delay))
#        self.zeta_delay = np.where(np.cos(self.theta_delay)<np.cos(self.mag_delay) * np.cos(np.radians(self.psi)),
#                                   np.pi - zeta, zeta)
        cos_zeta_delay  = (np.cos(self.theta_delay)-
                           np.cos(self.mag_delay)*np.cos(np.radians(self.psi)))/(np.sin(self.mag_delay)*np.sin(np.radians(self.psi)))
        self.zeta_delay = np.arccos(cos_zeta_delay)

##############################################################################################################################
        #=== LIGHT ABERRATION ===#
        zen = np.radians(self.data['theta'])
        # photon azimuth angle in the spin frame. Notice that self.spin_azim is between [0,pi].
        # The azimuth angle is redefined below between [-pi,pi], consistent with the
        # range of the output of the np.arctan2 function used later.
        azim = np.where((self.data['beta'] < np.pi - self.eta) |
                        (self.data['beta'] > 2.0 * np.pi - self.eta),
                        self.spin_azim, -self.spin_azim)
        self.spin_azim2 = azim
        # calculate v/c of the surface element. A correction due to gravitational redshit migth be needed
        self.beta_spin = self.rate * self.radius * np.sin(self.spin_colat) / c
        # aberration in the spin frame
        self.spin_zen_aberr = np.arccos(
            np.cos(zen) * np.sqrt(1 - self.beta_spin**2.0) /
            (1 - np.sin(zen) * np.sin(azim) * self.beta_spin))
        # arctan2 output between [-pi,pi]. the abs function narrow it to [0,pi]
        #        self.spin_azim_aberr = np.abs(np.arctan2(
        #            np.sin(zen) * np.sin(azim) - self.beta_spin,
        #            np.sin(zen) * np.cos(azim) * np.sqrt(1 - self.beta_spin**2.0)))

        self.spin_azim_aberr = (np.arctan2(
            np.sin(zen) * np.sin(azim) - self.beta_spin,
            np.sin(zen) * np.cos(azim) * np.sqrt(1 - self.beta_spin**2.0)))

        # calculate aberration back to the magnetic frame for different regions over the sphere
        mag_azim_aberr1 =  self.alpha + self.zeta_delay + self.spin_azim_aberr
        mag_azim_aberr2 = -self.alpha + self.zeta_delay + self.spin_azim_aberr
        mag_azim_aberr3 = -self.alpha - self.zeta_delay + self.spin_azim_aberr
        mag_azim_aberr4 =  self.alpha - self.zeta_delay + self.spin_azim_aberr
        mag_azim_aberr5 = -self.alpha + self.zeta_delay + self.spin_azim_aberr + 2.0 * np.pi
        mag_azim_aberr6 = -self.alpha - self.zeta_delay + self.spin_azim_aberr
        mag_azim_aberr7 = -self.alpha + self.zeta_delay + self.spin_azim_aberr
        mag_azim_aberr8 = -self.alpha - self.zeta_delay + self.spin_azim_aberr + 2.0 * np.pi

        #conditons to identify different regions over the sphere
        condition1 =  (self.data['beta'] > np.pi -  self.eta) & (self.data['beta'] < 2.0 * np.pi - self.eta)
        condition2 =   self.data['beta'] < np.pi + np.abs(self.eta_delay - self.eta)
        condition3 =  (self.data['beta'] < np.pi -  self.eta) & (self.data['beta'] > np.abs(self.eta_delay - self.eta))
        condition4 =  (self.spin_azim_aberr>0)
        condition5 =  np.radians(self.psi)>np.pi/2.0

        mag_azim_aberr = np.where(condition1,
                                  np.where(condition2,
                                           mag_azim_aberr1,
                                           mag_azim_aberr4),
                                  np.where(condition3,
                                           np.where(condition4,
                                                    mag_azim_aberr2,
                                                    np.where(condition5,mag_azim_aberr7,mag_azim_aberr5)),
                                           np.where(condition4,
                                                    mag_azim_aberr3,
                                                    np.where(condition5,mag_azim_aberr8,mag_azim_aberr6))))

        mag_azim_aberr = np.where(mag_azim_aberr > np.pi,
                                  mag_azim_aberr - 2.0 * np.pi,
                                  mag_azim_aberr)

        mag_azim_aberr = np.where(mag_azim_aberr < -np.pi,
                                  mag_azim_aberr + 2.0 * np.pi,
                                  mag_azim_aberr)

        self.mag_azim_aberr = mag_azim_aberr

        return self








'''
        #regions over the sphere
        #region1 = (self.data['beta'] > np.pi - self.eta) & (self.data['beta'] < 2.0 * np.pi - self.eta)
        #region2 =  self.data['beta'] < np.pi
        #region3 =  self.data['beta'] < np.pi - self.eta
        #region4 =  self.spin_azim_aberr>0
        #region5 =  np.radians(self.psi)>np.pi/2.0

        mag_azim_aberr = np.where(region1,
                                  np.where(region2,mag_azim_aberr1,mag_azim_aberr4),
                                  np.where(region3,
                                           np.where(region4,
                                                    mag_azim_aberr2,
                                                    np.where(region5,mag_azim_aberr7,mag_azim_aberr5)),
                                           np.where(region4,
                                                    mag_azim_aberr3,
                                                    np.where(region5,mag_azim_aberr8,mag_azim_aberr6))))

        mag_azim_aberr = np.where(np.abs(mag_azim_aberr) > np.pi,
                                  np.where(self.data['beta'] < np.pi, mag_azim_aberr - 2.0 * np.pi, 2.0 * np.pi + mag_azim_aberr),
                                  mag_azim_aberr)

        mag_azim_aberr = np.where(mag_azim_aberr>np.pi,
                                  mag_azim_aberr - 4.0*np.pi,
                                  mag_azim_aberr)

#       condition4 = (self.data['beta'] < np.pi - self.eta) | (self.data['beta'] > 2.0 * np.pi - self.eta)
#       condition5 =  self.data['beta'] < np.pi - self.eta
#       condition6 =  self.data['beta'] > 2.0 * np.pi - self.eta
#        mag_azim_aberr = np.where(condition1,
#                                  np.where(condition2,mag_azim_aberr1,mag_azim_aberr4),
#                                  np.where(condition3,mag_azim_aberr2,mag_azim_aberr3))

#        mag_azim_aberr = np.where(condition4,
#                                  np.where(self.spin_azim_aberr < 0,
#                                           np.where(self.data['beta'] < np.pi,mag_azim_aberr5,mag_azim_aberr6),
#                                           mag_azim_aberr),
#                                  mag_azim_aberr)
#        mag_azim_aberr = np.where(condition5,
#                                  np.where(self.spin_azim_aberr < 0,
#                                           np.where(np.radians(self.psi) > np.pi / 2.0,mag_azim_aberr7,mag_azim_aberr),
#                                           mag_azim_aberr),
#                                  mag_azim_aberr)
#        mag_azim_aberr = np.where(condition6,
#                                  np.where(self.spin_azim_aberr < 0,
#                                           np.where(np.radians(self.psi) > np.pi / 2.0,mag_azim_aberr8,mag_azim_aberr), 
#                                           mag_azim_aberr),
#                                  mag_azim_aberr)
'''