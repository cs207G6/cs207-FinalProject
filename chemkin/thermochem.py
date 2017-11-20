import numpy as np


class ThermochemRXNSetWrapper:
    def __init__(self, nasa):
        self.nasa7_coeffs = nasa


class ThermoChem:
    """Methods for calculating the backward reaction rate.

    Cp_over_R: Returns specific heat of each specie given by
           the NASA polynomials.
    H_over_RT:  Returns the enthalpy of each specie given by
        the NASA polynomials.
    S_over_R: Returns the entropy of each specie given by
          the NASA polynomials.
    backward_coeffs:  Returns the backward reaction rate
              coefficient for reach reaction.

    Please see the notes in each routine for clarifications and
    warnings.  You will need to customize these methods (and
    likely the entire class) to suit your own code base.
    Nevertheless, it is hoped that you will find these methods
    to be of some use.
    """

    def __init__(self, rxnset, T):
        self.rxnset = rxnset
        self.p0 = 1.0e+05  # Pa
        self.R = 8.3144598  # J / mol / K
        self.h_rt = self.H_over_RT(T)
        self.s_rt = self.S_over_R(T)
        self.T = T

    def H_over_RT(self, T):
        # WARNING:  This line will depend on your own data structures!
        # Be careful to get the correct coefficients for the appropriate
        # temperature range.  That is, for T <= Tmid get the low temperature
        # range coeffs and for T > Tmid get the high temperature range coeffs.
        a = self.rxnset.nasa7_coeffs

        H_RT = (a[:, 0] + a[:, 1] * T / 2.0 + a[:, 2] * T ** 2.0 / 3.0
                + a[:, 3] * T ** 3.0 / 4.0 + a[:, 4] * T ** 4.0 / 5.0
                + a[:, 5] / T)

        return H_RT

    def S_over_R(self, T):
        # WARNING:  This line will depend on your own data structures!
        # Be careful to get the correct coefficients for the appropriate
        # temperature range.  That is, for T <= Tmid get the low temperature
        # range coeffs and for T > Tmid get the high temperature range coeffs.
        a = self.rxnset.nasa7_coeffs

        S_R = (a[:, 0] * np.log(T) + a[:, 1] * T + a[:, 2] * T ** 2.0 / 2.0
               + a[:, 3] * T ** 3.0 / 3.0 + a[:, 4] * T ** 4.0 / 4.0 + a[:, 6])

        return S_R

    def backward_coeffs(self, nuij, kf):
        gamma = np.sum(nuij, axis=0)

        # Change in enthalpy and entropy for each reaction

        assert (len(nuij) == len(self.h_rt))
        assert (len(nuij) == len(self.s_rt))

        delta_H_over_RT = np.dot(nuij, self.h_rt)
        delta_S_over_R = np.dot(nuij, self.s_rt)

        # Negative of change in Gibbs free energy for each reaction
        delta_G_over_RT = delta_S_over_R - delta_H_over_RT

        # Prefactor in Ke
        fact = self.p0 / self.R / self.T

        # Ke
        kb = fact ** gamma * np.exp(delta_G_over_RT)

        return kf / kb

        # def Cp_over_R(self, T):
        #     # WARNING:  This line will depend on your own data structures!
        #     # Be careful to get the correct coefficients for the appropriate
        #     # temperature range.  That is, for T <= Tmid get the low temperature
        #     # range coeffs and for T > Tmid get the high temperature range coeffs.
        #     a = self.rxnset.nasa7_coeffs
        #
        #     Cp_R = (a[:, 0] + a[:, 1] * T + a[:, 2] * T ** 2.0
        #             + a[:, 3] * T ** 3.0 + a[:, 4] * T ** 4.0)
        #
        #     return Cp_R
