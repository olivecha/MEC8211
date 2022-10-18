import numpy as np
from diffusion import Diffusion
from scipy.interpolate import InterpolatedUnivariateSpline

class NearProblem(Diffusion):
    """
    Near Problem class derived from the Diffusion FDM class
    """

    def __init__(self, n_nodes, scheme=2, order=1):
        """ Constructor for the Near Problem class"""

        # initialize the Diffusion problem
        super().__init__(n_nodes, scheme=scheme)
        self.order=order

        # Solve the diffusion problem to obtain a reference solution
        C, R = self.solve_for_n_nodes()
        self.C_ref, self.R_ref = C.copy(), R.copy()
        # Fit a spline to the solution
        self.spl, self.dspl, self.d2spl = None, None, None
        self.fit_solution()

    def fit_solution(self, k=5):
        """ Fit the reference solution to a univariate spline """ 
        # Compute the fitted univariate spline 
        self.spl = InterpolatedUnivariateSpline(self.R_ref.copy(), self.C_ref.copy(), k=k)
        # Get the derivatives of the spline
        self.dspl = self.spl.derivative(1)
        self.d2spl = self.spl.derivative(2)

    def residual_source_term(self, r):
        """ Method returning the residual source therm from spline interpolation """
        r = np.array(r)
        r[np.isclose(r, 0.0)] = 1e-40
        e = -r *  ((-self.Deff) * self.dspl(r) - r * self.Deff * self.d2spl(r) + r * self.k * self.spl(r))
        return e/r

    def solve_with_residual(self, **kwargs):
        kwargs['source_term'] = self.residual_source_term
        return self.solve_for_n_nodes(**kwargs)
    

        


