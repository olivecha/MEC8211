import numpy as np
from diffusion import Diffusion
from scipy.interpolate import InterpolatedUnivariateSpline

class NearProblem(Diffusion):
    """
    Near Problem class derived from the Diffusion FDM class
    """

    def __init__(self, n_nodes, spline_k=5):
        """ Constructor for the Near Problem class"""

        # initialize the Diffusion problem
        super().__init__(n_nodes)

        # Solve the diffusion problem to obtain a reference solution
        C, R = self.solve(verbose=False)
        self.C_ref, self.R_ref = C.copy(), R.copy()
        # Define the spline and its derivatives
        self.spl, self.dspl, self.d2spl = None, None, None
        # fit the spline to the solution
        self.fit_solution(spline_k)

    def fit_solution(self, k):
        """ Fit the reference solution to a univariate spline """ 
        # Compute the fitted univariate spline 
        self.spl = InterpolatedUnivariateSpline(self.R_ref, self.C_ref, k=k)
        # Get the derivatives of the spline
        self.dspl = self.spl.derivative(1)
        self.d2spl = self.spl.derivative(2)

    def residual_source_term(self, r):
        """ Method returning the residual source therm from spline interpolation """
        r = np.array(r)
        # Sanitize close to zero values
        r[np.isclose(r, 0.0)] = 1e-40
        e = -r *  ((-self.Deff) * self.dspl(r) - r * self.Deff * self.d2spl(r) + r * self.k * self.spl(r))
        return e/r

    def solve_with_residual(self, **kwargs):
        """ Solve the Diffusion problem with the spline residual as a source term """
        kwargs['source_term'] = self.residual_source_term
        return self.solve(**kwargs)
    

class ManufacturedSol(Diffusion):
    """
    Method of Manufactured Solution class for validation derived from the Diffusion
    FDM class
    """

    def __init__(self, n_nodes, source_term):
        """ Constructor for the Manufactured Solution class"""

        # initialize the Diffusion problem
        super().__init__(n_nodes, default_source=False)
        self.source_term = source_term


    def solve_with_residual(self, **kwargs):
        """ Solve the Diffusion problem with an additional source term """
        kwargs['source_term'] = self.source_term
        return self.solve(**kwargs)
        


