import numpy as np
import time


class Diffusion(object):
    """
    Diffusion problem Finite Difference Method Solver
    """

    def __init__(self, n_nodes, Deff=1e-10, k=4e-9, Ce=10, R=0.5, scheme=1):
        """
        Constructor for the Diffusion problem solver
        1D equally spaced finite difference grid
        ______
        n_nodes: number of nodes in the domain
        Deff: Effective diffusion coefficient (m2/s)
        k: Reaction constant (4e-9 s-1)
        Ce: Dirichlet boundary value (10 mol/m3)
        R: Domain size (0.5 m)
        """
        # store the attributes
        self.n_nodes = n_nodes
        self.Deff = Deff
        self.k = k
        self.Ce = Ce
        self.R = R
        self.t = 0.0  # start at t = 0
        self.scheme = scheme
        # Create the grid and the grid size
        self.R_values = np.linspace(0, R, self.n_nodes)
        self.dr = self.R_values[1]
        # Initial values for C
        self.C_values = np.ones_like(self.R_values) * 0
        self.C_values[-1] = self.Ce
        # Set the system coefficients method
        if scheme == 1:
            self.system_coefficients = self.system_coefficients_1
        
        elif scheme == 2:
            self.system_coefficients = self.system_coefficients_2

    def step(self, dt, order=0, S=1e-8):
        """
        Step the equation forward in time with a value of dt
        """
        A = self.assemble(dt)
        b = self.leftside(dt, order=order, S=S)
        x = np.linalg.solve(A, b)
        self.t += dt
        self.C_values = x
        return self.t

    def system_coefficients_1(self, dt, r):
        """
        Coefficients of the linear equation system derived from the first
        differentiation scheme
        dt : time step (s)
        r : position in the domain (m)
        """
        C1 = (-dt*self.Deff)/(self.dr ** 2)
        C2 = 1.0 + (2 * dt * self.Deff)/(self.dr ** 2 ) + (dt * self.Deff)/(r * self.dr)
        C3 = (-dt * self.Deff)/(self.dr ** 2) - (dt * self.Deff) / (r * self.dr)
        return C1, C2, C3

    def system_coefficients_2(self, dt, r):
        """
        Coefficients of the linear equation system derived from the SECOND
        differentiation scheme
        dt : time step (s)
        r : position in the domain (m)
        """
        C1 = (dt * self.Deff) / (2 * r * self.dr) - (dt * self.Deff) / (self.dr ** 2)
        C2 = 1.0 + (2 * dt * self.Deff) / (self.dr ** 2)
        C3 = (-dt * self.Deff) / (self.dr ** 2) - (dt * self.Deff) / (2 * r * self.dr)
        return C1, C2, C3


    def assemble(self, dt):
        """
        Assemble the linear system for a timestep
        """
        # create the system matrix
        A = np.zeros((self.n_nodes, self.n_nodes))
        for i in np.arange(1, self.n_nodes-1):
            # Compute the t+1 terms coefficients
            r = self.R_values[i]
            # Construct the linear equation system
            C1, C2, C3 = self.system_coefficients(dt, r)
            A[i, i-1] = C1
            A[i, i] = C2
            A[i, i+1] = C3

        # Boundary conditions
        if self.scheme == 1:
            C1, C2, C3 = self.system_coefficients(dt, 1e-60)
            A[0, 0] = C2
            A[0, 1] = C1 + C3
            # A[0, 0] = 1.0 - (dt * self.Deff) / (1e-60 * self.dr)
            # A[0, 1] =  - (dt * self.Deff) / (1e-60 * self.dr)
            #A[0, 0] = 1.0 + (2 * dt * self.Deff) / (self.dr ** 2) + (dt * self.Deff) / (1e-60 * self.dr)
            #A[0, 1] = (-2.0 * dt * self.Deff) / (self.dr ** 2) - (dt * self.Deff) / (1e-60 * self.dr)

        elif self.scheme == 2:
            A[0, 0] = 1 + (2 * dt * self.Deff) / (self.dr ** 2)
            A[0, 1] = (-2 * dt * self.Deff) / (self.dr ** 2)

        A[-1, -1] = 1
        return A

    def leftside(self, dt, order=0, S=1e-8):
        """
        Compute the left side of the linear system
        order : order of the S term
        S : constant S term value
        """
        b = self.C_values

        if order == 0:
            b[:-1] -= dt*S

        elif order == 1:
            b[:-1] *= (1 - dt * self.k)
        
        return b

    def solve_for_n_nodes(self, n, scheme=None, order=0, step_size=1e8, max_solve_time=10.0):
        """
        Solve the diffusion problem with a grid containing n nodes
        """
        # Re instanciate the problem class
        self.__init__(n, scheme=scheme)
        # Convergence variables
        change = 100
        now = time.time()
        # Solve until the change is almost zero
        while ~np.isclose(change, 0):
            # store the Concentration array
            C_vals = self.C_values.copy()
            # Step forward in time
            _ = self.step(step_size, order=0)
            # Compute the change (L1 norm)
            change = np.sum(np.abs(C_vals - self.C_values))
            # Check the time-out criterion
            if (time.time() - now) > max_solve_time :
                print(f'Timed out for n = {n}')
                break

        return self.C_values, self.R_values

    @staticmethod
    def analytical_solution(r):
        """ Analytical solution for validation """
        C = 0.25 * 1e-8/1e-10  * 0.5**2 * ((r/0.5)**2 - 1) + 10
        return C
    
