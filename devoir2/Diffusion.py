import numpy as np
import time
from IPython.display import clear_output



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

    def step(self, dt, order=0, S=1e-8, source_term=None):
        """
        Step the equation forward in time with a value of dt
        """
        A = self.assemble(dt)
        b = self.leftside(dt, order=order, S=S, source_term=source_term)
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

        elif self.scheme == 2:
            A[0, 0] = 1 + (2 * dt * self.Deff) / (self.dr ** 2)
            A[0, 1] = (-2 * dt * self.Deff) / (self.dr ** 2)

        A[-1, -1] = 1
        return A

    def leftside(self, dt, order=0, S=1e-8, source_term=None):
        """
        Compute the left side of the linear system
        order : order of the S term
        S : constant S term value
        """
        b = self.C_values

        if source_term is None:
            source_term = lambda r:np.zeros_like(r)

        b += source_term(self.R_values)

        if order == 0:
            b[:-1] -= dt*S

        elif order == 1:
            b[:-1] *= (1 - dt * self.k)
        
        return b

    def solve_for_n_nodes(self, n_nodes=None, scheme=None, order=None, 
                          step_size=1e8, max_solve_time=10.0, solve_tol=1e-8,
                          min_steps=10, source_term=None, verbose=False):
        """
        Solve the diffusion problem to convergence with a grid containing n nodes
        """
        # Define the differentiation scheme
        if scheme is None:
            scheme = self.scheme
        # Define the number of nodes
        if n_nodes is None:
            n_nodes = self.n_nodes
            self.C_values = np.zeros_like(self.R_values)
            self.C_values[-1] = self.Ce
        else:
            self.__init__(n_nodes, scheme=scheme)
            if verbose:
                print('Re-initialized the Diffusion class')
        if order is None:
            if hasattr(self, 'order'):
                order = self.order
            else:
                order=0
        # Convergence variables
        change = 100
        now = time.time()
        counter = 0
        self.history = []
        # Solve until the change is almost zero
        while True:
            # store the Concentration array
            #if verbose:
            #    print(self.C_values)
            C_vals = self.C_values.copy()
            #if counter % 10 == 0:
            self.history.append(C_vals)
            # Step forward in time
            _ = self.step(step_size, order=order, source_term=source_term)
            # Compute the change (L2 norm)
            change = np.sqrt(np.sum((C_vals - self.C_values)**2))
            # Verbose
            if verbose:
                print(f'Concentration at half domain : {self.C_values[n_nodes//2]}')
                print(f'Iteration : {counter}')
                print(f'Relative change : {100*change}')
            counter += 1
            # Check the time-out criterion
            if (time.time() - now) > max_solve_time :
                print(f'Timed out for n = {counter}') 
                break
            elif np.isclose(change, 0.0, atol=solve_tol):
                if counter >= min_steps:
                    break
            if verbose:
                clear_output(wait=True)


        return self.C_values, self.R_values

    @staticmethod
    def analytical_solution(r):
        """ Analytical solution for validation """
        C = 0.25 * 1e-8/1e-10  * 0.5**2 * ((r/0.5)**2 - 1) + 10
        return C
    
