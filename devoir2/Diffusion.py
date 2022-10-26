import numpy as np
import time
from IPython.display import clear_output



class Diffusion(object):
    """
    Diffusion problem Finite Difference Method Solver
    """

    def __init__(self, n_nodes, Deff=1e-10, k=4e-9, Ce=10, R=0.5, default_source=True):
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
        self.default_source = default_source
        self.t = 0.0  # start at t = 0
        # Create the grid and the grid size
        self.R_values = np.linspace(0, R, self.n_nodes)
        self.dr = self.R_values[1]
        # Initial values for C
        self.C_values = np.zeros_like(self.R_values)
        self.C_values[-1] = self.Ce
        # Set the system coefficients method

    def step(self, dt, S=1e-8, source_term=None):
        """
        Step the equation forward in time with a value of dt
        """
        A = self.assemble(dt)
        b = self.leftside(dt, S=S, source_term=source_term)
        x = np.linalg.solve(A, b)
        self.t += dt
        self.C_values = x
        return self.t

    def system_coefficients(self, dt, r):
        """
        Coefficients of the linear equation system        
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
        # Only assemble nodes from i=1 to i=n-1
        for i in np.arange(1, self.n_nodes-1):
            # Compute the t+1 terms coefficients
            r = self.R_values[i]
            # Construct the linear equation system
            C1, C2, C3 = self.system_coefficients(dt, r)
            A[i, i-1] = C1
            A[i, i] = C2
            A[i, i+1] = C3

        # Gear scheme for the null flux condition at r=0
        A[0, 0:3] = -3, 4, -1

        # Cn = Ce
        A[-1, -1] = 1

        return A

    def leftside(self, dt, S=1e-8, source_term=None):
        """
        Compute the left side of the linear system
        S : constant S term value
        """
        # Get the actual concentration values
        b = self.C_values

        # If a source term is provided add it to the concentration values
        if source_term is not None:
            b[:-1] += source_term(self.R_values)[:-1]*dt

        if self.default_source:
            b[:-1] *= (1 - dt * self.k)
        
        # Null flux condition at r=0
        b[0] = 0

        return b


    def solve(self, step_size=1e8, max_solve_time=10.0, solve_tol=1e-8,
                    min_steps=10, source_term=None, verbose=False, min_time=None):
        """
        Solve the diffusion problem to convergence with a grid containing n nodes
        param step_size : time step size in seconds
        param max_solve_time : Max real time the solver can run in seconds
        param solve_tol : Solver relative tolerance for convergence of the L2 norm
        param min_steps : Minimum number of time steps the solver does
        param source_term : Source term callable to be added to the left hand side
        param verbose : if True display how things are going
        """
        # Convergence variables
        if min_time is None:
            min_time = step_size
        change = 100
        now = time.time()
        counter = 0
        self.history = []
        while True:
            # Solve until the change is almost zero
            while True:
                # store the Concentration array
                C_vals = self.C_values.copy()
                self.history.append(C_vals)
                # Step forward in time
                _ = self.step(step_size, source_term=source_term)
                # Compute the change (L2 norm)
                change = np.sqrt(np.sum((C_vals - self.C_values)**2))
                # Verbose
                if verbose:
                    clear_output(wait=True)
                    print(f'Time step size : {step_size} seconds')
                    print(f'Concentration at half domain : {self.C_values[self.n_nodes//2]}')
                    print(f'Iteration : {counter}')
                    print(f'Relative change : {100*change}')
                counter += 1
                # Check the time-out criterion
                if (time.time() - now) > max_solve_time :
                    print(f'Timed out for n = {counter}') 
                    break
                elif np.isclose(change, 0.0, rtol=solve_tol):
                    break
            # Step back
            self.C_values = self.history[-2]
            step_size *= 0.1
            solve_tol *= 0.1
            if step_size <= min_time:
                if counter >= min_steps:
                    break

        return self.C_values, self.R_values

    @staticmethod
    def analytical_solution(r):
        """ Analytical solution for validation """
        C = 0.25 * 1e-8/1e-10  * 0.5**2 * ((r/0.5)**2 - 1) + 10
        return C
    
