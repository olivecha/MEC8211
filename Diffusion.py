import numpy as np


class Diffusion(object):
    """
    Diffusion problem Finite Difference Method Solver
    """
    def __init__(self, n_nodes, Deff=10e-10, k=4e-9, Ce=10, R=0.5):
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
        # Create the grid and the grid size
        self.R_values = np.linspace(0, R, self.n_nodes)
        self.dr = self.R_values[1]
        # Initial values for C
        self.C_values = np.zeros_like(self.R_values)
        self.C_values[-1] = self.Ce

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

    def assemble(self, dt):
        """
        Assemble the linear system for a timestep
        """
        # create the system matrix
        A = np.zeros((self.n_nodes, self.n_nodes))
        coeff1 = (-dt*self.Deff)/self.dr
        for i in np.arange(1, self.n_nodes-1):
            # Compute the t+1 terms coefficients
            r = self.R_values[i]
            coeff2 = 1.0 + (2 * dt * self.Deff)/(self.dr ** 2 ) + (dt * self.Deff)/(r * self.dr)
            coeff3 = (-dt * self.Deff)/(self.dr ** 2) - (dt * self.Deff) / (r * self.dr)
            # Construct the linear equation system
            A[i, i-1] = coeff1
            A[i, i] = coeff2
            A[i, i+1] = coeff3

        # Boudary conditions
        A[0, 0] = 1
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
