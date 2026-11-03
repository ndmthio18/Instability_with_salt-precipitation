"""Reference algebra for a conservative equal-height constraint.

This adapter is tested against synthetic residuals; it is not yet wired to
DuMuX FVAssembler/NewtonSolver and must not be advertised as a validated solver.
"""
import numpy as np

class EqualHeightConstraint:
    def __init__(self, heights, components=2):
        self.heights=np.asarray(heights);self.levels,self.index=np.unique(self.heights,return_inverse=True);self.components=components
        self.P=np.zeros((len(heights)*components,len(self.levels)*components))
        for i,j in enumerate(self.index):
            for c in range(components):self.P[i*components+c,j*components+c]=1
    def lift(self,u):return (self.P@np.asarray(u).ravel()).reshape(len(self.heights),self.components)
    def residual(self,u,full_residual):return self.P.T@np.asarray(full_residual(self.lift(u))).ravel()
    def jacobian(self,j):return self.P.T@j@self.P
    def project_initial(self,full):
        # Only accept already uniform initial data. No averaging an evolved state.
        a=np.asarray(full);u=np.array([a[self.index==j][0] for j in range(len(self.levels))])
        if not np.array_equal(self.lift(u),a):raise ValueError('initial state not exactly uniform')
        return u

    def solve(self, initial, assemble, fixed=None, atol=1e-12, maxiter=20):
        """Reduced Newton with exact prescribed variables and no post-step averaging.

        assemble(full_state) must return integrated full residual and its Jacobian,
        for ONE fixed accepted old state and ONE requested dt. Constitutive and
        flux evaluations remain the responsibility of the caller. This dense
        reference implementation is for verification, not a production backend.
        """
        u=self.project_initial(initial).ravel().copy();fixed=fixed or {}
        for i,value in fixed.items():u[i]=value
        for iteration in range(maxiter+1):
            full_r,full_j=assemble(self.lift(u));r=self.P.T@np.asarray(full_r).ravel();j=self.jacobian(full_j)
            for i,value in fixed.items():r[i]=u[i]-value;j[i,:]=0;j[i,i]=1
            if np.max(abs(r))<=atol:return self.lift(u),dict(iterations=iteration,reduced_residual_Linf=float(np.max(abs(r))))
            if iteration==maxiter:raise RuntimeError('reduced Newton did not converge; no accepted step')
            u-=np.linalg.solve(j,r)
            for i,value in fixed.items():u[i]=value
