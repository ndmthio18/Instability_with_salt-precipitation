"""Closed-form corrected law and an independent numerical BVP solve."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

import numpy as np
from scipy.integrate import quad, solve_bvp


@dataclass(frozen=True)
class Parameters:
    ceq: float = 0.264940869768
    rho: float = 1202.168577
    rho_cr: float = 2165.0
    eps: float = 0.43
    d0: float = 1.3e-9
    kr: float = 1.0e-3
    beta: float = 2.0
    grain: float = 25.0e-6
    j: float = 1.78e-5

    @property
    def deff(self) -> float:
        """Tortuosity-corrected diffusivity used by the ALE operator."""
        return np.sqrt(self.eps) * self.d0

    @property
    def av(self) -> float:
        return 36.0 * (1.0 - self.eps) / self.grain

    @property
    def mathcal_d(self) -> float:
        return self.eps * self.deff / (1.0 - self.ceq)

    @property
    def ell(self) -> float:
        return np.sqrt(self.mathcal_d / (self.av * self.kr))

    @property
    def ks(self) -> float:
        return self.beta * (1.0 - self.eps) * self.kr

    @property
    def da_star(self) -> float:
        return self.ks**2 / (self.mathcal_d * self.av * self.kr)

    @property
    def x(self) -> float:
        return np.sqrt(self.da_star)

    @property
    def s_ref(self) -> float:
        return self.ceq * self.j / (1.0 - self.ceq)

    def to_dict(self) -> dict[str, float]:
        out = asdict(self)
        out.update(
            deff=self.deff,
            av=self.av,
            mathcal_d=self.mathcal_d,
            ell=self.ell,
            ks=self.ks,
            da_star=self.da_star,
            x=self.x,
            s_ref=self.s_ref,
        )
        return out


def thickness_for_da(da: float, p: Parameters) -> float:
    if da <= 0:
        raise ValueError("Da must be positive")
    return p.ell * np.sqrt(da)


def analytical_profile(
    z: np.ndarray | Iterable[float], da: float, u: float, p: Parameters
) -> tuple[np.ndarray, np.ndarray]:
    """Return excess concentration c and dc/dz for the prescribed-flux BVP."""
    if u < 0:
        raise ValueError("u must be nonnegative")
    z = np.asarray(z, dtype=float)
    h = thickness_for_da(da, p)
    if z.min() < 0 or z.max() > h * (1.0 + 1e-12):
        raise ValueError("profile coordinates outside crust")
    hhat = np.sqrt(da)
    khat = p.rho * p.mathcal_d / p.ell
    b = u * p.s_ref / khat
    denominator = np.sinh(hhat) + p.x * np.cosh(hhat)
    a = (p.s_ref / khat) * (
        1.0 - u * (np.cosh(hhat) + p.x * np.sinh(hhat))
    ) / denominator
    zeta = z / p.ell
    c = a * np.cosh(zeta) + b * np.sinh(zeta)
    dcdz = (a * np.sinh(zeta) + b * np.cosh(zeta)) / p.ell
    return c, dcdz


def corrected_observables(da: float, u: float, p: Parameters) -> dict[str, float]:
    h = thickness_for_da(da, p)
    c, dcdz = analytical_profile(np.array([0.0, h]), da, u, p)
    screening = 1.0 - u / np.cosh(np.sqrt(da))
    g_top = p.s_ref * p.x * screening / (p.x + np.tanh(np.sqrt(da)))
    velocity = g_top / (p.rho_cr * (1.0 - p.eps))
    b_soil = u * p.ceq * p.j
    j_bottom = p.s_ref - p.rho * p.mathcal_d * dcdz[0]
    j_top_transport = p.s_ref - p.rho * p.mathcal_d * dcdz[-1]
    reaction_signed = j_bottom - j_top_transport
    return {
        "da": da,
        "u": u,
        "height_m": h,
        "bottom_excess": float(c[0]),
        "top_excess": float(c[-1]),
        "bottom_concentration": float(p.ceq + c[0]),
        "top_concentration": float(p.ceq + c[-1]),
        "b_soil_kg_m2_s": b_soil,
        "screening_ratio": screening,
        "top_flux_kg_m2_s": g_top,
        "velocity_m_s": velocity,
        "bottom_total_salt_flux_kg_m2_s": j_bottom,
        "top_transport_salt_flux_kg_m2_s": j_top_transport,
        "internal_signed_kg_m2_s": reaction_signed,
    }

def solve_independent_bvp(
    da: float, u: float, p: Parameters, points: int = 401, tol: float = 1e-10
) -> dict[str, object]:
    """Solve the linear BVP without using the closed-form profile as a guess."""
    h = thickness_for_da(da, p)
    hhat = np.sqrt(da)
    zeta = np.linspace(0.0, hhat, points)
    guess = np.zeros((2, points), dtype=float)

    def ode(_zeta: np.ndarray, y: np.ndarray) -> np.ndarray:
        return np.vstack((y[1], y[0]))

    b_soil = u * p.ceq * p.j

    def bc(ya: np.ndarray, yb: np.ndarray) -> np.ndarray:
        return np.array(
            [
                p.rho * p.eps * p.deff * ya[1] / p.ell - b_soil,
                p.s_ref
                - p.rho * p.mathcal_d * yb[1] / p.ell
                - p.rho * p.ks * yb[0],
            ]
        )

    solution = solve_bvp(ode, bc, zeta, guess, tol=tol, max_nodes=20000)
    if not solution.success:
        raise RuntimeError(f"solve_bvp failed: {solution.message}")
    y = solution.sol(zeta)
    z = zeta * p.ell
    c = y[0]
    dcdz = y[1] / p.ell
    reaction = p.rho * p.av * p.kr * c
    internal_signed = float(
        quad(
            lambda coordinate: p.rho
            * p.av
            * p.kr
            * float(solution.sol(coordinate)[0])
            * p.ell,
            0.0,
            hhat,
            epsabs=1e-16,
            epsrel=1e-12,
            limit=200,
        )[0]
    )
    g_top = float(p.rho * p.ks * c[-1])
    return {
        "z_m": z,
        "excess": c,
        "dcdz": dcdz,
        "reaction_kg_m3_s": reaction,
        "top_concentration": float(p.ceq + c[-1]),
        "bottom_concentration": float(p.ceq + c[0]),
        "top_flux_kg_m2_s": g_top,
        "velocity_m_s": g_top / (p.rho_cr * (1.0 - p.eps)),
        "internal_signed_kg_m2_s": internal_signed,
        "nodes": int(solution.x.size),
        "status": int(solution.status),
        "rms_residual_max": float(np.max(solution.rms_residuals)),
    }
