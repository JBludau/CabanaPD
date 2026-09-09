---
icon: lucide/sigma
---

# Theoretical framework

CabanaPD uses the peridynamic formulation[^1] [^2]. Given a bounded body $\mathcal{B}\subset \mathbb{R}^3$, the equation of motion for a material point ${\bf x} \in \mathcal{B}$ at time $t\geqslant 0$ is:

$$
\rho
\frac{\partial^2 {\bf u}}{\partial t^2}({\bf x},t)
= \int_{\mathcal{H}_{\bf x}} {\bf f}({\bf x}^\prime, {\bf x},t) dV_{\bf x^\prime} + {\bf b}({\bf x},t),
$$

where $\mathbf{x}$ denotes the position in the reference configuration, $\rho$ is the mass density, ${\bf u}$ is the displacement field, ${\bf f}$ is the bond force function, ${\bf b}$ is a prescribed body force density, and $\mathcal{H}_{\bf x}$ is the nonlocal neighborhood of ${\bf x}$ defined by
$$
\mathcal{H}_{\bf x}:=\lbrace {\bf x}^\prime \in \mathcal{B}: \lVert{\bf x}^\prime -{\bf x}\rVert\leqslant \delta \rbrace,
$$
where $\|\cdot\|$ denotes the Euclidean norm, $\delta>0$ is the maximum distance over which nonlocal interactions occur (called the horizon), and ${\bf x}^\prime\in\mathcal{H}_{\bf{x}}$ denotes a neighboring material point.

## References

[^1]: S.A. Silling, Reformulation of elasticity theory for discontinuities and long-range forces, *Journal of the Mechanics and Physics of Solids* **48**(1) (2000): 175-209. [doi:10.1016/S0022-5096(99)00029-0](https://doi.org/10.1016/S0022-5096(99)00029-0)

[^2]: S.A. Silling, M. Epton, O. Weckner, J. Xu, and E. Askari, Peridynamic states and constitutive modeling, *Journal of Elasticity* **88** (2007): 151–184 (2007). [doi:10.1007/s10659-007-9125-1](https://doi.org/10.1007/s10659-007-9125-1)
