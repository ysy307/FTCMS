# Mixed finite-element formulation for the freezing column

Implementation snapshot: 2026-09-21. This note describes the experimental mixed
column solver, not an assertion that every production FTCMS path uses this formulation.
The physical equations below use dimensional variables; implementation scaling is stated separately.

## 1. Unknowns, units, and constitutive maps

Let \(\Omega\subset\mathbb R^2\) have intrinsic coordinates \((x,z)\), with \(z\) upward.
The independent unknowns are

\[
\mathbf y=(T,\mathbf q_l,p_w),\qquad \mathbf g=(0,-g).
\]

| Symbol | Meaning | Unit |
| --- | --- | --- |
| \(T\) | Temperature, with a Celsius reference for the implemented enthalpy | degrees C |
| \(\mathbf q_l\) | Liquid volumetric Darcy flux | m/s |
| \(p_w\) | Independent pore-water pressure parameter of the phase closure | Pa |
| \(p_l=\mathcal P(T,p_w)\) | Liquid pressure returned by the freezing closure | Pa |
| \(W=\mathcal W(T,p_w)\) | Stored water equivalent per bulk volume | dimensionless |
| \(U=\mathcal U(T,p_w)\) | Enthalpy density | J/m3 |
| \(\lambda=\mathcal L(T,p_w)\) | Thermal conductivity | W/(m K) |
| \(\mathbf K_s\) | Saturated hydraulic conductivity tensor | m/s |
| \(k_r(T,p_w)\) | Relative conductivity returned by the hydraulic closure | dimensionless |
| \(\mathbf\Lambda=k_r\mathbf K_s/(\rho_0g)\) | Pressure mobility | m2/(Pa s) |
| \(s\) | Volumetric water source | 1/s |

The distinction between \(p_w\) and \(p_l\) is essential:
**the P0 unknown is \(p_w\), not a separately discretized P0 liquid pressure.**
Even when \(p_w\) is constant on a cell, \(p_l(T,p_w)\) can vary inside that cell.
Material indices are implicit in these constitutive maps.

The implemented storage is

\[
W=\theta_l+\frac{\rho_i}{\rho_w}\theta_i+\theta_v^\star+S_c.
\]

Here \(\theta_v^\star\) is the model's water-equivalent vapor content, not the gas volume fraction.
The compressive contribution \(S_c\) is supplied by the hydraulic model on its applicable saturated branch.
It must not be added to a Jacobian capacity while being omitted from the residual storage.
The frequently plotted quantity

\[
\theta_{\mathrm{eq}}=\theta_l+0.917\theta_i
\]

is a diagnostic for the present material constants; it is not the full definition of \(W\).
In particular, neither \(W_T=0\) nor \(W=\theta_{\mathrm{eq}}\) is assumed.

With the implemented phase properties, the enthalpy has the structure

\[
\begin{aligned}
U={}&(1-\phi)\rho_s c_sT+\rho_w c_w\theta_lT
 +\rho_i c_i\theta_iT-\rho_i L_f\theta_i\\
 &+\rho_w c_v\theta_v^\star T+\rho_wL_v\theta_v^\star.
\end{aligned}
\]

Properties and phase contents can themselves depend on state. Consequently,
\(U_T\) is not merely a fixed mixture heat capacity.

## 2. Continuous equations and cryosuction

The mixed column solves

\[
\begin{aligned}
\partial_t U-\nabla\cdot(\lambda\nabla T)&=0,\\
\mathbf\Lambda^{-1}\mathbf q_l+\nabla p_l-\rho_0\mathbf g&=\mathbf0,\\
\partial_tW+\nabla\cdot\mathbf q_l&=s.
\end{aligned}
\]

These are the equations of the implemented water-equivalent storage model.
Conservation claims below refer to this \(W\), rather than an independently reconstructed variable-density mass model.
Liquid heat advection is absent from the current mixed-column heat equation.
Vapor effects in the closure do not introduce an independent vapor-flux unknown here.

For intuition, the non-segregating freezing law contains the signed equilibrium suction

\[
\psi_f(T)=-\rho_w L_f\log\!\left(\frac{T_K}{T_{f,K}}\right),
\qquad T_K=T+273.15.
\]

The phase model combines matric and freezing effects through its effective-suction closure:

\[
p_l=-\psi_{\mathrm{eff}}(T,p_w).
\]

This notation does not replace the implemented smoothing or phase splitting by a hard maximum.
For a smooth continuous state, the explanatory chain rule is

\[
\nabla p_l=\mathcal P_T\nabla T+\mathcal P_p\nabla p_w.
\]

The mixed discretization uses \(p_l\) in an integrated pressure term; it does not compute
this expression by setting the elementwise gradient of P0 pressure to zero.
An additional thermal-suction flux must not be added on top of the same \(p_l\) driver.

For locally scalar mobility,

\[
q_z=-\Lambda\left(\partial_zp_l+\rho_0g\right).
\]

Upward flow therefore requires a sufficiently negative upward liquid-pressure gradient.
A large suction does not by itself imply a large flow: the mobility can be extremely small in ice-rich soil.

## 3. Function spaces and boundary conditions

Partition the hydraulic boundary into prescribed-flux \(\Gamma_q\) and prescribed-pressure \(\Gamma_p\).
Define

\[
\begin{aligned}
V_T&=\{\eta\in H^1(\Omega):\eta=0\text{ on thermal Dirichlet boundaries}\},\\
V_q&=\{\mathbf v\in H(\operatorname{div};\Omega):\mathbf v\cdot\mathbf n=0\text{ on }\Gamma_q\},\\
Q&=L^2(\Omega),\\
H(\operatorname{div};\Omega)&=\{\mathbf v\in[L^2(\Omega)]^2:
\nabla\cdot\mathbf v\in L^2(\Omega)\}.
\end{aligned}
\]

Trial fields satisfy the corresponding nonhomogeneous essential data.
The active Mizoguchi column has \(\mathbf q_l\cdot\mathbf n=0\) on the entire boundary.
A prescribed liquid pressure on \(\Gamma_p\) is included below to explain the general weak form;
it is not an extra boundary condition imposed on that closed column.

The heat-transfer boundary condition is

\[
-\lambda\nabla T\cdot\mathbf n=h(T-T_\infty).
\]

The present benchmark uses top \(h=28\), \(T_\infty=-6\), bottom \(h=3\), \(T_\infty=6.7\),
and insulated sides. This is a bottom-heat-leak configuration, not an insulated-bottom experiment.

## 4. Weak form and backward Euler

Let \(\Delta t_n=t^{n+1}-t^n>0\), and keep the accepted old state fixed during a nonlinear solve.
Write \((a,b)_\Omega=\int_\Omega ab\,d\Omega\), with the corresponding vector product when needed.
Find \((T^{n+1},\mathbf q_l^{n+1},p_w^{n+1})\) such that

\[
\begin{aligned}
R_T(\eta)={}&
\left(\frac{U^{n+1}-U^n}{\Delta t_n},\eta\right)_\Omega
 +(\lambda^{n+1}\nabla T^{n+1},\nabla\eta)_\Omega
 +\langle h(T^{n+1}-T_\infty),\eta\rangle_{\Gamma_h}=0,\\
R_q(\mathbf v)={}&
(\mathbf\Lambda^{-1}\mathbf q_l-\rho_0\mathbf g,\mathbf v)_\Omega
 -(p_l,\nabla\cdot\mathbf v)_\Omega
 +\langle\bar p_l,\mathbf v\cdot\mathbf n\rangle_{\Gamma_p}=0,\\
R_p(w)={}&
\left(\frac{W^{n+1}-W^n}{\Delta t_n},w\right)_\Omega
 +(\nabla\cdot\mathbf q_l^{n+1},w)_\Omega-(s,w)_\Omega=0.
\end{aligned}
\]

Unmarked constitutive quantities in the Darcy row are evaluated at \(n+1\).
The sign of the pressure term follows directly from

\[
(\nabla p_l,\mathbf v)_\Omega
=-(p_l,\nabla\cdot\mathbf v)_\Omega
 +\langle p_l,\mathbf v\cdot\mathbf n\rangle_{\partial\Omega}.
\]

Thus discontinuous pressure is admissible: its classical gradient need not exist cell-to-cell.

Backward Euler differences the stored quantities themselves:

\[
W^{n+1}-W^n=\mathcal W(T^{n+1},p_w^{n+1})-\mathcal W(T^n,p_w^n).
\]

Replacing this by a capacity times a pressure increment is generally a different residual.
First-order temporal accuracy requires suitable regularity and convergence;
implicit stability is not a large-step accuracy guarantee.

## 5. Discrete spaces

On a conforming mesh \(\mathcal T_h\), use

\[
V_{T,h}\subset H^1,\qquad V_{q,h}\subset H(\operatorname{div}),\qquad Q_h\subset L^2.
\]

| Cell | Temperature | Flux | Pressure | Local dimensions \((T,q,p)\) |
| --- | --- | --- | --- | --- |
| Triangle | continuous P1 | RT0 | discontinuous P0 | \((3,3,1)\) |
| Axis-aligned rectangle | continuous Q1 | rectangular RT0 | discontinuous Q0 | \((4,4,1)\) |

Q0 and P0 both mean a cell constant here. Current temperature elements are first order;
P2 temperature is not implied by the general \(H^1\) formulation.

For global basis functions \(N_a,\boldsymbol\phi_f,\chi_K\),

\[
T_h=\sum_a T_aN_a,\qquad
\mathbf q_h=\sum_f Q_f\boldsymbol\phi_f,\qquad
p_{w,h}=\sum_K P_K\chi_K.
\]

Choose one orientation \(\mathbf n_f\) per facet. The mathematical flux degree of freedom is

\[
Q_f=\int_f\mathbf q_h\cdot\mathbf n_f\,ds.
\]

PETSc basis normalization or scaled coefficients need not numerically equal this integral;
physical facet moments are reconstructed before balance checks.
With the same facet orientation on both sides,

\[
[\mathbf q_h\cdot\mathbf n_f]_f=0.
\]

For outward cell normals, this is instead
\(\mathbf q_{K^+}\cdot\mathbf n_{K^+}+\mathbf q_{K^-}\cdot\mathbf n_{K^-}=0\).
Tangential components need not be continuous.
These facet moments and the contravariant mapping are the defining RT construction
([DefElement](https://defelement.org/elements/raviart-thomas.html)).

## 6. Why triangular and rectangular RT0 differ

### 6.1 Triangle

On an affine triangle,

\[
\operatorname{RT}_0(K)=\{(a+bx,c+bz):a,b,c\in\mathbb R\},\qquad
\nabla\cdot\mathbf q_h=2b.
\]

The same coefficient controls both directional derivatives:

\[
\partial_xq_x=\partial_zq_z=b.
\]

A field \((0,dz+c)\) with \(d\ne0\) cannot belong to this space.
For nonzero divergence \(D\), the local horizontal component is \(a+Dx/2\).
For scalar inverse mobility \(r>0\), minimizing over \(a\) gives

\[
\min_a\int_Kr(a+Dx/2)^2\,dK
=\frac{D^2}{4}\int_Kr(x-\bar x_r)^2\,dK,\qquad
\bar x_r=\frac{\int_Krx\,dK}{\int_Kr\,dK}.
\]

This is a local representability argument, not a full nonlinear error estimate.
It explains how width and large inverse mobility can penalize a nominally vertical storage-changing flow.

### 6.2 Axis-aligned rectangle

On a rectangle,

\[
\operatorname{RT}_0(K)=Q_{1,0}\times Q_{0,1}
=\{(a+bx,c+dz):a,b,c,d\in\mathbb R\},\qquad
\nabla\cdot\mathbf q_h=b+d.
\]

Now \((0,c+dz)\) is representable. For a one-column closed strip,
zero normal flux at both vertical sides implies \(q_x\equiv0\) in every rectangle,
while \(\partial_zq_z\) can remain nonzero.

This removes that particular directional constraint. It does not eliminate vertical discretization error,
nonlinear error, underintegration, or temporal error.
The triangle-to-rectangle comparison also changes P1 temperature to Q1 and changes quadrature geometry.
Observed differences cannot all be attributed quantitatively to the flux space alone.

### 6.3 Explicit physical-cell basis and matrix

For \(K=[0,h_x]\times[0,h_z]\), an outward-unit-moment basis is

\[
\begin{aligned}
\boldsymbol\phi_L&=((x-h_x)/(h_xh_z),0),&
\boldsymbol\phi_R&=(x/(h_xh_z),0),\\
\boldsymbol\phi_B&=(0,(z-h_z)/(h_xh_z)),&
\boldsymbol\phi_T&=(0,z/(h_xh_z)).
\end{aligned}
\]

Each basis has unit outward integral on its named face, zero on the others, and divergence \(1/|K|\).
Consequently, with outward moments \(F_{Kf}\),

\[
\nabla\cdot\mathbf q_h|_K=\frac{F_{KL}+F_{KR}+F_{KB}+F_{KT}}{|K|}.
\]

For constant scalar inverse mobility \(r\), the element resistance matrix in this basis is

\[
A_K=\frac r6
\begin{pmatrix}
2h_x/h_z&-h_x/h_z&0&0\\
-h_x/h_z&2h_x/h_z&0&0\\
0&0&2h_z/h_x&-h_z/h_x\\
0&0&-h_z/h_x&2h_z/h_x
\end{pmatrix}.
\]

Local-to-global assembly applies the facet orientation signs to rows and columns.
For variable or anisotropic mobility, use the actual integral, not this constant-coefficient matrix.

## 7. Piola map, geometry, and divergence

Let \(F_K(\hat x)=J_K\hat x+b_K\) be an orientation-preserving affine map.
The physical RT field is

\[
\mathbf q(x)=\frac{1}{\det J_K}J_K\hat{\mathbf q}(\hat x),\qquad
\nabla_x\cdot\mathbf q=\frac{1}{\det J_K}\hat\nabla\cdot\hat{\mathbf q}.
\]

The mapped normal moment satisfies

\[
\int_f\mathbf q\cdot\mathbf n\,ds
=\int_{\hat f}\hat{\mathbf q}\cdot\hat{\mathbf n}\,d\hat s.
\]

The current column geometry gate accepts triangles and axis-aligned rectangles.
It does not claim support for arbitrary bilinear quadrilaterals.
A parallelogram is affine too; rejecting it is a narrower implementation policy, not a definition of affine geometry.
A warped quadrilateral has a varying Jacobian and needs matching geometry-aware reconstruction and integration.

## 8. Local and global conservation

Because the indicator \(\chi_K\) belongs to \(Q_h\), choose \(w=\chi_K\):

\[
\boxed{
\int_K(W_h^{n+1}-W_h^n)\,dK
+\Delta t_n\sum_{f\subset\partial K}F_{Kf}^{n+1}
-\Delta t_n\int_Ks\,dK=0.}
\]

This is the cell balance, not merely a global cancellation test.
Summing over cells cancels every internal normal moment and gives

\[
\int_\Omega(W_h^{n+1}-W_h^n)\,d\Omega
+\Delta t_n\int_{\partial\Omega}\mathbf q_h^{n+1}\cdot\mathbf n\,ds
-\Delta t_n\int_\Omega s\,d\Omega=0.
\]

In the closed, source-free column, the discrete stored-water inventory stays constant.
The implemented assertion is only as exact as the assembled residual solve and floating-point arithmetic.
Report the actual defect

\[
\epsilon_K=\mathcal Q_K[W^{n+1}-W^n]
+\Delta t_n\sum_fF_{Kf}-\Delta t_n\mathcal Q_K[s]
\]

alongside storage and flux magnitudes. Here \(\mathcal Q_K\) denotes the assembly quadrature.
The sum \(\sum_K\epsilon_K\) alone can hide opposite local errors.
In 2D the integrated water-equivalent inventory has units m2 per unit out-of-plane thickness;
facet transport has units m2/s.

## 9. Quadrature and constitutive history

For a reference-cell rule,

\[
\mathcal Q_K[f]=\sum_{\alpha=1}^{n_q}
\hat w_\alpha\det J_K(\hat x_\alpha)f(F_K(\hat x_\alpha)).
\]

All main fields and auxiliary history fields must use the same physical quadrature points and weights.
The column copies the temperature FE quadrature into the other main and auxiliary elements.
This includes volume and face rules; see the local PETSc `PetscFECopyQuadrature` implementation.
Matching only the number of points is insufficient.

At every point, evaluate

\[
U^n_\alpha=\mathcal U(T_h^n(x_\alpha),p_{w,h}^n(x_\alpha)),\qquad
W^n_\alpha=\mathcal W(T_h^n(x_\alpha),p_{w,h}^n(x_\alpha)).
\]

Old temperature retains P1/Q1 interpolation and old pressure retains P0 interpolation.
Neither old enthalpy nor old water is first replaced by a cell-centroid closure value.
A rejected step must leave these accepted-state fields unchanged.

Although pressure is P0, \(W(T_h,p_{w,h})\), \(p_l(T_h,p_{w,h})\), and \(k_r(T_h,p_{w,h})\)
are generally not polynomial cell constants. Thus

\[
|K|W(T_h(x_K),p_{w,h}|_K)\ne\mathcal Q_K[W]
\]

in general. A 1 cm slice aligned with cell boundaries does not make centroid water averages exact.
Nor does conservation of the quadrature inventory prove convergence of that inventory to the exact integral.

## 10. Newton linearization and blocks

Set \(r=\mathbf\Lambda^{-1}\), and let \(\delta y=(\delta T,\delta\mathbf q,\delta p_w)\).
At one nonlinear iterate, with fixed old state and fixed sources,

\[
\begin{aligned}
DR_T[\delta y](\eta)={}&
\left(\frac{U_T\delta T+U_p\delta p_w}{\Delta t_n},\eta\right)
 +(\lambda\nabla\delta T,\nabla\eta)\\
&+((\lambda_T\delta T+\lambda_p\delta p_w)\nabla T,\nabla\eta)
 +\langle h\delta T,\eta\rangle,\\
DR_q[\delta y](\mathbf v)={}&
(r\delta\mathbf q,\mathbf v)
 +((r_T\delta T+r_p\delta p_w)\mathbf q,\mathbf v)
 -(\mathcal P_T\delta T+\mathcal P_p\delta p_w,\nabla\cdot\mathbf v),\\
DR_p[\delta y](w)={}&
\left(\frac{W_T\delta T+W_p\delta p_w}{\Delta t_n},w\right)
 +(\nabla\cdot\delta\mathbf q,w).
\end{aligned}
\]

Boundary coefficients are fixed for this expression. The nonlinear update solves

\[
\begin{pmatrix}
J_{TT}&0&J_{Tp}\\
J_{qT}&J_{qq}&J_{qp}\\
J_{pT}&J_{pq}&J_{pp}
\end{pmatrix}
\begin{pmatrix}\delta\mathbf T\\\delta\mathbf Q\\\delta\mathbf P\end{pmatrix}
=-\begin{pmatrix}\mathbf R_T\\\mathbf R_q\\\mathbf R_p\end{pmatrix}.
\]

The zero \(J_{Tq}\) reflects the absence of heat advection in this implementation.
For fixed saturated tensor,

\[
r_T=-r\,\partial_T\log k_r,\qquad r_p=-r\,\partial_{p_w}\log k_r.
\]

The code registers these blocks explicitly; PETSc does not automatically differentiate the Fortran closure.
Liquid-pressure derivatives come from the phase model. Derivatives of water, enthalpy,
thermal conductivity, and log relative conductivity are obtained from branch-aware finite differences
of the actual residual quantities. At a constitutive switch this is a numerical/generalized tangent,
not a claim of a globally smooth exact Newton Jacobian.

### Isothermal linear limit and Schur complement

For \(p_l=p_w\), fixed \(r\), and nonnegative storage derivative, define

\[
A_{fg}=(r\boldsymbol\phi_g,\boldsymbol\phi_f),\quad
B_{Kf}=(\chi_K,\nabla\cdot\boldsymbol\phi_f),\quad
C_{KL}=\left(\frac{W_p}{\Delta t_n}\chi_L,\chi_K\right).
\]

The hydraulic block is

\[
\begin{pmatrix}A&-B^T\\B&C\end{pmatrix},\qquad
S=C+BA^{-1}B^T.
\]

Changing the sign of the conservation row gives a symmetric indefinite saddle-point form.
The full freezing Jacobian is generally nonsymmetric and has additional constitutive couplings;
its Schur complement is not automatically this simple positive operator.
In the incompressible, all-flux-boundary steady limit, pressure can have a constant nullspace;
a gauge or nullspace treatment is then needed. Do not impose a zero-mean gauge blindly on a
transient problem whose storage law fixes the absolute pressure parameter.

## 11. High contrast, stability, and limits of the guarantees

The resistance block is assembled as

\[
A_{fg}=\sum_K\int_K
\boldsymbol\phi_f^T\mathbf\Lambda^{-1}(x,T,p_w)\boldsymbol\phi_g\,dK.
\]

In the special one-dimensional steady, source-free scalar problem,

\[
p_l(H)-p_l(0)-\rho_0g_zH
=-q_z\int_0^H\Lambda^{-1}(z)\,dz,
\qquad
\Lambda_{\mathrm{eff}}=\frac{H}{\int_0^H\Lambda^{-1}(z)\,dz}.
\]

This recovers series resistance without manually replacing a general 2D tensor by a face harmonic average.
A general anisotropic or transient mixed system is not equivalent to a two-point flux formula.
Zero mobility is outside the inverse-mobility formulation and requires a separate impermeable-limit treatment.

For the standard mixed pair on an appropriate shape-regular affine mesh family,
stability is expressed by an inf-sup condition, with suitable pressure space/boundary normalization:

\[
\inf_{w_h\ne0}\sup_{\mathbf v_h\ne0}
\frac{(\nabla\cdot\mathbf v_h,w_h)}
{\|\mathbf v_h\|_{H(\mathrm{div})}\|w_h\|_{L^2}}
\ge\beta>0.
\]

One must not infer a contrast-independent condition number, monotone pressure,
nonlinear convergence, or an error bound independent of extreme element aspect ratio from this statement.
RT0 normal continuity and local conservation also do not by themselves prove a future displacement-pressure
Biot discretization locking-free; its full coupled spaces and parameter limits must be assessed separately.
No theoretical error bound for the complete nonlinear frozen-column implementation is established here.

## 12. Scaling and PETSc residual correspondence

The implementation stores scaled unknowns

\[
T=T_*\tilde T,\qquad p_w=p_*\tilde p_w,\qquad
\mathbf q=q_*\tilde{\mathbf q}.
\]

Coordinates remain physical. Rows are scaled as

\[
\tilde R_T=R_T/R_*,\qquad
\tilde R_q=R_q/p_*,\qquad
\tilde R_p=(L_*/q_*)R_p.
\]

The PETSc pointwise representation is \(R(v)=\int(v\cdot f_0+\nabla v:f_1)\).
In physical variables before row scaling:

| Row | \(f_0\) | \(f_1\) |
| --- | --- | --- |
| Heat | \((U-U^n)/\Delta t\) | \(\lambda\nabla T\) |
| Darcy | \(r\mathbf q-\rho_0\mathbf g\) | \(-p_lI\) |
| Water | \((W-W^n)/\Delta t+\nabla\cdot\mathbf q-s\) | \(0\) |

Scaling preserves residual zeros but changes the numerical meaning of absolute residual tolerances.
A tolerance proportional to \(\sqrt N\) has no universal physical interpretation without the row scaling and units.

## 13. Mesh and time-step comparisons

For the current rectangle of height \(H=0.2\) m and width \(L_x=0.02\) m,

\[
h_z=H/N_z,\quad h_x=L_x/N_x,\quad
N_K=N_xN_z,\quad N_v=(N_x+1)(N_z+1).
\]

| \(N_z\) | \(h_z\) [mm] | Cells at \(N_x=4\) |
| --- | --- | --- |
| 200 | 1 | 800 |
| 400 | 0.5 | 1600 |
| 800 | 0.25 | 3200 |
| 1600 | 0.125 | 6400 |
| 3200 | 0.0625 | 12800 |

Total facets are \(N_f=N_x(N_z+1)+(N_x+1)N_z\).
Before constraints the unknown count is \(N_v+N_f+N_K\).
Closed hydraulic boundaries remove \(2N_x+2N_z\) prescribed facet degrees of freedom.
Increasing \(N_z\) at fixed \(N_x\) also increases aspect ratio; formal uniform mesh-family estimates
must not be quoted without checking their geometry assumptions.

Different error sources should be kept conceptually separate:

\[
\text{observed error}\;\longleftarrow\;
\{\text{space},\text{time},\text{quadrature},\text{nonlinear/linear solve},\text{model/data}\}.
\]

This is a classification, not a proved additive error bound.
A maximum time step is not the actual accepted step, and iteration-count adaptivity is not a temporal-error estimator.
Different meshes can induce different accepted time grids even with identical control settings.
Changing intermediate output horizons can also change time grids when the integrator clips steps to hit them.
The requested runs use a maximum of 3600 s; time-step convergence validation is currently deferred by user choice.
Mesh comparisons should state this limitation rather than claim that all observed differences are purely spatial.

## 14. Source map and reading order

The equations above were checked against these local files. Paths are relative to this document.
Ongoing performance work can change implementation details; this note records the formulation rather than live run status.

- [Pointwise residuals and Jacobian blocks](../src/app/mixed_hydraulic/mixed_hydraulic_column_pointwise.F90)
- [Spaces, quadrature, boundary setup, time stepping, and reports](../src/app/mixed_hydraulic/mixed_hydraulic_column.F90)
- [Closure bridge and numerical tangents](../src/app/mixed_hydraulic/mixed_hydraulic_frozen.F90)
- [Flux reconstruction and balance audit](../src/app/mixed_hydraulic/mixed_hydraulic_balance.F90)
- [Hydraulic storage](../src/physics/governing/hydraulic/hydraulic_coefficients.F90)
- [Enthalpy](../src/physics/governing/thermal/thermal_coefficients.F90)
- [Freezing suction](../src/physics/constitutive/models/phase_change/liquid_solid/gcc_base.F90)
- [Cell tests and Mizoguchi driver](../test/test_mixed_hydraulic.F90)
- [Local PETSc quadrature-copy implementation](../third_party/petsc-src/src/dm/dt/fe/interface/fe.c)

For the mathematical RT family definition and mapping conventions, see
[DefElement: Raviart–Thomas](https://defelement.org/elements/raviart-thomas.html).
For a reader focused on cryosuction, start with sections 1–4 and 6;
for implementation review, continue with sections 8–12.
