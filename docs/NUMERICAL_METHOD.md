# Numerical method

## Original problem and scaling

The supported pencil is `K x + lambda G x = 0`, with real symmetric K/G and positive-definite K after constraints. G may be indefinite or singular. This sign convention is explicit in every bundle. A negative multiplier corresponds to a different load direction and is not turned into a positive multiplier by taking its absolute value.

With `D = diag(1/sqrt(diag(K)))`, solve the congruently scaled system `Ks = D K D`, `Gs = D G D`, and recover `x = D y`. Scaling does not change eigenvalues and is not artificial stabilization. Mechanisms and unsupported matrices are rejected, not repaired with small springs.

For small reference problems solve `-G y = mu K y` using the SPD K metric, then `lambda=1/mu` for positive reciprocal roots. Tiny reciprocal roots are discarded relative to the spectral scale; the dense reference routine is not a completeness certificate. Do not use `eigh(K,-G)` when -G is indefinite.

## Member predictor

Let P select all reduced equations supporting a member contribution. The elastic condensed stiffness is `(P.T K^-1 P)^-1`. Only selected inverse columns are computed through batched solves, sharing one elastic LU factorization. The full inverse is never formed.

Combine this stiffness with the target member's assembled G_i. This changes the surrounding frame to an elastic restraint and is therefore a predictor, not the original-frame eigenproblem. It guides the search only.

The alpha retains all member support equations and limits their count. The proposed interface-only condensation with interior member reconstruction is not implemented yet. There is no persistent cache across process runs or load cases. Large member supports require a reviewed shift override or a later interface-reduction implementation.

## Original-system searches

Predicted positive multipliers are clustered by a relative span. One sparse factorization of `Ks + sigma Gs` is used for each shared search. SciPy's buckling-mode Lanczos iteration uses A=Ks, M=-Gs, and OPinv=(Ks+sigma Gs)^-1. Seed vectors bias the iteration but do not impose a spatial constraint.

The shift is slightly detuned from the predictor to avoid an exact pole. Global random noise prevents an exactly invariant local starting subspace. Cluster-edge members without significant participation can trigger a singleton refinement within the remaining budget. Different shifts require different factorizations. This is not a guarantee of a fixed speed-up.

Native modes are first checked and reused. A member already represented by an energy candidate can skip further searching. This stopping rule does not guarantee that the candidate is LTB; therefore the result is explicitly not a lowest-LTB assessment.

## Eigenpair and attribution checks

An accepted numerical candidate must satisfy the original full-system equation within the scaled relative residual tolerance. Vectors are normalized by x.T K x. Member participation is `x_i.T K_i x_i / (x.T K x)`, using the member contribution rather than a global principal submatrix.

Observation operators can report member-local lateral and twist norms. Nonzero norms do not establish coupled LTB: noise, rigid transport, mixed modes and different physical units require a more sophisticated kinematic classifier. The current report deliberately uses MEMBER_MODE_CANDIDATE, not 'LTB verified'.

Nearly identical eigenvectors are deduplicated using their K-inner product. Independent vectors at the same eigenvalue are retained. However, member attribution across a repeated eigenspace is still basis-dependent; a basis-invariant subspace classifier is an open task.

## Limits and failure states

No Sturm/inertia completeness certification, general interval exhaustion, nonlinear equilibrium continuation, code resistance, or automatic Mcr conversion is provided. Residual validity, physical interpretation, lowest-mode coverage and structural resistance are separate questions.

Configured caps bound number of shifts, iteration count, local predictor size and RHS batch size. Time limits inside the core are soft and checked between operations. Sparse factorization fill and retained mode vectors can still dominate memory. The Mechanical worker watchdog adds sampled process-level limits, but does not make a difficult full-frame factorization intrinsically cheap.

A given load pattern may have genuinely coupled member modes rather than an independent LTB eigenvector for each selection. UNRESOLVED and CANDIDATE_REQUIRES_REVIEW are expected outcomes, not conditions to hide or convert into a green check.

Reference: [SciPy eigsh buckling mode](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html).
