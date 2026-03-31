# Piecewise-Deterministic Markov Processes
## Theory and Methods for Stationary Distributions

**Andres Ferragut — Universidad ORT Uruguay**

---

## 1. Motivation

Davis (1984) observed that continuous-time stochastic models consist of three ingredients: diffusion (Brownian noise), deterministic motion, and random jumps. While diffusion theory is highly unified via Itô calculus, non-diffusion models — combining deterministic motion and random jumps — were a heterogeneous collection of special cases. The goal of PDMPs is to provide a canonical class covering virtually all such models, with a unified stochastic calculus.

**Classical examples:** countable-state Markov processes, M/G/1 and GI/G/1 queues (virtual waiting time), dam and storage models.

**Modern applications:** TCP/IP congestion control, reliability and maintenance models, biochemical reaction networks, insurance ruin theory.

---

## 2. Definition of the PDMP

A PDMP on an open set $E \subseteq \mathbb{R}^d$ is determined by three local characteristics:

**(1) The flow $\varphi$.** Unique solution map of $\dot{x} = f(x)$, where $f$ is Lipschitz continuous:
$$\varphi(0,x) = x, \qquad \varphi(t+s,x) = \varphi(t,\varphi(s,x)).$$

**(2) The jump rate $\lambda$.** A measurable function $\lambda: E \to \mathbb{R}_+$. The **integrated hazard along the trajectory** is:
$$\Lambda(t,x) = \int_0^t \lambda(\varphi(s,x))\,ds.$$
This is the cumulative intensity seen by a particle starting at $x$ and following the ODE flow for time $t$. The first jump time satisfies:
$$P(T_1 > t \mid X_0 = x) = e^{-\Lambda(t,x)}.$$

**(3) The jump kernel $Q$.** A probability measure $Q(x,\cdot)$ on $E$. Post-jump locations $Z_i \sim Q(\varphi(T_i, x_0), \cdot)$ are drawn independently of the jump times.

### Reachable boundary and Davis's assumption

The **reachable boundary** is $\partial^* E = \{z \in \partial E : \varphi(-t,z) \in E\ \text{for all}\ t \in (0,\varepsilon)\}$ — those boundary points the flow can reach from the interior. The **boundary exit time** is $t^*(x) = \inf\{t > 0 : \varphi(t,x) \in \partial^* E\}$.

> **Davis's assumption.** If $t^*(x) < \infty$, then $\Lambda(t^*(x), x) = \infty$. Consequence: $P(T_1 < t^*(x)) = 1$, so the Poisson clock always fires before the boundary is reached. Forced boundary jumps are excluded. *Costa (1990) relaxes this assumption.*

### Sample path construction

Starting from $x_0 \in E$, one inter-jump cycle is constructed as follows:

1. **Draw the alarm level.** Sample $U \sim \text{Uniform}(0,1)$ and set $s^* = -\log U$.
2. **Run the flow and clock jointly.** Integrate $\dot{x} = f(x)$ from $x_0$ while accumulating $\Lambda(t, x_0)$. Stop at $T_1 = \inf\{t > 0 : \Lambda(t,x_0) \geq s^*\}$.
3. **Draw the post-jump location.** Sample $Z_1 \sim Q(\varphi(T_1, x_0), \cdot)$ independently of $T_1$.
4. **Restart** from $X_{T_1} = Z_1$ and repeat.

> **Key point.** Steps 1 and 2 cannot be separated: $T_1$ is the *output* of the joint integration, not an input to the evolution.

### Markov property and non-explosion

The semigroup property of $\varphi$ ensures the process is **strong Markov**. Non-explosion requires $E[N_t] < \infty$ for all $t$, equivalently $T_n \to \infty$ a.s. This must be verified for each model.

---

## 3. Example: TCP/IP Congestion Control

The TCP/IP AIMD model (additive increase, multiplicative decrease) is a natural PDMP.

| Object | Value |
|--------|-------|
| State space | $E = \mathbb{R}_+$, $x$ = sending rate |
| Flow | $f(x) = c > 0$, so $\varphi(t,x) = x + ct$ |
| Jump rate | $\lambda(x) = \alpha x$ |
| Jump kernel | $Q(x,\cdot) = \delta_{x/2}$ |

No reachable boundary ($t^*(x) = \infty$ for all $x$), so Davis's assumption holds trivially.

The generator for test functions $\psi \in \mathcal{D}(\mathcal{A})$ is:
$$\mathcal{A}\psi(x) = c\,\psi'(x) + \alpha x\left[\psi\!\left(\tfrac{x}{2}\right) - \psi(x)\right].$$

---

## 4. The Extended Generator

**Theorem (Davis 1984, Theorem 5.5).** For $\psi \in \mathcal{D}(\mathcal{A})$, the extended generator of the PDMP is:
$$\mathcal{A}\psi(x) = \underbrace{f(x)\cdot\nabla\psi(x)}_{\text{ODE drift}} + \underbrace{\lambda(x)\int_E[\psi(y)-\psi(x)]\,Q(x,dy)}_{\text{average jump contribution}}.$$

The domain $\mathcal{D}(\mathcal{A})$ consists of functions $\psi$ satisfying:
1. $t \mapsto \psi(\varphi(t,x))$ is absolutely continuous on $[0,t^*(x))$ for each $x \in E$.
2. **Boundary condition:** $\psi(x) = \int_E \psi(y)\,Q(x,dy)$ for $x \in \partial^* E$ (vacuous under Davis's assumption).
3. **Integrability:** $\int_E[\psi(y)-\psi(x)]\,Q(x,dy) \in L^{1,\text{loc}}(p)$.

### Dynkin formula

For $\psi \in \mathcal{D}(\mathcal{A})$, the process
$$C_t^\psi = \psi(X_t) - \psi(X_0) - \int_0^t \mathcal{A}\psi(X_s)\,ds$$
is a local $\mathcal{F}_t$-martingale. Note that $\mathcal{A}$ is **not a local operator** — the integral term couples $\psi(x)$ to all states reachable by a jump. Unlike diffusion, the calculus here is ordinary calculus with a stochastic interpretation — no Itô correction terms.

---

## 5. Kolmogorov Forward Equation

Let $p(x,t)$ denote the density of $X_t$. Taking the $L^2$-adjoint of $\mathcal{A}$:

$$\partial_t p(x,t) = \underbrace{-\nabla\cdot[f(x)\,p(x,t)]}_{\text{transport}} \underbrace{-\,\lambda(x)\,p(x,t)}_{\text{loss via jumps}} + \underbrace{\int_E \lambda(y)\,p(y,t)\,Q(y,dx)}_{\text{gain via jumps}}.$$

This is a **first-order hyperbolic PDE** (not Fokker-Planck — there is no diffusion term) coupled to a non-local integral operator. The gain term couples the density at all states $y$ jumping into $x$.

**Stationary condition** (set $\partial_t p = 0$): a density $\pi(x)$ is stationary if and only if:
$$\nabla\cdot[f(x)\,\pi(x)] + \lambda(x)\,\pi(x) = \int_E \lambda(y)\,\pi(y)\,Q(y,dx).$$

This is an integro-differential equation, generally hard to solve analytically.

### TCP/IP stationary equation

With $f(x) = c$, $\lambda(x) = \alpha x$, $Q(x,\cdot) = \delta_{x/2}$, changing variables $y = 2x$ in the gain term:

$$c\,\pi'(x) + \alpha x\,\pi(x) = 4\alpha x\,\pi(2x).$$

This couples $\pi(x)$ at scale $x$ to $\pi(2x)$ at scale $2x$ — a **functional differential equation**. Standard ODE methods do not apply directly. The natural tool is the **Mellin transform** (Baccelli et al., 2002).

---

## 6. Ergodicity

**Definition.** A PDMP with semigroup $(P_t)$ is **ergodic** if it admits a unique invariant probability measure $\pi$ satisfying $\int_E P_t\psi\,d\pi = \int_E \psi\,d\pi$ for all bounded measurable $\psi$ and $t \geq 0$, and $\mu_t \to \pi$ weakly for all initial conditions $x \in E$.

**Ergodic theorem.** Under ergodicity, for all $x \in E$ and $g \in L^1(\pi)$:
$$\frac{1}{T}\int_0^T g(X_t)\,dt \xrightarrow{T\to\infty} \int_E g\,d\pi \quad \text{a.s.}$$
This justifies Monte Carlo estimation of $\pi$ from a single trajectory.

**Stronger notions:**
- **Positive recurrence:** $E_x[\tau_C] < \infty$ for all $x$ and some petite set $C$.
- **Geometric ergodicity:** $\|\mu_t - \pi\|_V \leq C_0\,V(x_0)\,e^{-\varrho t}$ for some $V \geq 1$.

### Conditions for ergodicity

**Condition 1: $\nu$-irreducibility** (Meyn & Tweedie 1993, Ch. 4). The PDMP is $\nu$-irreducible with respect to a $\sigma$-finite measure $\nu$ if $\int_0^\infty P_t(x,A)\,dt > 0$ for all $x \in E$ and $\nu(A) > 0$. For a PDMP this requires: (i) ODE reachability — from $x$, the flow must reach a neighborhood of some $y$ from which a jump into $A$ is possible; and (ii) jump support — $Q(y,A) > 0$ for some such $y$.

> **Warning.** If the ODE has an invariant set $\mathcal{L}$ and $Q(x,\mathcal{L}) = 1$ for all $x \in \mathcal{L}$, the process is trapped and irreducibility fails. Irreducibility requires jumps to escape invariant sets of the flow.

Irreducibility is **necessary but not sufficient** for ergodicity.

**Condition 2: Harris recurrence** (Meyn & Tweedie 1993, Ch. 9). A $\nu$-irreducible PDMP is **Harris recurrent** if for every $A$ with $\nu(A) > 0$:
$$P_x\!\left(\int_0^\infty \mathbf{1}_{X_t \in A}\,dt = \infty\right) = 1 \quad \text{for all } x \in E.$$
Harris recurrence ensures every positive-measure set is visited infinitely often from every starting point.

**Condition 3: Petite sets** (Meyn & Tweedie 1993, Ch. 5). A set $C \subset E$ is **petite** if there exist $T > 0$, $\varepsilon > 0$ and a probability measure $\nu_C$ such that $P^T(x,\cdot) \geq \varepsilon\,\nu_C(\cdot)$ for all $x \in C$.

Every visit to $C$ carries a probability $\varepsilon$ of *refreshing* the distribution to $\nu_C$, breaking dependence on the initial condition. Without petite sets, two trajectories could visit $C$ infinitely often yet never couple. A compact set $C$ is petite for a PDMP if: $\lambda(x) \geq \lambda_{\min} > 0$ on $C$; $Q(x,\cdot)$ has a density bounded below on $C$; and the ODE is locally controllable within $C$.

**Condition 4: Positive recurrence** (Meyn & Tweedie 1993, Theorem 10.0.1). A Harris recurrent PDMP is **positively recurrent** if for some petite set $C$ with $\nu(C) > 0$: $E_x[\tau_C] < \infty$ for all $x \in E$, where $\tau_C = \inf\{t > 0 : X_t \in C\}$. If this holds for one petite set it holds for all. Moreover:
$$\pi(A) \propto E_x\!\left[\int_0^{\tau_C} \mathbf{1}_{X_t \in A}\,dt\right],$$
i.e. $\pi$ is proportional to the expected time spent in $A$ during one excursion from $C$.

> **Summary.** $\nu$-irreducibility + Harris recurrence + positive recurrence $\Rightarrow$ ergodicity.

### Sufficient conditions

**Foster's criterion** (Meyn & Tweedie 1993, Theorem 3.2; Costa & Dufour 1999, Theorem 3.1). Suppose the PDMP is $\nu$-irreducible. If there exist $V \geq 1$, a petite set $C$, and constants $\varepsilon > 0$, $b < \infty$ such that:
$$\mathcal{A}V(x) \leq -\varepsilon + b\,\mathbf{1}_C(x),$$
then the PDMP is Harris recurrent, positively recurrent, and admits a unique invariant distribution $\pi$. Applying the Dynkin formula and optional stopping:
$$E_x[\tau_C] \leq \frac{V(x) - 1 + b}{\varepsilon} < \infty.$$

**Lyapunov criterion** (Meyn & Tweedie 1993; Durmus, Guillin & Monmarché 2021). If instead $\mathcal{A}V(x) \leq -c\,V(x) + b\,\mathbf{1}_C(x)$ for some $c > 0$, then the PDMP is geometrically ergodic.

---

## 7. TCP/IP: Ergodicity Analysis

**Irreducibility.** The TCP/IP PDMP is $\nu$-irreducible with respect to Lebesgue measure. From any $x > 0$, the flow reaches any $y > x$ in time $(y-x)/c$; a jump then sends the process to $y/2$. By choosing $y = 2z$ for any target $z$, any interval $(a,b)$ is reachable in one flow segment plus one jump.

**Petite set.** For any $0 < a < b < \infty$, the compact set $C = [a,b]$ is petite. Since $Q(x,\cdot) = \delta_{x/2}$ is a point mass, one jump alone gives no absolutely continuous component. However, after the first jump the process lands at $x/2$ and flows for a random time before the second jump. Since the flow time varies continuously, the second post-jump location has an absolutely continuous distribution, giving the required minorization.

**Foster's criterion with $V(x) = x$.**
$$\mathcal{A}V(x) = c \cdot 1 + \alpha x\left[\frac{x}{2} - x\right] = c - \frac{\alpha x^2}{2}.$$
Choose $\delta > 0$ small and $C = [\delta,\, \sqrt{2c/\alpha} + \delta]$. For $x > \sqrt{2c/\alpha} + \delta$:
$$\mathcal{A}V(x) \leq c - \frac{\alpha}{2}\!\left(\sqrt{\frac{2c}{\alpha}} + \delta\right)^2 = -\alpha\delta\sqrt{\frac{c}{2\alpha}} - \frac{\alpha\delta^2}{2} =: -\varepsilon < 0.$$
Inside $C$: $\mathcal{A}V(x) \leq c =: b$. Foster's criterion holds, giving Harris recurrence, positive recurrence, and a unique stationary distribution $\pi$.

**Second moment from the generator.** Since $\pi$ is stationary, $\int \mathcal{A}\psi\,d\pi = 0$ for all $\psi \in \mathcal{D}(\mathcal{A})$. Choosing $\psi(x) = x$:
$$c - \frac{\alpha}{2}E_\pi[X^2] = 0 \implies \boxed{E_\pi[X^2] = \frac{2c}{\alpha}}.$$

**Bound on the mean.** By Jensen's inequality ($x^2$ is convex):
$$E_\pi[X]^2 \leq E_\pi[X^2] = \frac{2c}{\alpha} \implies E_\pi[X] \leq \sqrt{\frac{2c}{\alpha}}.$$
The exact value of $E_\pi[X]$ requires the Mellin transform (Baccelli et al., 2002).

---

## References

- M. H. A. Davis. *Piecewise-deterministic Markov processes: A general class of non-diffusion stochastic models.* J. R. Statist. Soc. B **46**(3):353–388, 1984.
- M. H. A. Davis. *Markov Models and Optimization.* Chapman & Hall, 1993.
- O. L. V. Costa. *Stationary distributions for piecewise-deterministic Markov processes.* J. Appl. Prob. **27**:60–73, 1990.
- R. L. Tweedie. *Sufficient conditions for ergodicity and recurrence of Markov chains on a general state space.* Stoch. Proc. Appl. **3**:385–403, 1975.
- S. P. Meyn, R. L. Tweedie. *Markov Chains and Stochastic Stability.* Springer-Verlag, 1993. (2nd ed.: Cambridge University Press, 2009.)
- S. P. Meyn, R. L. Tweedie. *Stability of Markovian processes III: Foster–Lyapunov criteria for continuous-time processes.* Adv. Appl. Prob. **25**(3):518–548, 1993.
- F. Dufour, O. L. V. Costa. *Stability of piecewise deterministic Markov processes.* SIAM J. Control Optim. **37**(5):1483–1502, 1999.
- A. Durmus, A. Guillin, P. Monmarché. *Piecewise deterministic Markov processes and their invariant measures.* Ann. Inst. H. Poincaré **57**(3), 2021.
- F. Baccelli, D. McDonald, J. Reynier. *A mean-field model for multiple TCP connections through a buffer implementing RED.* Performance Evaluation **49**(1–4):77–97, 2002.
