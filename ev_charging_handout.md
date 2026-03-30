# Stationary Analysis of a Battery-Buffered EV Charging Station
## A Piecewise-Deterministic Approach — Lecture Handout

---

## 1. Model

### 1.1 System Description

We study a single charging station with one charger, a battery buffer of capacity $B \in (0, \infty]$, and a grid connection. The station operates as an **Erlang loss system**: a car arriving when the charger is busy is turned away.

**Parameters**

| Symbol | Meaning | Units |
|--------|---------|-------|
| $\lambda$ | Poisson arrival rate of cars | cars/hr |
| $Y \sim F$ | Energy demand per car, mean $1/\mu$ | kWh |
| $p_0$ | Battery charging rate from grid (idle) | kW |
| $p_1 > p_0$ | Grid supply rate when serving | kW |
| $p_2 > p_1$ | Fast charging rate (grid + battery) | kW |
| $c = p_2 - p_1$ | Net battery discharge rate | kW |
| $B$ | Battery capacity | kWh |

**State:** $(n(t),\, x(t))$ where $n \in \{0,1\}$ is the number of cars present and $x \in [0, B]$ is the battery charge level.

### 1.2 Two-Phase Service

When a car is present, service proceeds in up to two phases depending on the battery level:

- **Phase 1 — Fast charging** ($x > 0$): the battery is non-empty; the charger delivers power at rate $p_2$; the battery depletes at rate $c$, so $\dot{x} = -c$.
- **Phase 2 — Slow charging** ($x = 0$): the battery is empty; the charger delivers only the grid rate $p_1 < p_2$; the battery remains at zero, $\dot{x} = 0$.

> **Key observation.** If $Y \sim \mathrm{Exp}(\mu)$, the memoryless property implies the residual demand at any time is again $\mathrm{Exp}(\mu)$. The departure rate is therefore **Markovian** and depends only on the current state $(n, x)$:
> $$\delta(x) = \begin{cases} \mu p_2 & x > 0 \\ \mu p_1 & x = 0 \end{cases}$$
> For a **general** distribution $F$, the residual demand must be tracked explicitly — see Section 5.

### 1.3 The PDMP Formulation

The process $Z(t) = (n(t), x(t))$ on $E = \{0,1\} \times [0,B]$ is a **piecewise-deterministic Markov process** (Davis, 1984) specified by three objects:

**Flow** $\dot{x} = r(n,x)$:
$$r(n,x) = \begin{cases} +p_0 & n=0,\; x < B \\ 0 & n=0,\; x = B \\ -c & n=1,\; x > 0 \\ 0 & n=1,\; x = 0 \end{cases}$$

**Jump rate** $\lambda(n,x)$:
$$\lambda(n,x) = \begin{cases} \lambda & n = 0 \\ \mu p_2 & n=1,\; x > 0 \\ \mu p_1 & n=1,\; x = 0 \end{cases}$$

**Transition kernel** $Q$:
$$Q\bigl((n,x),\,\cdot\bigr) = \delta_{(1-n,\; x)}$$
Only the discrete component $n$ flips at each jump; $x$ is unchanged.

**Boundary behavior:**
- $x = 0$, $n = 1$: reflecting boundary — slow-charging mode, departure rate switches to $\mu p_1$.
- $x = B$, $n = 0$: reflecting boundary — battery full, $\dot{x} = 0$.
- No probability mass accumulates at $(n=0, x=0)$ or $(n=1, x=B)$.

---

## 2. Exponential Demands, $B = \infty$

### 2.1 Stationary Distribution Structure

For $B = \infty$ the state space is $E = \{0,1\} \times [0,\infty)$. There is no upper boundary and no point mass $\pi_0^{(B)}$. The stationary distribution consists of:
- densities $f_0(x)$, $f_1(x)$ on $(0,\infty)$
- a point mass $\pi_1^{(0)} = \Pr(n=1, x=0)$ at the slow-charging boundary.

> **Stability is not automatic.** For $B < \infty$ the state space is compact and a stationary distribution always exists. For $B = \infty$ the battery can drift to $+\infty$.

### 2.2 Stability Condition

Define the **key eigenvalue**:
$$\theta_1 := \frac{\mu p_2}{c} - \frac{\lambda}{p_0}$$

**Theorem (Stability, $B=\infty$).** *A proper stationary distribution exists if and only if*
$$\boxed{\theta_1 < 0 \quad \iff \quad \lambda > \lambda_c := \frac{\mu p_2 p_0}{c}}$$

**Physical interpretation.** The condition $\lambda > \lambda_c$ says that the mean energy drained per busy period ($c/(\mu p_2)$) exceeds the mean energy added per idle period ($p_0/\lambda$). The station must be busy enough to drain what the grid puts in.

**Regime diagram.** The sign of $\theta_1$ determines the shape of the stationary battery density:

| $\theta_1 > 0$ | $\theta_1 = 0$ | $\theta_1 < 0$ |
|---|---|---|
| Density increasing toward $B$ | Uniform density | Density decreasing toward 0 |
| Battery tends to be **full** | Critical balance | Battery tends to be **empty** |
| Fast charging almost always available | — | Slow charging prevalent |

### 2.2a Proof of the Stability Condition via Foster–Lyapunov

We prove that $\theta_1 < 0$ is sufficient for positive recurrence by constructing an explicit Lyapunov function for the generator $\mathcal{L}$ of the PDMP. Recall that $\mathcal{L}$ acts on functions $V(n,x)$ as:
$$\mathcal{L}V(n,x) = r(n,x)\,\partial_x V(n,x) + \lambda(n,x)\,[V(1-n,x) - V(n,x)]$$
where the first term is the flow part and the second is the jump part.

The **Foster–Lyapunov criterion** says: if there exists $V \geq 0$ with $V \to \infty$, a petite set $C$ (here, any compact set suffices), and constants $\varepsilon > 0$, $b < \infty$ such that
$$\mathcal{L}V(n,x) \leq -\varepsilon + b\,\mathbf{1}_C(n,x)$$
then the process is positive recurrent.

**First attempt: $V(n,x) = x$.** The jump part vanishes since $x$ is unchanged at jumps. The flow part gives:
$$\mathcal{L}V(0,x) = p_0 > 0, \qquad \mathcal{L}V(1,x) = -c < 0.$$
This fails: the drift is positive in state $n=0$ for all $x$, including arbitrarily large $x$ outside any compact set.

**Second attempt: $V(n,x) = x + \alpha n$**, $\alpha > 0$. Now the jump part contributes $\lambda(n,x)\cdot\alpha(1-2n)$:
$$\mathcal{L}V(0,x) = p_0 + \lambda\alpha > 0, \qquad \mathcal{L}V(1,x) = -c - \mu p_2\alpha < 0.$$
Still fails: no choice of $\alpha$ makes both values simultaneously negative.

**The fix: add a corrector $h_n$.** Try $V(n,x) = x + h_n$ where $h_n$ is a constant depending only on $n$. Since $\partial_x V = 1$ and the jump part now contributes $\lambda(n,x)(h_{1-n} - h_n)$:
$$\mathcal{L}V(0,x) = p_0 + \lambda(h_1 - h_0)$$
$$\mathcal{L}V(1,x) = -c + \mu p_2(h_0 - h_1)$$

Both are constants (independent of $x$). Set $d := h_0 - h_1$ and require both expressions to equal the same constant $-\varepsilon$:
$$p_0 - \lambda d = -\varepsilon \quad \text{and} \quad -c + \mu p_2 d = -\varepsilon.$$
Equating gives:
$$p_0 - \lambda d = -c + \mu p_2 d \quad\Longrightarrow\quad d = \frac{p_0 + c}{\lambda + \mu p_2}$$

Substituting back:
$$-\varepsilon = p_0 - \lambda\cdot\frac{p_0+c}{\lambda+\mu p_2} = \frac{p_0(\lambda+\mu p_2) - \lambda(p_0+c)}{\lambda+\mu p_2} = \frac{p_0\mu p_2 - \lambda c}{\lambda+\mu p_2}$$

This is negative — i.e., $\varepsilon > 0$ — if and only if $\lambda c > p_0\mu p_2$, which is exactly $\theta_1 < 0$.

**Theorem (Foster–Lyapunov Certificate).** *Take $h_0 = (p_0+c)/(\lambda+\mu p_2)$ and $h_1 = 0$, so that $V(n,x) = x + h_n$. Then for all $x > 0$:*
$$\mathcal{L}V(n,x) = \frac{p_0\mu p_2 - \lambda c}{\lambda+\mu p_2} = -\varepsilon$$
*uniformly in both states $n \in \{0,1\}$. When $\theta_1 < 0$ this constant is negative, and the Foster–Lyapunov criterion with any compact petite set $C$ establishes positive recurrence.* $\square$

**Verification.** Let us check both states explicitly. With $d = h_0 - h_1 = (p_0+c)/(\lambda+\mu p_2)$:

State $n=0$:
$$\mathcal{L}V(0,x) = p_0 - \lambda d = p_0 - \frac{\lambda(p_0+c)}{\lambda+\mu p_2} = \frac{p_0\mu p_2 - \lambda c}{\lambda+\mu p_2} = -\varepsilon \checkmark$$

State $n=1$, $x>0$:
$$\mathcal{L}V(1,x) = -c + \mu p_2 d = -c + \frac{\mu p_2(p_0+c)}{\lambda+\mu p_2} = \frac{-c(\lambda+\mu p_2)+\mu p_2(p_0+c)}{\lambda+\mu p_2} = \frac{p_0\mu p_2 - \lambda c}{\lambda+\mu p_2} = -\varepsilon \checkmark$$

Both states give exactly $-\varepsilon$, with no remainder, for all $x > 0$.

**Remark.** The corrector $h_n$ has a clear interpretation: it makes $V$ slightly larger in state $n=0$ (where $\dot{x} = p_0 > 0$ pushes $V$ upward) than in state $n=1$ (where $\dot{x} = -c < 0$ pulls $V$ downward). The amount $d = (p_0+c)/(\lambda+\mu p_2)$ is exactly what is needed to equalize the drift across states. The common drift value $-\varepsilon$ is the **weighted average** of the per-state drifts $p_0$ and $-c$ under the stationary distribution $(\pi_0, \pi_1) = (\mu p_2, \lambda)/(\lambda+\mu p_2)$ of the fast chain:
$$-\varepsilon = \pi_0 \cdot p_0 + \pi_1 \cdot (-c) = \frac{p_0\mu p_2 - \lambda c}{\lambda+\mu p_2}$$

### 2.3 Balance Equations

For $x > 0$:
$$p_0\, f_0'(x) = -\lambda\, f_0(x) + \mu p_2\, f_1(x) \tag{E0}$$
$$-c\, f_1'(x) = \phantom{-}\lambda\, f_0(x) - \mu p_2\, f_1(x) \tag{E1}$$

At $x = 0$:
$$p_0\, f_0(0^+) = \mu p_1\, \pi_1^{(0)} = c\, f_1(0^+)$$

In matrix form $\mathbf{f}' = A\mathbf{f}$ with $A = \mathbf{D}^{-1}Q$, $\mathbf{D} = \mathrm{diag}(p_0, -c)$:
$$A = \begin{pmatrix} -\lambda/p_0 & \mu p_2/p_0 \\ \lambda/c & -\mu p_2/c \end{pmatrix}$$

The characteristic equation is $\theta(\theta - \theta_1) = 0$, with eigenvalues $0$ and $\theta_1$.

### 2.4 Solving ($B = \infty$, $\theta_1 < 0$)

The general solution is:
$$\begin{pmatrix}f_0(x)\\f_1(x)\end{pmatrix} = \alpha_0\begin{pmatrix}\mu p_2\\\lambda\end{pmatrix} + \alpha_1\begin{pmatrix}c\\p_0\end{pmatrix}e^{\theta_1 x}$$

For $B = \infty$ with $\theta_1 < 0$, normalizability on $[0,\infty)$ requires $\alpha_0 = 0$ (the constant mode is not integrable).

**Theorem (Stationary Distribution, $B=\infty$).** *Under $\theta_1 < 0$, with normalization constant $K_\infty = (p_0+c)/|\theta_1| + p_0 c/(\mu p_1)$:*
$$f_0(x) = \frac{c}{K_\infty}e^{\theta_1 x}, \qquad f_1(x) = \frac{p_0}{K_\infty}e^{\theta_1 x}, \qquad \pi_1^{(0)} = \frac{p_0 c}{\mu p_1 K_\infty}$$

**Corollary (Proportionality / Decoupling).** *For all $x > 0$:*
$$\frac{f_1(x)}{f_0(x)} = \frac{p_0}{c} = \text{const}$$

*The conditional probability of a car being present, given battery level $x > 0$, is independent of $x$:*
$$\Pr(n=1 \mid x) = \frac{p_0}{p_0+c}, \qquad \Pr(n=0 \mid x) = \frac{c}{p_0+c}$$

The full distribution is characterized by the single marginal battery density $f(x) = f_0(x) + f_1(x) = \frac{p_0+c}{K_\infty}e^{\theta_1 x}$.

---

## 3. Exponential Demands, $B < \infty$

### 3.1 What Changes

For finite $B$, the upper boundary $x = B$ is reflecting when $n = 0$, producing a new point mass:
$$\pi_0^{(B)} = \Pr(n=0,\; x=B)$$

**Additional boundary conditions at $x = B$:**
$$\lambda\,\pi_0^{(B)} = p_0\,f_0(B^-), \qquad c\,f_1(B^-) = \lambda\,\pi_0^{(B)}$$

The interior balance equations (E0)–(E1) are identical to the $B=\infty$ case. The boundary conditions now determine both $\alpha_0$ and $\alpha_1$.

### 3.2 Stationary Distribution

The flux condition at $x = B$ gives $\alpha_0(c\lambda - p_0\mu p_2) = 0$. Since $\theta_1 \neq 0$ implies $c\lambda \neq p_0\mu p_2$, we again conclude $\alpha_0 = 0$.

**Theorem (Stationary Distribution, $B < \infty$, $\theta_1 \neq 0$).**
$$f_0(x) = \frac{c}{K}e^{\theta_1 x}, \qquad f_1(x) = \frac{p_0}{K}e^{\theta_1 x}, \qquad x \in (0,B)$$
$$\pi_1^{(0)} = \frac{p_0 c}{\mu p_1 K}, \qquad \pi_0^{(B)} = \frac{p_0 c}{\lambda K}e^{\theta_1 B}$$

where
$$K = \frac{p_0+c}{\theta_1}\bigl(e^{\theta_1 B}-1\bigr) + \frac{p_0 c}{\mu p_1} + \frac{p_0 c}{\lambda}e^{\theta_1 B}$$

> **Remark.** The single-exponential form and proportionality $f_1/f_0 = p_0/c$ hold for *all* $\theta_1 \neq 0$, whether positive or negative. The $B = \infty$ case is recovered as $B \to \infty$ with $\theta_1 < 0$. The critical case $\theta_1 = 0$ gives a uniform interior density $f_n(x) = \text{const}$.

---

## 4. Performance Metrics (Exponential Demands)

Throughout, $A = c/[(1-P_B)K]$ normalizes the battery density seen by admitted customers, and $g_B = \pi_0^{(B)}/(1-P_B)$ is the point mass at $x = B$ in the admitted distribution.

### 4.1 Blocking Probability

By **PASTA** (Poisson Arrivals See Time Averages), arriving cars observe the stationary distribution, so $P_B = \Pr(n=1)$.

**Theorem (Blocking Probability).**
$$P_B = \frac{p_0}{K}\!\left[\frac{e^{\theta_1 B}-1}{\theta_1} + \frac{c}{\mu p_1}\right]$$

**Limits:**
$$P_B \xrightarrow{B \to 0} \frac{\lambda}{\lambda + \mu p_1} \quad \text{(M/M/1/1 with slow rate)} \qquad P_B \xrightarrow{B \to \infty} \frac{\lambda}{\lambda + \mu p_2} \quad \text{(M/M/1/1 with fast rate)}$$

> Since $\pi_0^{(0)} = 0$, every admitted car finds $x > 0$ and begins fast charging immediately.

### 4.2 Battery Depletion Probability

**Theorem (Depletion Probability).** *An admitted car with demand $y$ depletes the battery if and only if $y > p_2 x/c$. Averaging over the admitted battery density:*
$$P_D(y) = \begin{cases} \dfrac{A}{\theta_1}\!\left(e^{\theta_1 cy/p_2} - 1\right) & y < p_2 B/c \\[6pt] 1 & y \geq p_2 B/c \end{cases}$$

### 4.3 Charging Time

**Theorem (Charging Time Decomposition).** *Given admission at battery level $x$ with demand $y$:*
$$T(x,y) = \frac{y}{p_2} + \underbrace{\left(\frac{1}{p_1} - \frac{1}{p_2}\right)\max\!\left(y - \frac{p_2 x}{c},\, 0\right)}_{\text{penalty for battery depletion}}$$

The first term is the ideal fast-charging time; the second is the extra time incurred when the battery runs out mid-service.

### 4.4 Mean Charging Efficiency

**Theorem (Mean Charging Efficiency).** *Let $\gamma := \theta_1 - \mu p_2/c = -\lambda/p_0 < 0$ (always). Then:*
$$\frac{\mathbb{E}[T]}{\mathbb{E}[Y]} = \frac{1}{p_2} + \left(\frac{1}{p_1} - \frac{1}{p_2}\right)\!\left[\frac{A(e^{\gamma B}-1)}{\gamma} + g_B\,e^{-\mu p_2 B/c}\right]$$

> **Universal convergence.** Since $\gamma = -\lambda/p_0 < 0$ *always*, regardless of all other parameters:
> $$\frac{\mathbb{E}[T]}{\mathbb{E}[Y]} \;\xrightarrow{B \to \infty}\; \frac{1}{p_2}$$
> Enlarging the battery always improves efficiency, converging to the ideal fast-charge rate.

### 4.5 System-Level Metrics

**Theorem (Mode Fractions, Mean Battery, Mean Grid Power).**

*Stationary time fractions:*
$$\rho_{\mathrm{slow}} = \frac{p_0 c}{\mu p_1 K}, \qquad \rho_{\mathrm{fast}} = \frac{p_0(e^{\theta_1 B}-1)}{\theta_1 K}, \qquad \rho_{\mathrm{idle}} = \frac{c(e^{\theta_1 B}-1)}{\theta_1 K} + \frac{p_0 c}{\lambda K}e^{\theta_1 B}$$

satisfying $\rho_{\mathrm{slow}} + \rho_{\mathrm{fast}} + \rho_{\mathrm{idle}} = 1$ and $\rho_{\mathrm{slow}} + \rho_{\mathrm{fast}} = P_B$.

*Mean battery level:*
$$\mathbb{E}[x] = \frac{1}{K}\!\left[(p_0+c)\!\left(\frac{Be^{\theta_1 B}}{\theta_1} - \frac{e^{\theta_1 B}-1}{\theta_1^2}\right) + \frac{p_0 cB}{\lambda}e^{\theta_1 B}\right]$$

*Mean grid power (bounded in $[p_0, p_1]$):*
$$\mathbb{E}[P_{\mathrm{grid}}] = p_1 - (p_1 - p_0)\,\rho_{\mathrm{idle}}$$

A larger battery increases idle time and thus lowers mean grid draw. As $B \to \infty$, $\rho_{\mathrm{slow}} \to 0$ and slow charging disappears entirely.

### 4.6 Summary Table

| Metric | Formula | Limit $B \to \infty$ |
|--------|---------|---------------------|
| $P_B$ | $\frac{p_0}{K}\!\left[\frac{e^{\theta_1 B}-1}{\theta_1}+\frac{c}{\mu p_1}\right]$ | $\frac{\lambda}{\lambda+\mu p_2}$ |
| $P_D(y)$ | $\frac{A}{\theta_1}(e^{\theta_1 cy/p_2}-1)$ | $1 - e^{\theta_1 cy/p_2}$ |
| $\mathbb{E}[T]/\mathbb{E}[Y]$ | closed form above | $1/p_2$ |
| $\mathbb{E}[x]$ | closed form above | grows $\propto B$ |
| $\mathbb{E}[P_{\mathrm{grid}}]$ | $p_1-(p_1-p_0)\rho_{\mathrm{idle}}$ | $p_1-(p_1-p_0)(1-P_B)$ |
| $\rho_{\mathrm{fast}}$ | $\frac{p_0(e^{\theta_1 B}-1)}{\theta_1 K}$ | $P_B - \rho_{\mathrm{slow}}$ |

All metrics are **explicit** functions of $(p_0, p_1, p_2, \mu, \lambda, B)$ stemming from the single-exponential structure of $f_0$ and $f_1$.

---

## 5. General Energy Demands

### 5.1 Why Exponentiality Was Special

For $Y \sim \mathrm{Exp}(\mu)$, the memoryless property collapses the service history into a scalar and makes $(n, x)$ Markov, leading to an ODE system in $x$ with closed-form solutions.

For a general distribution $F$:
- The residual demand at time $t$ depends on both elapsed time and the charging rate history.
- The departure rate is not a function of $(n, x)$ alone.
- $(n, x)$ is **not** Markov.

**Augmented state.** Let $r(t)$ denote the residual energy demand of the car in service (with $r = 0$ when $n = 0$). Then $(n(t), x(t), r(t))$ is Markov and defines a PDMP on:
$$E = \{0\} \times [0,B] \times \{0\} \;\cup\; \{1\} \times [0,B] \times [0,\infty)$$

### 5.2 Augmented PDMP

**Flow:**
$$\begin{pmatrix}\dot{x}\\\dot{r}\end{pmatrix} = \begin{cases} (p_0,\; 0)^\top & n=0 \\[4pt] (-c,\; -p_2)^\top & n=1,\; x>0 \\[4pt] (0,\; -p_1)^\top & n=1,\; x=0 \end{cases}$$

**Jumps:**
- *Arrivals* (rate $\lambda$, when $n=0$): $(0,x,0) \to (1,x,r')$ with $r' \sim f_Y$.
- *Departures* (when $r$ hits 0): $(1,x,0) \to (0,x,0)$. These are **forced jumps** at the boundary $r = 0$ in the sense of Davis (1984).

Two regions of interest:
- **Region I:** $n=1$, $x>0$, $r>0$ — fast charging, drift $(-c, -p_2)$.
- **Region II:** $n=1$, $x=0$, $r>0$ — slow charging, drift $(0, -p_1)$.

### 5.3 Main Theorem

**Theorem (Stationary Distribution for General $F$).** *Let $f_0(x)$ be the stationary density for $n=0$, $x>0$, and $M_Y(t) = \mathbb{E}[e^{tY}]$ the moment generating function of $Y$.*

**(i) Characteristic solution.** In Region I, the joint density $f_1(x,r)$ for $n=1$, $x>0$, $r>0$ satisfies the transport PDE with source:
$$c\,\partial_x f_1 + p_2\,\partial_r f_1 = \lambda\,f_0(x)\,f_Y(r)$$

The characteristics are lines $p_2 x - cr = \mathrm{const}$. Integrating along them from the entry point of each trajectory (a car arriving at battery level $x + c(R-r)/p_2$ with demand $R$):
$$\boxed{f_1(x,r) = \frac{\lambda}{p_2}\int_r^\infty f_0\!\left(x + \frac{c(R-r)}{p_2}\right)f_Y(R)\,dR}$$

**(ii) Laplace transform of $f_0$.** The departure flux $p_2 f_1(x,0^+)$ feeds back into the $f_0$ balance equation $p_0 f_0' = -\lambda f_0 + p_2 f_1(\cdot, 0^+)$. Taking the Laplace transform in $x$ and using the convolution structure:
$$\boxed{\hat{f}_0(s) = \frac{p_0\,f_0(0^+)}{p_0 s + \lambda\!\left[1 - M_Y\!\left(\dfrac{cs}{p_2}\right)\right]}}$$

where $f_0(0^+)$ is determined by normalization, and
$$\psi(s) := p_0 s + \lambda\!\left[1 - M_Y\!\left(\frac{cs}{p_2}\right)\right]$$
is the **characteristic function** of the problem.

### 5.4 Proof Sketch

**Step 1 — Transport PDE with source.** The balance equation for $f_1(x,r)$ in Region I is
$$c\,\partial_x f_1 + p_2\,\partial_r f_1 = \lambda\,f_0(x)\,f_Y(r)$$
Arrivals inject probability mass at rate $\lambda f_0(x) f_Y(r)$; there are no interior sinks (departures only occur at $r = 0$).

**Step 2 — Integrate along characteristics.** The characteristic through $(x,r)$ was seeded by a car arriving at battery level $X = x + c(R-r)/p_2$ with demand $R > r$. Each trajectory passes through $(x,r)$ exactly once at speed $p_2$ in the $r$-direction, giving the integral formula for $f_1$.

**Step 3 — Close via the $f_0$ equation.** Evaluate the departure flux $p_2 f_1(x,0^+)$ from the characteristic formula, substitute into $p_0 f_0' = -\lambda f_0 + p_2 f_1(\cdot, 0^+)$, and take the Laplace transform. The convolution in $x$ transforms to a product involving $M_Y(cs/p_2)$, giving the algebraic equation for $\hat{f}_0(s)$. $\square$

### 5.5 Recovery of the Exponential Case

For $Y \sim \mathrm{Exp}(\mu)$: $M_Y(t) = \mu/(\mu - t)$, so
$$1 - M_Y\!\left(\frac{cs}{p_2}\right) = \frac{-cs/p_2}{\mu - cs/p_2}$$

Substituting:
$$\hat{f}_0(s) = \frac{p_0\,f_0(0^+)}{s - \theta_1} \;\;\xrightarrow{\mathcal{L}^{-1}}\;\; f_0(x) = p_0\,f_0(0^+)\,e^{\theta_1 x} \qquad (\theta_1 < 0)$$

This recovers the single exponential derived by the ODE method. $\checkmark$

**Structure for general $F$.** The poles of $\hat{f}_0(s)$ — i.e., the zeros of $\psi(s)$ with negative real part — determine the exponential modes of $f_0(x)$:

| Distribution $F$ | $M_Y(cs/p_2)$ | Modes of $f_0$ |
|---|---|---|
| $\mathrm{Exp}(\mu)$ | Möbius, degree 1 | 1 exponential |
| $\mathrm{Erlang}(n,\mu)$ | Rational, degree $n$ | $n$ exponentials |
| Hyperexponential | Rational | Finite sum |
| Phase-type $\mathrm{PH}(\boldsymbol{\alpha}, T)$ | Rational | $|T|$ exponentials |
| Log-normal, Weibull | Non-rational | Infinite sum |

### 5.6 Characteristic Equation and Stability

The characteristic function $\psi(s) = p_0 s + \lambda[1 - M_Y(cs/p_2)]$ satisfies:
- $\psi(0) = 0$ always (trivial zero, corresponds to normalization).
- $\psi'(0) = p_0 - \lambda c/(\mu p_2)$, which is negative if and only if $\theta_1 < 0$.

> **Stability condition is distribution-free.**
> The condition $\theta_1 < 0$, equivalently $\lambda > \lambda_c = \mu p_2 p_0/c$, depends on $F$ **only through its mean** $1/\mu$. The same critical load $\lambda_c$ governs stability for any demand distribution with mean $1/\mu$.

### 5.7 Performance Metrics for General $F$

**Blocking probability** ($B = \infty$). Since $\psi(0) = 0$, L'Hôpital applied to $\hat{f}_0(0) = p_0 f_0(0^+)/\psi(0)$ gives:
$$\hat{f}_0(0) = \frac{p_0\,f_0(0^+)}{\psi'(0)}, \qquad 1 - P_B = \lambda\,\hat{f}_0(0)$$

$$\boxed{P_B = \frac{\lambda c/(\mu p_2 p_0)}{1 + \lambda c/(\mu p_2 p_0) + \mu p_1\,\pi_1^{(0)}/(p_0 f_0(0^+))}}$$

To leading order (large $B$, small $\pi_1^{(0)}$): $P_B \approx (\lambda/\lambda_c)/(1 + \lambda/\lambda_c)$.

**The blocking probability depends on $F$ only through its mean** $1/\mu$ at leading order.

---

**Mean battery level** ($B = \infty$). Using $\mathbb{E}_0[x] = -\hat{f}_0'(0)$ and differentiating $\hat{f}_0(s) = p_0 f_0(0^+)/\psi(s)$, L'Hôpital yields:
$$\mathbb{E}_0[x] = \frac{p_0\,f_0(0^+)\,\psi''(0)}{2\,\psi'(0)^2}$$

Computing $\psi''(0) = -\lambda c^2 \mathbb{E}[Y^2]/p_2^2$:
$$\boxed{\mathbb{E}_0[x] = \frac{\lambda c^2\,\mathbb{E}[Y^2]}{2p_2^2} \cdot \frac{p_0\,f_0(0^+)}{\psi'(0)^2}}$$

**The mean battery level depends on $F$ through its second moment $\mathbb{E}[Y^2]$.** Higher demand variance implies wider fluctuations in the drain cycle and therefore a larger mean battery level.

---

**Mean charging time** (general $F$, any $B$). Using the stop-loss identity $\mathbb{E}[\max(Y-t, 0)] = \int_t^\infty \bar{F}(y)\,dy$:
$$\mathbb{E}[T] = \frac{\mathbb{E}[Y]}{p_2} + \left(\frac{1}{p_1} - \frac{1}{p_2}\right) \int_0^\infty g(x_0) \int_{p_2 x_0/c}^\infty \bar{F}(y)\,dy\; dx_0$$

where $g(x_0)$ is the battery density seen at admission and $\bar{F}(y) = \Pr(Y > y)$.

For $Y \sim \mathrm{Exp}(\mu)$: $\bar{F}(y) = e^{-\mu y}$, the inner integral is explicit, and the formula reduces to Theorem 4.4 above.

As $B \to \infty$ for any $F$ with finite mean: $g(x_0) \to 0$ for small $x_0$ and the depletion integral vanishes, so $\mathbb{E}[T]/\mathbb{E}[Y] \to 1/p_2$ — the universal convergence holds for all $F$.

---

### 5.8 Summary: Exponential vs General Demands

| Feature | Exponential $Y$ | General $Y$ |
|---|---|---|
| State space | $(n, x)$ | $(n, x, r)$ |
| Departure mechanism | Poisson, rate $\delta(x)$ | Deterministic hitting $r=0$ |
| Balance equations | ODE system in $x$ | PDE system in $(x,r)$ |
| Solution method | Eigenvalue decomposition | Characteristics + Laplace |
| $\hat{f}_0(s)$ | $A/(s - \theta_1)$ | $p_0 f_0(0^+)/\psi(s)$ |
| $f_0(x)$ | Single exponential | Sum of exponentials |
| Proportionality $f_1/f_0$ | Constant in $x$ | $x$-dependent |
| Stability condition | $\lambda > \lambda_c$ | Same $\lambda_c$ (mean only) |
| $P_B$ | Explicit closed form | Via $\psi'(0)$; depends on mean only |
| $\mathbb{E}[x]$ | Explicit closed form | Via $\psi''(0)$; needs $\mathbb{E}[Y^2]$ |
| $\mathbb{E}[T]/\mathbb{E}[Y]$ | Explicit closed form | Via stop-loss transform of $F$ |

---

## 6. Open Problems

- **Multi-server extension** ($s > 1$, shared grid): eigenvalue structure $\theta_k^* = k\mu p_2/(p_0 + kp_2 - p_1) - \lambda/p_0$, independent of $s$.
- **Finite $B$ with general $F$**: the slow-charging boundary (Region II) couples back to the fast-charging region through additional boundary conditions not present when $B = \infty$.
- **Non-loss systems**: waiting room ($s$-server queue with finite or infinite waiting space).
- **Phase-type demands**: matrix-analytic approach exploiting the rational structure of $M_Y$.
- **Optimal design**: joint optimization of $B$, $p_0$, $p_2$ given cost and revenue structure.

---

## References

- Davis, M. H. A. (1984). Piecewise-deterministic Markov processes: a general class of non-diffusion stochastic models. *J. Royal Statistical Society B*, 46(3), 353–388.
- Wolff, R. W. (1982). Poisson arrivals see time averages. *Operations Research*, 30(2), 223–231.
- Anick, D., Mitra, D., and Sondhi, M. M. (1982). Stochastic theory of a data-handling system with multiple sources. *Bell System Technical Journal*, 61(8), 1871–1894.
- Asmussen, S. (1995). Stationary distributions for fluid flow models with or without Brownian noise. *Stochastic Models*, 11(1), 21–49.
- Akar, N. and Sohraby, K. (2004). Infinite- and finite-buffer Markov fluid queues: a unified analysis. *Journal of Applied Probability*, 41(2), 557–569.
- Latouche, G. and Ramaswami, V. (1999). *Introduction to Matrix Analytic Methods in Stochastic Modeling*. SIAM.
