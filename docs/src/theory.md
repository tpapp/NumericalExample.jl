## Theory

The objective of the exercise is solving a simple Ramsey model with taxes in the time domain. The setup and characterization is of course standard, we write it out here in detail to make the presentation self-contained.

### Model setup

The representative agent maximizes lifetime utility
```math
\sum_{t=0}^\infty \beta^t \bigl[u(c_t, \ell_t) + v(g_t)\bigr]
```
where $c_t$ and $\ell_t$ are consumption and labor in the time period $t$, and $g_t$ is government expenditure, financed from taxes.

The two production factors labor ``\ell_t`` and capital ``k_t``, and the production function is
```math
F(k, \ell) = k^\alpha \ell^{1-\alpha}
```
where we have normalized TFP to ``1`` without loss of generality. Capital is held by the agents, depreciates at rate ``\delta``,  and follows the law of motion
```math
k_{t+1} = (1-\tau_k) r_t k_t + (1-\delta) k_t + (1-\tau_\ell) w_t \ell_t - c_t
```
where ``\tau_k`` and ``\tau_\ell`` are capital and labor tax rates, which finance government expenditure
```math
g_t = \tau_k r_t k_t + \tau_\ell w_t \ell_t
```
The functional form of ``v`` is irrelevant, as utility is separable. For ``u``, we use
```math
u(c, \ell) = \theta \log(c) + (1-\theta) \log(1-\ell)
```
ie the total labor endownment is normalized to ``1``.

We solve for the *competitive equilibrium*, ie a sequence of ``\{ c_t, g_t, \ell_t, k_t, r_t, w_t\}`` for which the choices of the agents are optimal given the prices, the resource constraints hold, and markets clear.

We assume the usual conditions to get a well-behaved solution, including sufficient strict concavity where it matters.

### Equilibrium characterization

Firms maximize profits with a homothetic production function, so prices equal marginal product, specifically
```math
r_t = \alpha k_t^{\alpha -1 }\ell_t^{1-\alpha}
```
and
```math
w_t = (1-\alpha) k_t^{\alpha }\ell_t^{-\alpha}
```

Setting up the Lagrangean as

```math
\mathcal{L} = \sum_{t=0}^\infty \beta^t \Bigl[ u(c_t, \ell_t) +
\mu_t \bigl( (1 - \tau_k) r_t k_t + (1-\delta) k_t + (1-\tau_\ell) w_t \ell_t - c_t - k_{t+1} \bigr)
\Bigr]
```

and obtain the FOC

```math
\frac{\partial \mathcal{L}}{\partial c_t} = \beta^t \left( \frac{\partial}{\partial c} u(c_t, \ell_t) - \mu_t \right) = 0 \\
\frac{\partial \mathcal{L}}{\partial \ell_t} = \beta^t \left( \frac{\partial}{\partial \ell} u(c_t, \ell_t) - \mu_t (1 -\tau_\ell) w_t \right) = 0\\
\frac{\partial \mathcal{L}}{\partial k_{t+1}} = \beta^{t+1} \mu_{t+1} \bigl(  (1-\tau_k)r_{t+1} + 1- \delta \bigr) - \beta^t \mu_t  = 0
```

which we combine to the usual intratemporal

```math
\frac{\partial}{\partial \ell} u(c_t, \ell_t) + \frac{\partial}{\partial c} u(c_t, \ell_t) (1-\tau_\ell) w_t = 0
```

and intertemporal FOC

```math
 \frac{\partial}{\partial c} u(c_t, \ell_t) = \beta \frac{\partial}{\partial c} u(c_{t+1}, \ell_{t+1}) \bigl( (1-\tau_k) r_{t+1} + 1 - \delta \bigr)
```

### Steady state

It is important to solve for the *steady state*: first, this is what the model converges to, second, it provides a good initial guess. For this model, we can do this analytically, but for less tractable models you may need a numerical rootfinder.

Let ``\bar{c}, \bar{k}, \dots`` denote the steady state variables. Also introduce ``\lambda = \bar{\ell} / \bar{k}`` for the steady state labor/capital ratio.  Then

```math
\bar{r} = \alpha \bar{k}^{\alpha - 1} \bar{\ell}^{1 - \alpha} = \alpha \lambda^{1-\alpha}
```

and

```math
\bar{w} = (1-\alpha) \bar{k}^\alpha \bar{\ell}^{-\alpha} = (1-\alpha) \lambda^{-\alpha}
```

Then the intertemporal condition in the steady state can be written as

```math
1 = \beta \bigl( (1-\tau_k) \alpha \lambda^{1-\alpha} + 1 - \delta \bigr)
```

which can easily be solved for ``\lambda``. The budget constraint of the agent can be divided by capital in the steady state to yield

```math
1 = (1-\tau_k) \alpha \lambda^{1-\alpha} + 1 - \delta + (1-\tau_\ell) (1-\alpha) \lambda^{1-\alpha} - \frac{\bar{c}}{\bar{k}}
```

which yields ``\bar{c}/\bar{k}``. Finally, we can write the intratemporal FOC as
```math
\frac{1-\theta}{1/\bar{k}-\lambda} = \frac{\theta}{\bar{c}/\bar{k}} ( 1-\tau_\ell) (1-\alpha) \lambda^{-\alpha}
```

which we can solve for ``\bar{k}``, and then obtain ``\lambda`` and ``\bar{c}``.
