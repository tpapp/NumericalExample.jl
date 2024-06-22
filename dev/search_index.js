var documenterSearchIndex = {"docs":
[{"location":"theory/#Theory","page":"Theory","title":"Theory","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"The objective of the exercise is solving a simple Ramsey model with taxes in the time domain. The setup and characterization is of course standard, we write it out here in detail to make the presentation self-contained.","category":"page"},{"location":"theory/#Model-setup","page":"Theory","title":"Model setup","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"The representative agent maximizes lifetime utility","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"sum_t=0^infty beta^t biglu(c_t ell_t) + v(g_t)bigr","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where c_t and ell_t are consumption and labor in the time period t, and g_t is government expenditure, financed from taxes.","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"The two production factors labor ell_t and capital k_t, and the production function is","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"F(k ell) = k^alpha ell^1-alpha","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where we have normalized TFP to 1 without loss of generality. Capital is held by the agents, depreciates at rate delta,  and follows the law of motion","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"k_t+1 = (1-tau_k) r_t k_t + (1-delta) k_t + (1-tau_ell) w_t ell_t - c_t","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where tau_k and tau_ell are capital and labor tax rates, which finance government expenditure","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"g_t = tau_k r_t k_t + tau_ell w_t ell_t","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"The functional form of v is irrelevant, as utility is separable. For u, we use","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"u(c ell) = theta log(c) + (1-theta) log(1-ell)","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"ie the total labor endownment is normalized to 1.","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"We solve for the competitive equilibrium, ie a sequence of  c_t g_t ell_t k_t r_t w_t for which the choices of the agents are optimal given the prices, the resource constraints hold, and markets clear.","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"We assume the usual conditions to get a well-behaved solution, including sufficient strict concavity where it matters.","category":"page"},{"location":"theory/#Equilibrium-characterization","page":"Theory","title":"Equilibrium characterization","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"Firms maximize profits with a homothetic production function, so prices equal marginal product, specifically","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"r_t = alpha k_t^alpha -1 ell_t^1-alpha","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"and","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"w_t = (1-alpha) k_t^alpha ell_t^-alpha","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Setting up the Lagrangean as","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"mathcalL = sum_t=0^infty beta^t Bigl u(c_t ell_t) +\nmu_t bigl( (1 - tau_k) r_t k_t + (1-delta) k_t + (1-tau_ell) w_t ell_t - c_t - k_t+1 bigr)\nBigr","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"and obtain the FOC","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"fracpartial mathcalLpartial c_t = beta^t left( fracpartialpartial c u(c_t ell_t) - mu_t right) = 0 \nfracpartial mathcalLpartial ell_t = beta^t left( fracpartialpartial ell u(c_t ell_t) - mu_t (1 -tau_ell) w_t right) = 0\nfracpartial mathcalLpartial k_t+1 = beta^t+1 mu_t+1 bigl(  (1-tau_k)r_t+1 + 1- delta bigr) - beta^t mu_t  = 0","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"which we combine to the usual intratemporal","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"fracpartialpartial ell u(c_t ell_t) + fracpartialpartial c u(c_t ell_t) (1-tau_ell) w_t = 0","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"and intertemporal FOC","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":" fracpartialpartial c u(c_t ell_t) = beta fracpartialpartial c u(c_t+1 ell_t+1) bigl( (1-tau_k) r_t+1 + 1 - delta bigr)","category":"page"},{"location":"theory/#Steady-state","page":"Theory","title":"Steady state","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"It is important to solve for the steady state: first, this is what the model converges to, second, it provides a good initial guess. For this model, we can do this analytically, but for less tractable models you may need a numerical rootfinder.","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Let barc bark dots denote the steady state variables. Also introduce lambda = barell  bark for the steady state labor/capital ratio.  Then","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"barr = alpha bark^alpha - 1 barell^1 - alpha = alpha lambda^1-alpha","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"and","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"barw = (1-alpha) bark^alpha barell^-alpha = (1-alpha) lambda^-alpha","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Then the intertemporal condition in the steady state can be written as","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"1 = beta bigl( (1-tau_k) alpha lambda^1-alpha + 1 - delta bigr)","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"which can easily be solved for lambda. The budget constraint of the agent can be divided by capital in the steady state to yield","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"1 = (1-tau_k) alpha lambda^1-alpha + 1 - delta + (1-tau_ell) (1-alpha) lambda^1-alpha - fracbarcbark","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"which yields barcbark. Finally, we can write the intratemporal FOC as","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"frac1-theta1bark-lambda = fracthetabarcbark ( 1-tau_ell) (1-alpha) lambda^-alpha","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"which we can solve for bark, and then obtain lambda and barc.","category":"page"},{"location":"numerical/#Numerical-approximation","page":"Numerical approximation","title":"Numerical approximation","text":"","category":"section"},{"location":"numerical/","page":"Numerical approximation","title":"Numerical approximation","text":"We set up a container for the model parameters, but it is not exported from the module.","category":"page"},{"location":"numerical/","page":"Numerical approximation","title":"Numerical approximation","text":"NumericalExample.ModelParameters","category":"page"},{"location":"numerical/#NumericalExample.ModelParameters","page":"Numerical approximation","title":"NumericalExample.ModelParameters","text":"A container for the model parameters. Use model_parameters to create, which performs some simple checks.\n\nθ: weight on (log) consumption in period utility\nβ: discount factor\nα: capital share\nδ: depreciation\nτ_k: tax rate on capital returns\nτ_ℓ: tax rate on wages\n\n\n\n\n\n","category":"type"},{"location":"numerical/","page":"Numerical approximation","title":"Numerical approximation","text":"We use a function to instantiate it, which is better style in Julia and allows for seamless extension of the interface later.","category":"page"},{"location":"numerical/","page":"Numerical approximation","title":"Numerical approximation","text":"model_parameters","category":"page"},{"location":"numerical/#NumericalExample.model_parameters","page":"Numerical approximation","title":"NumericalExample.model_parameters","text":"model_parameters(; θ, β, α, δ, τ_k, τ_ℓ)\n\n\nCreate a value that holds the parameters of the model.\n\nKeyword arguments\n\nθ: the parameter θ log(c) + (1-θ) log(1-ℓ) in the utility function\nβ: the discount factor\nα: the capital share, in the production function k^α ℓ^1-α\nδ: the depreciation rate\nτ_k, τ_ℓ: tax rates for capital and labor, respectively\n\n\n\n\n\n","category":"function"},{"location":"numerical/","page":"Numerical approximation","title":"Numerical approximation","text":"We code the elements of the model in simple functions, which are easier to check and debug.","category":"page"},{"location":"numerical/","page":"Numerical approximation","title":"Numerical approximation","text":"NumericalExample.rent_and_wage\nNumericalExample.tax_revenue\nNumericalExample.period_budget","category":"page"},{"location":"numerical/#NumericalExample.rent_and_wage","page":"Numerical approximation","title":"NumericalExample.rent_and_wage","text":"rent_and_wage(model, k, ℓ)\n\n\nCalculate the rent and the wage in the model, using capital k and labor ℓ.\n\n\n\n\n\n","category":"function"},{"location":"numerical/#NumericalExample.tax_revenue","page":"Numerical approximation","title":"NumericalExample.tax_revenue","text":"tax_revenue(model, k, ℓ)\n\n\nCalculate the tax revenue in the model, using capital k and labor ℓ.\n\n\n\n\n\n","category":"function"},{"location":"numerical/#NumericalExample.period_budget","page":"Numerical approximation","title":"NumericalExample.period_budget","text":"period_budget(model, k, ℓ)\n\n\nCalculate the budget of the agent each period (can be spent on consumption and next-period capital).\n\n\n\n\n\n","category":"function"},{"location":"numerical/","page":"Numerical approximation","title":"Numerical approximation","text":"NumericalExample.euler_residual\nNumericalExample.labor_FOC_residual","category":"page"},{"location":"numerical/#NumericalExample.euler_residual","page":"Numerical approximation","title":"NumericalExample.euler_residual","text":"euler_residual(model; c, k′, ℓ′, c′)\n\n\nCalculate the residual of the intertemporal equilibrium condition (Euler equation) at a given point in time.\n\nc is consumption for this period, k′, ℓ′, c′ are capital, labor, and consumption next period, respectively.\n\n\n\n\n\n","category":"function"},{"location":"numerical/#NumericalExample.labor_FOC_residual","page":"Numerical approximation","title":"NumericalExample.labor_FOC_residual","text":"labor_FOC_residual(model; k, ℓ, c)\n\n\nCalculate the residual of the intratemporal first order condition.\n\nk, ℓ, c are capital, labor, and consumption this period, respectively.\n\n\n\n\n\n","category":"function"},{"location":"application/#Simple-comparative-statics","page":"Simple comparative statics","title":"Simple comparative statics","text":"","category":"section"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"Again, this model is much to simplified to give us any realistic insights into economic policy, so we just demonstrate the comparison of transients.","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"First, we set up the work environment. Accessors.jl is a neat library for modifying elements of a container, here our model.","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"using NumericalExample, Accessors, Miter","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"We pick some reasonable model parameters (again, this is not calibrated to anything realistic, but parameters broadly match usual RBC targets), and solve the model.","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"model = model_parameters(; α = 0.33, β = 0.95, θ = 0.5, δ = 0.05, τ_k = 0.1, τ_ℓ = 0.4)\nk0 = 1.2\nsol = solve_model(model, k0)","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"Then experiment with raising capital taxes from 10% to 30%.","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"model2 = @set model.τ_k = 0.3\nsol2 = solve_model(model2, k0)","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"For plotting, a simple utility function we wrote up in this package. It uses Miter.jl, but any Julia plotting solution would work fine. We use a range of evenly spaced values 20 years, and label our models:","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"ts = range(0, 20; length = 20)\ngraph_labels = [\"baseline\", math\"\\mathrm{capital\\ tax}\\uparrow\"]","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"Let's plot the evolution of capital first. ","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"plot_vs_time(t -> [sol.k̃(t), sol2.k̃(t)], ts, \"capital\"; graph_labels);\nMiter.save(\"capital_vs_time.svg\", ans); nothing # hide","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"With higher taxes, the steady state is lower, so capital converges to a lower value.","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"(Image: )","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"If we look at labor instead, it is apparent that labor is subtituted for capital in production.","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"plot_vs_time(t -> [sol.ℓ̃(t), sol2.ℓ̃(t)], ts, \"labor\"; graph_labels);\nMiter.save(\"labor_vs_time.svg\", ans); nothing # hide","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"(Image: )","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"However, consumption is lower when capital taxes are higher, since the economy produces less.","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"plot_vs_time(t -> [sol.c̃(t), sol2.c̃(t)], ts, \"consumption\"; graph_labels);\nMiter.save(\"consumption_vs_time.svg\", ans); nothing # hide","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"(Image: )","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"However, as expected, tax revenue increases. Since we did not define the utility for the government expenditure, we do not compare welfare.","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"plot_vs_time(t -> [tax_revenue(model, sol.k̃(t), sol.ℓ̃(t)),\n                   tax_revenue(model2, sol2.k̃(t), sol2.ℓ̃(t))],\n             ts, \"tax revenue\"; graph_labels);\nMiter.save(\"tax_revenue_vs_time.svg\", ans); nothing # hide","category":"page"},{"location":"application/","page":"Simple comparative statics","title":"Simple comparative statics","text":"(Image: )","category":"page"},{"location":"#Overview","page":"Overview","title":"Overview","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"This is a self-contained worked example of solving a very simple macroeconomic model in Julia and doing some comparative analysis with the results.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"Since this is just a demonstration, the model is of course not realistic, more like something you would see in a homework problem in an introductory graduate course. The model would need to be extended with a more detailed production structure, capital and labor frictions, nontrivial consumers, etc, and estimated based on data, to be considered useful for actual policy analysis.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"The code is kept simple and it is well-documented. An important part of this demonstration is how all components are tested in detail to ensure correctness.","category":"page"}]
}
