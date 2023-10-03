# ConductionRadiationCoupling
Python code that solves the hyperbolic thermal transport equation in a 1-D slab coupled to the near-field radiative heat transfer from a neighboring slab. In this work we describe the thermalization dynamics and temperature distributions in bodies
composed of hyperbolic materials described by the Cattane-Vernotte equation,  exchanging energy via near-field thermal radiation. As opposed
to Fourier materials Non-Fourier hyperbolic materials do not assume an instant response between
the temperature gradiend and heat flux. The Cattaneo-Vernotte heat conduction equation can be written as

$$-\kappa \nabla  T = \vec{q} + \tau \frac{\partial \vec{q}}{\partial t} $$

where $\kappa$ is material thermal conductivity, $\vec{q}$ the heat flux, $ \tau $ the material response time and where the second term on the right hand side is a first order approximation to the time lag in heat conduction which  depends on the
given material. This time delay can alter the temperature distribution
inside the bodies and slow down the thermalization dynamic. Our results can have important
repercusions and applications in hyperthermia photothermal cancer treatment, since it has been
shown that hyperbolic equations better describe the energy transfer in this type of medium. We solve the 1-D case in which the energy balance
equation can be written as
$$ \frac{ \partial u}{\partial t}  = \Phi_c +  \Phi_r$$

where $u$ is the energy density in the body, $\Phi_c$ and $\Phi_r$ are the heat flux contribution to conduction and radiation respectively.
## Mode Analysis
The solution to the temperature dynamics  in each body can be founbd as a sum of thermal modes 
$$ T_1(z,t) = T_{20} + \sum_{n=0}^{\infty}\left( a_n e^{r^+_{1n}t} + b_n e^{r^-_{1n}t} \right) \cos \left( x_n \frac{(z+L)}{L}\right) $$ 

The python code carries out an analysis of the mode lifetimes for a given set of geometrical and dielectrical parameters
## Reference
For further discussion see _Transient effects in the coupling of thermal radiation and non-Fourier heat transport at the nano-scale_
