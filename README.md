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

$$  \frac{ \partial u}{\partial t}  = \Phi_c +  \Phi_r$$

where $u$ is the energy density in the body, $\Phi_c$ and $\Phi_r$ are the heat flux contribution to conduction and radiation respectively. Further references can be found on the soon to be published work :Thermalization of Non-Fourier Hyperbolic Bodies in the
Near-field. 
