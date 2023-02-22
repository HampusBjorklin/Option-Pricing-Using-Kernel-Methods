# Option-Pricing-Using-Kernel-Methods
## abstract
This study takes on the task of pricing a high dimensional basket option using deterministic
and efficient kernel method. The original method is proposed in a pre-print by Christian Rieger
and Holger Wendland. The method is used to approximate high dimensional functions by their
lower dimensional parts in a Reproducing Kernel Hilbert Space. The high dimensional function
is solved on a set of anchor points using freezing projections. To approximate the function, the
result from this set of anchor points is subsequently interpolated onto the rest of the domain. It
can be shown that under certain assumptions on the smoothness on a closed Hilbert subspace,
the curse of dimensionality can be circumvented with the method. This method is used for
approximating the solution to the Black-Scholes PDE using a multiquadric radial basis function
as kernel inside the reproducing kernel proposed by Rieger and Wendland. The method is
tailored to the Black-Scholes PDE by deriving the first and second order derivatives to the
reproducing kernel as well as performing a coordinate system shift in order to capture the
directions of highest importance in the solution domain. The implementation result in
reasonable errors at polynomial complexity instead of the exponential complexity faced by
conventional methods.
