Boundary Condition
==================

:download:`Download Minchen's derivation <../assets/Boundary_Conditions.pdf>`

Change of variable
^^^^^^^^^^^^^^^^^^^^^^^^

A boundary condition defined on node :math:`x_i` along directon :math:`n_i^1` can be expressed as :math:`(n_i^1)^T x_i = \hat{z}_{i}^1`.

This node can be parameterized as :math:`x_i = B_i z_i`, where :math:`B_i = [n_i^1, n_i^2, n_i^3]` is an othogonal matrix extended from :math:`n_i^1`. Now the constraint becomes :math:`z_i^1 = \hat{z}_{i}^1`.

Note that if the node :math:`x_i` is prescribed on two orthogonal directons :math:`n_i^1` and :math:`n_i^2`, we only need to extend :math:`n_i^3`, and the constraint will be on :math:`z_i^1` and :math:`z_i^2`.

The position vector :math:`\mathbf{x} = [x_1; x_2;...;x_N]` can be parameterized as :math:`\mathbf{x} = B\mathbf{z}`, where :math:`B=\begin{bmatrix}B_1 & & & \\ & B_2 & &\\ && ...&\\&&&B_N\end{bmatrix}` is a diagonal block orthogonal matrix.

The original function to be optimized is now :math:`\hat{e}(\mathbf{z}) = e(B\mathbf{z})`. The derivatives on :math:`\mathbf{z}` are :math:`\nabla_\mathbf{z} \hat{e} = B^T \nabla_\mathbf{x} e` and :math:`\nabla^2_\mathbf{z} \hat{e}(\mathbf{z}) = B^T (\nabla^2_\mathbf{x} e) B`. The constraints are just to prescribe some bits of :math:`z`.


Augmented Lagrangian Method
^^^^^^^^^^^^^^^^^^^^^^^^
When the boundary conditions are all satisfied, we can use subspace technique to reduce system's DOFs. Otherwise, we need augmented Lagrangian method to pull the constrained node to its target along some direction.

In :math:`z` domain, it's trivial to apply the augmented Lagrangian method. The function to be optimized is:

.. math::

    L(\mathbf{z}) = \hat{e}(\mathbf{z}) - \sum_{i\in \mathcal{I}}\mathbf{\lambda}^i (\mathbf{z}^i - \hat{\mathbf{z}}^i)+ \frac{\kappa}{2} \sum_{i\in \mathcal{I}}(\mathbf{z}^i - \hat{\mathbf{z}}^i)^2



