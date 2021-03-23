Diff Test
=========

This is :download:`a good method<https://www.cs.ucr.edu/~craigs/papers/2016-derivatives/paper.pdf>` summarized by
Professor Craig Schroeder.

We follow Craig's note to test our derivatives. 
To check whether the implementation of :math:`f'(x)` is correct, we take both :math:`f'` 
and :math:`f'`, evaluate them at the following points:

.. math::

    \frac{f(x+\Delta x)-f(x-\Delta x)}{2 \Delta x} - \frac{f'(x+\Delta x)+f'(x-\Delta x)}{2} = O(\Delta x^2 )

For multivariable function :math:`g(\mathbf{x})`, we have:

.. math::

    \frac{g(\mathbf{x} + \delta \mathbf{x}) - g(\mathbf{x} - \delta \mathbf{x}) - (g'(\mathbf{x} + \delta \mathbf{x}) + g'(\mathbf{x} - \delta \mathbf{x}))\cdot \delta \mathbf{x}}{\delta} = O(\delta^2) 

Here :math:`\delta \mathbf{x}` is is a small incremental value, 
whose entries are chosen uniformly from :math:`[-\delta, \delta]`.
The above formula should give a :math:`O(\delta ^2 )` error for correct :math:`g'`, while having a :math:`O(1)` error for incorrect implementation. 
In practice, we choose :math:`\delta = \epsilon ^{1/3}`, where :math:`\epsilon` is the machine precision.

Furthermore, If we replace scalar function :math:`g` by a vector function :math:`\mathbf{g}`, its derivative :math:`\mathbf{g}'` becomes a matrix.

.. math::

    \frac{\mathbf{g}(\mathbf{x} + \delta \mathbf{x}) - \mathbf{g}(\mathbf{x} - \delta \mathbf{x}) - (\mathbf{g}'(\mathbf{x} + \delta \mathbf{x}) + \mathbf{g}'(\mathbf{x} - \delta \mathbf{x}))\cdot \delta \mathbf{x}}{\delta} = O(\delta^2) 

There are several ways to check whether the above condition is satisfied. 
One way to do so is to randomly select :math:`\delta \mathbf{x}` within :math:`[-\delta, \delta]`.
Compute the error with the correct gradient :math:`g'` and a fake gradient, e.g. :math:`2g'`. Check whether the ratio between fake error and real error is sufficiently small. Use the following code to check the implementation of gradient:

.. code-block:: C++

    const auto f = [&](const Eigen::VectorXd& s) -> double {
        // implementation of f
        // return f(s)
    };
    const auto g = [&](const Eigen::VectorXd& s, Eigen::VectorXd& grad) {
        // implementation of f'
        // set grad as f'(x)
    };
    Eigen::VectorXd s(N);
    CHECK(Bow::FiniteDiff::check_gradient(s, f, g, 1e-4, 1e-3));

Similarly, we can check the implementation of hessian:

.. code-block:: C++

    const auto g = [&](const Eigen::VectorXd& X, Eigen::VectorXd& grad) {
        // implementation of g
        // set grad as g(x)
    };
    const auto h = [&](const Eigen::VectorXd& X, Eigen::SparseMatrix<double>& hess) {
        // implementation of g'
        // set hess as g'(x)
    };
    Eigen::VectorXd X(N);
    CHECK(Bow::FiniteDiff::check_jacobian(X, g, h, 1e-6, 1e-3));

Another way to check the implementation is to check the convergence rate. shrink :math:`\delta` by half and checks whether the log error is decreasing by a certain step.
Here's an example for checking the differentiation:

.. code-block:: C++

    const auto f = [&](const Eigen::VectorXd& x_vec) -> double {
        Bow::Field<Bow::Vector<T, dim>> x = Bow::to_field<dim>(x_vec);
        this->m_energy_terms[0]->precompute(x);
        return this->m_energy_terms[0]->energy(x);
    };
    const auto g = [&](const Eigen::VectorXd& x_vec, Eigen::VectorXd& grad_vec) {
        Bow::Field<Bow::Vector<T, dim>> x = Bow::to_field<dim>(x_vec);
        Bow::Field<Bow::Vector<T, dim>> grad;
        this->m_energy_terms[0]->precompute(x);
        this->m_energy_terms[0]->gradient(x, grad);
        grad_vec = Bow::to_vec(grad);
    };
    const auto h = [&](const Eigen::VectorXd& x_vec, Eigen::SparseMatrix<T>& hess) {
        Bow::Field<Bow::Vector<T, dim>> x = Bow::to_field<dim>(x_vec);
        this->m_energy_terms[0]->precompute(x);
        this->m_energy_terms[0]->hessian(x, hess, false);
    };
    const auto project = [&](Eigen::VectorXd& step_vec) {
        Bow::Field<Bow::Vector<T, dim>> step = Bow::to_field<dim>(step_vec);
        grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
            TV x_n = node.template cast<T>() * dx;
            for (int d = 0; d < BC_order[g.idx]; ++d) {
                TV n = BC_basis[g.idx].col(d);
                step[g.idx] -= step[g.idx].dot(n) * n;
            }
        });
        step_vec = to_vec(step);
    };
    const auto transform = [&](Eigen::VectorXd& step_vec) {
        Bow::Field<Bow::Vector<T, dim>> step = Bow::to_field<dim>(step_vec);
        grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
            step[g.idx] = BC_basis[g.idx].transpose() * step[g.idx];
        });
        step_vec = to_vec(step);
    };
    FiniteDiff::ziran_check_false(x, f, g, h, project, transform);
