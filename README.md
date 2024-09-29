SymDR
========
SymDR is a library made to automate the process of finding dispersion relations.

It was created in The Great Mathematical Workshop 2024 by the command of young ambitious scientists.

Features
===========
Finding dispersion relation in:
- Equations
  ```math
  u_t + u u_x = A^2 u_xx
  ```
- System of equations
```math
\begin{equation*}
\left\{
 \begin{array}{lcl}
%\begin{array}{l}
\displaystyle v_t+vv_x+ \frac{v-u}{\tau} = 0,\\
\displaystyle u_t + uu_x -\nu u_{xx} + \frac{\varepsilon(u-v)}{\tau} = 0.
\end{array}
\right.
\end{equation*}
```

- Discrete equations
  ```math
  [u_t+u]_a+\frac{u_{a+1}-u_{a-1}}{2h}=0
  ```
- Systems of discrete equations

```math
\begin{equation*}
\left\{
 \begin{array}{lcl}
\displaystyle v_t+v \dfrac{v(x+h)-v(x-h)}{2h} + \frac{v-u}{\tau} = 0,\\
\displaystyle u_t + u \dfrac{u(x+h)-u(x-h)}{2h} -\nu u_{xx} + \frac{\varepsilon(u-v)}{\tau} = 0.
\end{array}
\right.
\end{equation*}
```


Quickstart
=============
For the best experience, we highly recommend using [Jupyter Notebook](https://jupyter.org/) or [Google Collab](https://colab.research.google.com).

SymDR supports [Python](http://python.org/) >=3.7 and only depends on [Sympy](https://www.sympy.org).
Install the Python package using pip.
```
pip install symdr
```

Example
==========
- If you are new to symbolic mathematics in Python, read SymPy [introduction](https://docs.sympy.org/latest/tutorials/intro-tutorial/index.html) first.
- For the end-to-end example of equation analysis, see [notebook](https://github.com/symdr/symdr/blob/master/example.ipynb).

>!!!!!!! IMPORTANT !!!!!!!!
>Variables `x`, `t`, `w`, `k` must not be redefined. They are used by the algorithm. When working with discrete cases, the variables `h`, `tau`, `a` and `n` are also added to this list.


SymDR has own objects for grid functions. To turn on pretty-printing for these functions, you need to add:

```
init_printing(latex_printer=discrete_latex_printer)
```
Creation of a grid function is the same as standart function:

```
u = DiscreteGrid('u')
```

Consider a discrete analog of Korteweg–De Vries equation:
```math
u_{t}=6uu_{x}-u_{xxx}
```

For the third order derivative, we are going to use next scheme:
```math
f^{(3)}(x)=\frac{f(x+2h)-2f(x+h)+2f(x-h)-f(x-2h)}{2h^3} 
```
And we can write our equation
```
equation = u.diff(t) + (u.at_x(a+2) - 2 * u.at_x(a+1) + 2 * u.at_x(a-1) - u.at_x(a-2)) / (2 * h ** 3)
```
>The class allows us to denote a derivative at (a, n) with "diff" method exactly the same way as SymPy functions do.
>Shift at a point is done by any of three methods:
>- at_x - shift in space (i.e. ```u.at_x(a+2)``` means $u_{a+2}^n$)
>- at_t - shift in time (i.e. ```u.at_t(n-2)``` means $u_a^{n-2}$)
>- at - shift in both axes at once (i.e. ```u.at(a+2, n-1)``` means $u_{a+2}^{n-1}$)

>Technical note: instead of moving the function into a subtree of ```Derivative``` object, as it is with SymPy functions, differentiation of ```DiscreteGrid``` object simply returns an object of the same class, but with different arguments.

Finally, let's find dispersion relation:

```
d_equation_dr(equation)
```
Also you can rewrite it with Euler formula:

```
d_equation_dr(equation, trig_rewrite=True)
```

And here's the full code:

```
from sympy import *
from symdr import *

u = DiscreteGrid('u')
equation = u.diff(t) + (u.at_x(a+2) - 2 * u.at_x(a+1) + 2 * u.at_x(a-1) - u.at_x(a-2)) / (2 * h ** 3)

d_equation_dr(equation)
```
