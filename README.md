SymDR
========
SymDR - is the library made to automate the process of finding dispersion relations.

It was created in The Great Mathematical Workshop 2024 by the command of young ambitious scientists.

Features
===========
Finding dispersion relation in:
- Equations
$u_t + u u_x = A^2 u_xx$
- System of equations

- Discrete equation
$[u_t+u]_a+ \frac{u_{a+1}-u_{a-1}}{2h}=0$
- Systems of discrete equations


Quickstart
=============
For the best experience we are highly recommend to use [Jupyter Notebook](https://jupyter.org/) or [Google Collab](https://colab.research.google.com).

SymDR supports [Python](http://python.org/) >=3.7 and only depend on [Sympy](https://www.sympy.org).
Install the python package using pip
```
pip install symdr
```

Example
==========
- If you are new in symbolic mathematics in python read SymPy [introduction](https://docs.sympy.org/latest/tutorials/intro-tutorial/index.html) first.
- For the end to end example of equation analysis see [notebook](https://github.com/symdr/symdr/blob/master/example.ipynb)

>!!!!!!! IMPORTANT !!!!!!!!
>Variables `x`, `t`, `w`, `k` cannot be redefined. They are used for the algorithm. When working with discrete cases, the variables `h`, `tau`, `a` and `n` are also added to this list


SymDR have own objects for grid functions. To pretty printing this functions you need to write:

```
init_printing(latex_printer=discrete_latex_printer)
```
Creating the grid function is the same to standart function:

```
u = DiscreteGrid('u')
```

Let see discrete analog of Korteweg–De Vries equation:
$$ u_{t}=6uu_{x}-u_{xxx} $$

For the third order derivative we going to use next scheme:
$$ f^{(3)}(x)=\frac{f(x+2h)-2f(x+h)+2f(x-h)-f(x-2h)}{2h^3}$$

And we can write our equation
```
equation = u.diff(t) + (u.at_x(a+2) - 2 * u.at_x(a+1) + 2 * u.at_x(a-1) - u.at_x(a-2)) / (2 * h ** 3)
```
>The class allows us to denote a derivative at (a, n) with "diff" method exactly the same way as SymPy functions do.
>Shift at point is done by any of three methods:
>- at_x - shift in space (i.e. ```u.at_x(a+2)``` means $u_{a+2}^n$)
>- at_t - shift in time (i.e. ```u.at_t(n-2)``` means $u_a^{n-2}$)
>- at - shift in both axes at once (i.e. ```u.at(a+2, n-1)``` means $u_{a+2}^{n-1}$)

>Technical note: instead of moving function into a subtree of "Derivative" object, as it is with SymPy functions, differentiation of "DiscreteGrid" object simply returns an object of the same class, but with different arguments.

Finally let's find disperssion relation:

```
d_equation_dr(equation)
```
Also you can rewrite it with Euler formula :

```
d_equation_dr(equation, trig_rewrite=True)
```

And full code:

```
from sympy import *
from symdr import *

u = DiscreteGrid('u')
equation = u.diff(t) + (u.at_x(a+2) - 2 * u.at_x(a+1) + 2 * u.at_x(a-1) - u.at_x(a-2)) / (2 * h ** 3)

d_equation_dr(equation)
```
