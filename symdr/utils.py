from sympy import Function, symbols, exp, I, Integer, Symbol
from .discrete_funcs import DiscreteGrid, DiscreteGridBase

x, t = symbols("x t")  # Equation parameters
k, w = symbols("k omega")  # Wave parameter
h, tau = symbols("h tau")

def _get_function_list(system):
    funcs = set()
    for equation in system:
        funcs = funcs.union(equation.atoms(Function), {i.to_grid() for i in equation.atoms(DiscreteGridBase)})

    return list(funcs)

def _get_deriv(term):
    for i in term:
        if i._is_diff or i._is_shifted:
            return [i]
    return []

def is_continuous(expr):
    if isinstance(expr, DiscreteGrid):
        return expr.is_continious
    elif isinstance(expr, (Integer, Symbol)):
        return True
    else:
        return all(is_continuous(arg) for arg in expr.args)
