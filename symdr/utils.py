from sympy import Function, symbols, exp, I, Integer, Symbol
from .discrete_funcs import DiscreteGrid, DiscreteGridBase

x, t = symbols("x t")  # Equation parameters

k, w, U, C = symbols("k omega U C")  # Wave parameter

h, tau = symbols("h tau")

wave = U * exp(I * (k * x - w * t))  # The wave


def _get_C_from_deriv(deriv):
    
    diffs = dict(deriv.args[1:])
    a = diffs[x] if x in diffs else 0
    b = diffs[t] if t in diffs else 0
    
    return (I * k) ** a * (-I * w) ** b

def _linearise_product(prod, func):
    for idx, factor in enumerate(prod.args):
        if factor.is_Derivative:
            f_ab = prod / prod.args[idx]
            return f_ab.subs(func, C) * _get_C_from_deriv(factor)
        
    return factor.diff(func).subs(func, C)



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
