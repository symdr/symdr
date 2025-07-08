from .utils import _get_C_from_deriv, _get_function_list, _get_deriv, is_continuous
from .utils import k, w, h, tau

from sympy import symbols, Function, exp, \
                  Integer, solveset, EmptySet, \
                  Eq, I, linear_eq_to_matrix, cos

from sympy import Derivative as D
from sympy import poly, LC

from .discrete_funcs import *


def equation_dr(expr):
    return system_dr([expr])[2]

"""
def system_dr(equations):
    fns = _get_function_list(equations)
    const = symbols(f"c0:{len(fns)}")
    amplitude = symbols(f"d0:{len(fns)}")
    prs   = list(zip(fns, const))

    def linearize(terms):
        expr_ln = len(terms)
        linerized_terms = []
        amplss = list(zip(fns, amplitude))
        prs   = list(zip(fns, const))

        ampls = {i[0]: i[1] for i in amplss}

        for term in range(expr_ln):
          d = list(terms[term].atoms(D))
          if (len(d) != 0):
            f = LC(poly(terms[term], d[0]))


            linerized_terms.append(
            ((I*k)**(d[0].variables.count(x)))*((-I*w)**d[0].variables.count(t)) \
            * f.subs(prs)* ampls[d[0].expr])

          else:

            linerized_terms.append(sum([ampls[fns[i]]*D(terms[term], fns[i]).subs(prs).doit() for i in range(len(fns))]))

        return sum(linerized_terms)
    lin_sys = []
    for equation in equations:
        terms = (equation.expand()).as_ordered_terms()
        lin_sys.append(linearize(terms))

    DR = linear_eq_to_matrix(lin_sys, amplitude)[0].det(method="lu")
    return (prs, linear_eq_to_matrix(lin_sys, amplitude)[0], DR)
"""

def system_dr(equations):
    func_list = _get_function_list(equations)
    const_list = symbols(f"c0:{len(func_list)}")
    amp_list = symbols(f"d0:{len(func_list)}")
    func_values = list(zip(func_list, const_list))
    func_amps = dict(zip(func_list, amp_list))

    def linearise(addends):
        addends_amount = len(addends)
        linearised_addends = []

        for addend in addends:
            derivatives = list(addend.atoms(D))
            has_derivative = len(derivatives) != 0

            if has_derivative:
                derivative = derivatives[0]
                factor = LC(poly(addend, derivative))
                x_order, t_order = derivative.variables.count(x), derivative.variables.count(t)
                factor_at_sol = factor.subs(func_values)
                func = derivative.expr
                addend = (I * k) ** x_order * (-I * w) ** t_order * factor_at_sol * func_amps[func]
            else:
                addend = sum(func_amps[func] * D(addend, func).subs(func_values).doit() for func in func_list)

            linearised_addends.append(addend)

        return sum(linearised_addends)

    linearised_equations = [linearise(equation.expand().as_ordered_terms()) for equation in equations]
    matrix = linear_eq_to_matrix(linearised_equations, amp_list)[0]
    disp_rel = matrix.det(method="lu")
    return func_values, matrix, disp_rel

def d_equation_dr(expr):
    return d_system_dr([expr])[2].rewrite(exp, cos).expand()

def d_system_dr(systems):
    order = len(systems)
    dfns = _get_function_list(systems)
    dconst     = symbols(f"c0:{len(dfns)}")
    damplitude = symbols(f"d0:{len(dfns)}")
    dprs   = list(zip(dfns, dconst))

    def linearize(terms):
        expr_ln = len(terms)
        linerized_terms = []
        amplss = list(zip(dfns, damplitude))
        ampls = {i[0]: i[1] for i in amplss}
        dc = {i[1]:i[0] for i in dprs}
        dct = {i[0]:i[1] for i in dprs}

        for i in range(len(terms)):
          if is_continuous(terms[i]):
            f = _get_deriv(terms[i].atoms(DiscreteGrid))

            if len(f):
              q = LC(poly(terms[i], f[0]))

              count_x, count_t = f[0].args[-2][1], f[0].args[-1][1]

              linerized_terms.append(
                      ((I*k)**count_x)*((-I*w)**count_t)*q.subs(dprs)*ampls[f[0].args[0].to_grid()]
                      )
            else:
                #return (terms[i], dprs)
                m = terms[i].subs(dprs)
                linerized_terms.append(sum([ampls[dc[k]]*m.diff(k) for k in dconst]))

          else:
            f = _get_deriv(terms[i].atoms(DiscreteGrid))
            q = LC(poly(terms[i], f[0]))


            ak = f[0].args[1] - a
            nk = f[0].args[2] - n

            linerized_terms.append(
                (q.subs(dprs)*(exp(I*(k*ak*h - w*nk*tau))) * \
                (ampls[f[0].args[0].to_grid()])).expand())




        return sum(linerized_terms)

    lin_sys = []
    for equation in systems:

        terms = (equation.expand()).as_ordered_terms()
        lin_sys.append(linearize(terms))


    DR = linear_eq_to_matrix(lin_sys, damplitude)[0].det(method="lu")
    return (dprs, linear_eq_to_matrix(lin_sys, damplitude)[0], DR)
