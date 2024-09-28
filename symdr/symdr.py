from .utils import _get_C_from_deriv, _linearise_product, _get_function_list, _get_deriv, is_continuous
from .utils import k, w, U, C, h, tau, wave 

from sympy import symbols, Function, exp, \
                  Integer, solveset, EmptySet, \
                  Eq, I, linear_eq_to_matrix

from sympy import Derivative as D
from sympy import poly, LC

from .discrete_funcs import *


def equation_dr(expr):
    func = _get_function_list([expr])[0]
    
    if solveset(Eq(expr.subs(func, U).doit(), 0), U) == EmptySet:
        raise ValueError("Couldn't find constant solutions. More advanced solution check is not supported.")
    
    expr = expr.expand()
    linearised = Integer(0)

    for addend in expr.as_ordered_terms():
        if addend.is_Derivative:
            # print(addend, "->", get_C_from_deriv(addend))
            linearised += _get_C_from_deriv(addend)
        elif addend.is_Mul:
            linearised += _linearise_product(addend, func)
        else:
            linearised += addend.diff(u).subs(u, C)

    return linearised.expand()

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

            linerized_terms.append(
                sum([ampls[fns[i]]*D(terms[term], fns[i]).subs(prs).doit() for i in range(len(fns))]))

        return sum(linerized_terms)
    lin_sys = []
    for equation in equations:
        terms = (equation.expand()).as_ordered_terms ()
        lin_sys.append(linearize(terms))


    DR = linear_eq_to_matrix(lin_sys, amplitude)[0].det(method="lu")
    return (prs, linear_eq_to_matrix(lin_sys, amplitude)[0], DR)

def d_equation_dr(expr, trig_rewrite=False):
    functions = [base.to_grid() for base in expr.atoms(DiscreteGridBase)]
    if len(functions) != 1:
        raise ValueError("Equation has an invalid amount of functions")
    u = functions[0]

    linearised = 0

    for addend in expr.as_ordered_terms():
        if not is_continuous(addend):
            f_u, shifted = Integer(1), Integer(1)
            for factor in addend.as_ordered_factors():
                if factor.atoms(DiscreteGrid) == {u}:
                    f_u *= factor
                else:
                    shifted *= factor
            f_u = f_u.subs(u, C)
            shifted = shifted.subs(u, C + wave)
            linearised += f_u * shifted
        elif addend.atoms(DiscreteGrid) == {u}:
            linearised += addend.subs(u, C).diff(C) * wave
        else:
            for factor in addend.as_ordered_factors():
                if isinstance(factor, DiscreteGrid) and factor._is_diff:
                    C_ab = (I * k) ** factor.args[3][1] * (-I * w) ** factor.args[4][1]
                    linearised += C_ab * (addend / factor).subs(u, C) * wave
                    #print(C_ab * (addend / factor).subs(u, UC) * wave, 3)
                    break

    if trig_rewrite:
        return (linearised / wave).expand().rewrite(exp, cos).expand()
    else:
        return (linearised / wave).expand()

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
