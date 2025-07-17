from .utils import _get_function_list, _get_deriv, is_continuous
from .utils import k, w, h, tau

from sympy import symbols, Function, exp, \
                  Integer, solveset, EmptySet, \
                  Eq, I, linear_eq_to_matrix, cos, Symbol, latex

from sympy import Derivative as D
from sympy import poly, LC

from .discrete_funcs import *
from .sph_sums import SPHInclusiveSum
from .sph import system_sph_dr, system_sph_dr_dimless


def equation_dr(expr):
    return system_dr([expr])[2]

def system_dr(equations):
    func_list = _get_function_list(equations)
    const_list = [Symbol(f"C_{func.func}") for func in func_list]
    amp_list = [Symbol(f"\\hat{{{latex(func.func)}}}") for func in func_list]
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
                factor_at_sol = factor.subs(func_values).doit()
                func = derivative.expr
                lin_addend = (I * k) ** x_order * (-I * w) ** t_order * factor_at_sol * func_amps[func]
            else:
                lin_addend = sum(func_amps[func] * D(addend, func).subs(func_values).doit() for func in func_list)

            linearised_addends.append(lin_addend)

        return sum(linearised_addends)

    linearised_equations = [linearise(equation.expand().as_ordered_terms()) for equation in equations]
    matrix = linear_eq_to_matrix(linearised_equations, amp_list)[0]
    disp_rel = matrix.det(method="lu")
    return func_values, matrix, disp_rel

def d_equation_dr(expr):
    return d_system_dr([expr])[2]

def d_system_dr(equations):
    grid_list = _get_function_list(equations)
    value_list = [Symbol(f"C_{grid.name}") for grid in grid_list]
    amp_list = [Symbol(f"\\hat{{{latex(grid.func)}}}") for grid in grid_list]
    grid_values = list(zip(grid_list, value_list))
    grid_amps = dict(zip(grid_list, amp_list))

    def linearise(addends):
        linearised_addends = []
        const_to_ampl = dict(zip(value_list, amp_list))

        for addend in addends:
            if is_continuous(addend):
                derivatives = _get_deriv(addend.atoms(DiscreteGrid))
                has_derivative = len(derivatives) != 0

                if has_derivative:
                    derivative = derivatives[0]
                    factor = LC(poly(addend, derivative))
                    order_x, order_t = derivative.args[-2][1], derivative.args[-1][1]
                    factor_at_sol = factor.subs(grid_values)
                    grid = derivative.args[0].to_grid()
                    lin_addend = (I * k) ** order_x * (-I * w) ** order_t * factor_at_sol * grid_amps[grid]
                else:
                    addend_at_sol = addend.subs(grid_values)
                    lin_addend = sum(const_to_ampl[const] * addend_at_sol.diff(const) for const in value_list)

            else:
                derivative = _get_deriv(addend.atoms(DiscreteGrid))[0]
                factor = LC(poly(addend, derivative))
                space_shift = (derivative.args[1] - a) * h
                time_shift = (derivative.args[2] - n) * tau
                factor_at_sol = factor.subs(grid_values)
                grid = derivative.args[0].to_grid()
                lin_addend = factor_at_sol * exp(I * (k * space_shift - w * time_shift)) * grid_amps[grid]

            linearised_addends.append(lin_addend)

        return sum(linearised_addends)

    linearised_equations = [linearise(equation.expand().as_ordered_terms()) for equation in equations]
    matrix = linear_eq_to_matrix(linearised_equations, amp_list)[0]
    disp_rel = matrix.det(method="lu").rewrite(exp, cos).expand()
    return grid_values, matrix, disp_rel
    
