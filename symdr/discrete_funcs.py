from sympy import *
from sympy.printing.latex import LatexPrinter
from sympy.core.expr import Basic, Expr
from sympy.core.containers import Tuple
from sympy import Idx



a, n = symbols("a n", cls=Idx)
x, t = symbols("x t")
tau, h = symbols("tau h")


def is_continuous(expr):
    if isinstance(expr, DiscreteGrid):
        return expr.is_continious
    elif isinstance(expr, (Integer, Symbol)):
        return True
    else:
        return all(is_continuous(arg) for arg in expr.args)


class DiscreteLatexPrinter(LatexPrinter):
    def __init__(self, settings=None):
        LatexPrinter.__init__(self, settings)
        self.ignore_grid = False

    def _print_Pow(self, expr):
        if isinstance(expr.base, DiscreteGrid):
            return "({%s})^{%s}" % (self._print(expr.base), self._print(expr.exp))
        return super()._print_Pow(expr)

    def _print_DiscreteGrid(self, expr):
        base_name = self.doprint(expr.args[0])
        x_displ = self.doprint(expr.args[1])
        t_displ = self.doprint(expr.args[2])
        x_diff_order = expr.args[3]
        t_diff_order = expr.args[4]

        diff = self.doprint(x) * x_diff_order[1] + \
               self.doprint(t) * t_diff_order[1]  # ???
        if expr._is_shifted:
            return f"[{base_name}]" + "_{" + x_displ + "}^{" + t_displ + "}"
        if not self.ignore_grid:
            if expr._is_diff:
                return f"[{base_name}_" + "{" + diff + "}]_a^n"
            else:
                return f"[{base_name}]_a^n"
        else:
            if expr._is_diff:
                return f"{base_name}_" + "{" + diff + "}"
            else:
                return base_name

    def _print_Add(self, expr):
        if not expr.atoms(DiscreteGrid) or self.ignore_grid:
            return super()._print_Add(expr)

        cont_add_list = []
        disc_add_list = []

        for arg in expr.args:
            # print(arg)
            if is_continuous(arg) and arg.atoms(DiscreteGrid):
                # print("continuous")
                cont_add_list.append(arg)
            else:
                # print("discrete")
                disc_add_list.append(arg)

        cont_add = Add(*cont_add_list, evaluate=False)
        disc_add = Add(*disc_add_list, evaluate=False)

        self.ignore_grid = True
        cont_str = f"\left[{super()._print_Add(cont_add)}\\right]_a^n"
        self.ignore_grid = False
        disc_str = super()._print_Add(disc_add)

        if len(cont_add_list) == 1 and isinstance(cont_add_list[0], DiscreteGrid):
            cont_str = self._print_DiscreteGrid(cont_add_list[0])

        if not cont_add_list:
            return disc_str
        elif not disc_add_list:
            return cont_str
        else:
            if disc_str[0] == '-':
                return cont_str + disc_str
            return cont_str + ' + ' + disc_str

    def _print_Mul(self, expr):
        if not expr.atoms(DiscreteGrid) or self.ignore_grid:
            return super()._print_Mul(expr)

        cont_mul_list = []
        disc_mul_list = []

        for arg in expr.args:
            # print(arg)
            if is_continuous(arg):
                # print("continuous")
                cont_mul_list.append(arg)
            else:
                # print("discrete")
                disc_mul_list.append(arg)

        cont_mul = Mul(*cont_mul_list, evaluate=False)
        disc_mul = Mul(*disc_mul_list, evaluate=False)

        self.ignore_grid = True
        cont_str = f"\left[{super()._print_Mul(cont_mul)}\\right]_a^n"
        self.ignore_grid = False
        disc_str = super()._print_Mul(disc_mul)

        if not cont_mul.atoms(DiscreteGrid):
            return super()._print_Mul(expr)

        if len(cont_mul_list) == 1 and isinstance(cont_mul_list[0], Add):
            cont_str = f"\left({cont_str}\\right)"
        if len(disc_mul_list) == 1 and isinstance(disc_mul_list[0], Add):
            disc_str = f"\left({disc_str}\\right)"

        if not cont_mul_list:
            return disc_str
        elif not disc_mul_list:
            return cont_str
        else:
            return cont_str + disc_str

def discrete_latex_printer(expr, **kwargs):
    return DiscreteLatexPrinter(kwargs).doprint(expr)

#init_printing(latex_printer=discrete_latex_printer)
# TODO: add a latex printing check

class DiscreteGridBase(Basic):
    def __new__(cls, name):
        return super().__new__(cls, Symbol(name))

    def _latex(self, printer):
        return printer.doprint(self.args[0])

    def to_grid(self):
        return DiscreteGrid(self, a, n, (x, 0), (t, 0))

class DiscreteGrid(Expr):
    is_commutative = True

    def __new__(cls, desig, x_shift=a, t_shift=n,
                x_diff_order=Tuple(x,0), t_diff_order=Tuple(t,0)):
        base = DiscreteGridBase(desig) if isinstance(desig, str) else desig
        return super().__new__(cls, base, x_shift, t_shift,
                               x_diff_order, t_diff_order)

    @property
    def _is_diff(self):
        return self.args[3][1] != S.Zero or self.args[4][1] != S.Zero

    @property
    def _is_shifted(self):
        return self.args[1] != a or self.args[2] != n

    @property
    def is_continious(self):
        return not self._is_shifted

    @property
    def is_discrete(self):
        return self._is_shifted

    def diff(self, *diffs):
        if self._is_shifted:
            raise TypeError("Cannot differentiate a shifted grid")

        x_count = 0
        t_count = 0

        for differ in diffs:
            if differ == x:
                x_count += 1
            elif differ == t:
                t_count += 1
            else:
                raise ValueError(f"Unknown differetial {differ}")

        new_arg3 = Tuple(x, self.args[3][1] + x_count)
        new_arg4 = Tuple(t, self.args[4][1] + t_count)

        return DiscreteGrid(self.args[0], self.args[1], self.args[2],
                            new_arg3, new_arg4)

    def at(self, x_shift, t_shift):
        if self._is_diff:
            raise TypeError("Cannot shift a differentiated grid")

        return DiscreteGrid(self.args[0], x_shift, t_shift,
                            self.args[3], self.args[4])

    def at_x(self, x_shift):
        if self._is_diff:
            raise TypeError("Cannot shift a differentiated grid")

        return DiscreteGrid(self.args[0], x_shift, self.args[2],
                            self.args[3], self.args[4])

    def at_t(self, t_shift):
        if self._is_diff:
            raise TypeError("Cannot shift a differentiated grid")

        return DiscreteGrid(self.args[0], self.args[1], t_shift,
                            self.args[3], self.args[4])

    _t = at_t
    _x = at_x

    def _eval_subs(self, old, new):
        if self._is_shifted:
            shifted_x = x + (self.args[1] - a) * h
            shifted_t = t + (self.args[2] - n) * tau
            return new.subs({x: shifted_x, t: shifted_t})
        elif self._is_diff:
            return Derivative(new, self.args[3], self.args[4])
        else:
            return new

