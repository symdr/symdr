from sympy.utilities.iterables import sift
from sympy.core.mul import Mul
from sympy.core.add import Add
from sympy.core.expr import Expr
from sympy.matrices.matrixbase import MatrixBase  # might be overkill
from sympy.printing.precedence import PRECEDENCE_TRADITIONAL

class SPHInclusiveSum(Expr):
    is_commutative = True
    
    def __new__(cls, expr, index):
        return super().__new__(cls, expr, index)

    @property
    def index(self):
        return self.args[1]

    @property
    def function(self):
        return self.args[0]

    @property
    def kind(self):
        return self.expr.kind

    @property
    def free_symbols(self):
        return self.function.free_symbols - {self.index}

    def _eval_factor(self, **hints):
        summand = self.function.factor(**hints)
        if summand.is_Mul:
            output = sift(summand.args, lambda w: w.is_commutative and not self.index in w.free_symbols)
            return Mul(*output[True]) * self.func(Mul(*output[False]), self.index)

        return self

    def _eval_expand_basic(self, **hints):
        summand = self.function.expand(**hints)
        if summand.is_Add:
            return Add(*[self.func(addend, self.index) for addend in summand.args])
        elif isinstance(summand, MatrixBase):
            return summand.applyfunc(lambda elem: self.func(elem, self.index))
        elif summand != self.function:
            return self.func(summand, self.index)
        return self

    def _latex(self, printer):
        tex_func = printer._print(self.function)
        if self.function.is_Add:
            tex_func = f"\left({tex_func}\right)"
        return r"\sum_{} {}".format(self.index, tex_func)  # TEST LATER

class SPHExclusiveSum(Expr):
    is_commutative = True
    
    def __new__(cls, expr, index, excluded):
        return super().__new__(cls, expr, index, excluded)

    @property
    def index(self):
        return self.args[1]

    @property
    def excluded(self):
        return self.args[2]

    @property
    def function(self):
        return self.args[0]

    @property
    def kind(self):
        return self.expr.kind

    @property
    def free_symbols(self):
        return self.function.free_symbols - {self.index}

    def _eval_factor(self, **hints):
        summand = self.function.factor(**hints)
        if summand.is_Mul:
            output = sift(summand.args, lambda w: w.is_commutative and not self.index in w.free_symbols)
            return Mul(*output[True]) * self.func(Mul(*output[False]), self.index, self.excluded)

        return self

    def _eval_expand_basic(self, **hints):
        summand = self.function.expand(**hints)
        if summand.is_Add:
            return Add(*[self.func(addend, self.index, self.excluded) for addend in summand.args])
        elif isinstance(summand, MatrixBase):
            return summand.applyfunc(lambda elem: self.func(elem, self.index, self.excluded))
        elif summand != self.function:
            return self.func(summand, self.index, self.excluded)
        return self

    def _latex(self, printer):
        tex_func = printer._print(self.function)
        if self.function.is_Add:
            tex_func = f"\left({tex_func}\right)"
        return r"\sum_{{} \ne {}} {}".format(self.index, self.excluded, tex_func)  # TEST LATER

PRECEDENCE_TRADITIONAL["SPHInclusiveSum"] = PRECEDENCE_TRADITIONAL["Sum"]
PRECEDENCE_TRADITIONAL["SPHExclusiveSum"] = PRECEDENCE_TRADITIONAL["Sum"]
