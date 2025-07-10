from sympy.utilities.iterables import sift
from sympy.core.mul import Mul
from sympy.core.add import Add
from sympy.core.expr import Expr
from sympy.matrices.matrixbase import MatrixBase

class SPHInclusiveSum(Expr):
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
        return self.function.free_symbols

    def _eval_factor(self, **hints):
        summand = self.function.factor(**hints)
        if summand.is_Mul:
            output = sift(summand.args, lambda w: w.is_commutative and not self.index in w.free_symbols)
            return Mul(*output[True]) * self.func(self.index, Mul(*output[False]))

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

