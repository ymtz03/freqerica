import sympy

class WrappedExpr(complex):
    def __new__(cls, expr):
        return super().__new__(cls, 'nan') # dummy
        #return super().__new__(cls, 0) # dummy
    
    def __init__(self, expr):
        if isinstance(expr, WrappedExpr):
            self.expr = expr.expr
        else:
            self.expr = sympy.Symbol(expr) if isinstance(expr, str) else expr

    def __repr__(self):
        return self.expr.__repr__()

    def __str__(self):
        return self.expr.__str__()

    def __bytes__(self):
        return self.expr.__bytes__()

    def __format__(self, format_spec):
        return self.expr.__format__(format_spec)

    def __bool__(self):
        try:
            return self.expr.__bool__()
        except TypeError:
            return False

    def subs(self, *args, **kwargs):
        return self.expr.subs(*args, **kwargs)
    

    def __lt__(self, other): # <
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__lt__(rhs))

    def __le__(self, other): # <=
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__le__(rhs))
    
    def __eq__(self, other): # ==
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__eq__(rhs))
    
    def __ne__(self, other): # !=
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__ne__(rhs))
    
    def __gt__(self, other): # >
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__gt__(rhs))
    
    def __ge__(self, other): # >=
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__ge__(rhs))

    
    def __add__(self, other):
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__add__(rhs))

    def __sub__(self, other):
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__sub__(rhs))

    def __mul__(self, other):
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__mul__(rhs))

    def __matmul__(self, other):
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__matmul__(rhs))

    def __truediv__(self, other):
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__truediv__(rhs))

    def __floordiv__(self, other):
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__floordiv__(rhs))

    def __mod__(self, other):
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__mod__(rhs))

    def __pow__(self, other, *modulo):
        rhs = other.expr if isinstance(other, WrappedExpr) else other
        if modulo and isinstance(modulo[0], WrappedExpr): modulo = (modulo[0].expr,)
        return WrappedExpr(self.expr.__pow__(rhs, *modulo))

    
    def __radd__(self, other):
        lhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__radd__(lhs))

    def __rsub__(self, other):
        lhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__rsub__(lhs))

    def __rmul__(self, other):
        lhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__rmul__(lhs))

    def __rmatmul__(self, other):
        lhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__rmatmul__(lhs))

    def __rtruediv__(self, other):
        lhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__rtruediv__(lhs))

    def __rfloordiv__(self, other):
        lhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__rfloordiv__(lhs))

    def __rmod__(self, other):
        lhs = other.expr if isinstance(other, WrappedExpr) else other
        return WrappedExpr(self.expr.__rmod__(lhs))

    
    def __neg__(self):
        return WrappedExpr(self.expr.__neg__())

    def __pos__(self):
        return WrappedExpr(self.expr.__pos__())

    def __abs__(self):
        return WrappedExpr(self.expr.__abs__())

    
import copy, copyreg
copyreg.pickle(WrappedExpr, lambda w: (WrappedExpr, (w.expr,)))

from sympy.core.sympify import converter
converter[WrappedExpr] = lambda x:x


def _test():
    s1 = sympy.Symbol('spam')
    s2 = sympy.Symbol('egg')
    w1 = WrappedExpr(s1)
    w2 = WrappedExpr(s2)
    
    print(isinstance(w1, float))
    print(hash(s1))
    print()
    
    print(w1+w2)
    print(type(w1+w2))
    print()

    print(w1+100)
    print(type(w1+100))
    print(type((w1+100).expr))
    print()

    print(100+w1)
    print(type(100+w1))
    print()

    print(100.+w1)
    print(type(100.+w1))
    print()

    c = 100+2j
    print(type(c))
    print(c+w1)
    print(type(c+w1))

    print(w1-w2)
    print(w1-w1)
    print('aaa', pow(w1, w2))
    print('aaa', pow(w1, w2, 2))
    print('aaa', pow(w1, w2, w1))


    print(copy.copy(w1))
    print(copy.deepcopy(w1))


    import openfermion
    
    fop = openfermion.FermionOperator('1^ 2', w1)
    print(fop)
    print(fop+fop)
    fop+=fop*2+openfermion.FermionOperator('',1)
    print(fop)
