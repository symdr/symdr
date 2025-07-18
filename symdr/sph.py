from .sph_sums import SPHInclusiveSum, SPHExclusiveSum
from .utils import k, w
from functools import reduce
import sympy

def _add_perturbations_and_expand(normal_funcs, stat_funcs, pert, equation, x_func, kernel_func, variable, global_idx):
    t = variable
    W = kernel_func
    
    lin_eq = 0

    time_diff_sub = {func(global_idx, t).diff(t): (-sympy.I * w) * pert[func](global_idx, t) for func in normal_funcs | {x_func}}
    global_lin_subs = {func(global_idx, t): stat_funcs[func] for func in normal_funcs}
    
    for addend in equation.as_ordered_terms():
        if addend.has(sympy.Derivative) and t in tuple(addend.atoms(sympy.Derivative))[0].variables:
            lin_addend = addend.subs(time_diff_sub)
        elif addend.has(SPHInclusiveSum, SPHExclusiveSum):
            order_wild, W_coeff_wild, W_arg_wild = sympy.Wild("lol"), sympy.Wild("loler"), sympy.Wild("lolest")

            sph_sum = tuple(addend.atoms(SPHInclusiveSum, SPHExclusiveSum))[0]
            idx = sph_sum.index
            coeff = addend.as_coefficient(sph_sum)
            excluded = sph_sum.excluded if sph_sum.func == SPHExclusiveSum else None

            under_sph_sum = coeff * sph_sum.function
            #addend_match = under_sph_sum.match(W_coeff_wild * sympy.Derivative(W_ab, (x_func(global_idx, t), order_wild)).doit())
            addend_match = under_sph_sum.match(W_coeff_wild * kernel_func(order_wild, W_arg_wild))
            
            order, W_coeff, W_arg = addend_match[order_wild], addend_match[W_coeff_wild], addend_match[W_arg_wild]

            dr_ab = pert[x_func](global_idx, t) - pert[x_func](idx, t)
            #new_W_ab = sympy.Derivative(W_ab, (x_func(global_idx, t), order)) + sympy.Derivative(W_ab, (x_func(global_idx, t), order + 1)) * dr_ab
            new_W_ab = kernel_func(order, W_arg) + kernel_func(order + 1, W_arg) * dr_ab

            lin_subs = global_lin_subs | {func(idx, t): stat_funcs[func] for func in normal_funcs}

            new_coeff = W_coeff.subs(lin_subs)
            new_coeff += sum(W_coeff.diff(func(global_idx, t)).subs(lin_subs) * pert[func](global_idx, t) for func in normal_funcs)
            new_coeff += sum(W_coeff.diff(func(idx, t)).subs(lin_subs) * pert[func](idx, t) for func in normal_funcs)

            if excluded is not None:
                lin_addend = sph_sum.func(new_coeff * new_W_ab, idx, excluded).expand()
            else:
                lin_addend = sph_sum.func(new_coeff * new_W_ab, idx).expand()
        else:
            new_addend = addend.subs(global_lin_subs)
            new_addend += sum(addend.diff(func(global_idx, t)).subs(global_lin_subs) * pert[func](global_idx, t) for func in normal_funcs)

            lin_addend = new_addend
            
        lin_eq += lin_addend

    return lin_eq

def _linearise_perturbations(normal_funcs, pert, lin_eq, x_func, kernel_func, variable, global_idx):
    t = variable
    W = kernel_func
    
    lin_eq2 = 0
    for addend in lin_eq.as_ordered_terms():
        if addend.has(SPHInclusiveSum, SPHExclusiveSum):
            sph_sum = tuple(addend.atoms(SPHInclusiveSum, SPHExclusiveSum))[0]
            expr, idx = sph_sum.function, sph_sum.index
            deltas = [pert[func](global_idx, t) for func in normal_funcs | {x_func}] + [pert[func](idx, t) for func in normal_funcs | {x_func}]
            poly = expr.as_poly(*deltas)
            lin_poly = sum(poly.coeff_monomial(delta) * delta for delta in deltas)
            if lin_poly != 0:
                lin_eq2 += addend.replace(expr, lin_poly)
        else:
            deltas = [pert[func](global_idx, t) for func in normal_funcs | {x_func}]
            poly = addend.as_poly(*deltas)
            lin_poly = sum(poly.coeff_monomial(delta) * delta for delta in deltas)
            lin_eq2 += lin_poly.as_expr()

    return lin_eq2

def _reduce_to_global_idx(normal_funcs, pert, lin_eq2, x_func, kernel_func, variable, global_idx):
    t = variable
    W = kernel_func

    lin_eq3 = 0
    for addend in lin_eq2.as_ordered_terms():
        if addend.has(SPHInclusiveSum, SPHExclusiveSum):
            sph_sum = tuple(addend.atoms(SPHInclusiveSum, SPHExclusiveSum))[0]
            expr, idx = sph_sum.function, sph_sum.index
            comp_exp = sympy.cos(k * (x_func(idx, t) - x_func(global_idx, t))) + sympy.I * sympy.sin(k * (x_func(idx, t) - x_func(global_idx, t)))
            sub_dict = {delta(idx, t): comp_exp * delta(global_idx, t) for delta in pert.values()}
            lin_eq3 += addend.subs(sub_dict).expand()
        else:
            lin_eq3 += addend

    return lin_eq3

def parity_wrt(expr, var):  # A hack that should be good enough for SPH
    # -1 if odd, 0 if neither odd nor even, 1 if even
    if expr.is_Add:  # REALLY hoping that nothing bad comes out like f(x) + f(-x)
        parities = set(parity_wrt(arg, var) for arg in expr.args)
        if parities == {-1} or parities == {1}:
            return tuple(parities)[0]
        else:
            return 0
    if expr.is_Mul:
        parities = tuple(parity_wrt(arg, var) for arg in expr.args)
        parity = reduce(lambda a, b: a * b, parities)
        return parity
    if expr.func == sympy.cos:
        return 1
    if expr.func == sympy.sin:
        return parity_wrt(expr.args[0], var)
    if expr.is_Pow:
        return 1 if expr.exp % 2 == 0 else parity_wrt(expr.base, var)
    if not expr.has(var):
        return 1
    if expr == var:
        return -1
    else:
        return 0

def _simplify_by_parity(lin_eq3, x_func, kernel_func, variable, global_idx):
    # if this will ever fail, substitute and integrate over symmetric domains.
    t = variable
    W = kernel_func
    
    lin_eq4 = 0
    for addend in lin_eq3.as_ordered_terms():
        if not addend.has(SPHInclusiveSum, SPHExclusiveSum):
            lin_eq4 += addend
            continue

        sph_sum = tuple(addend.atoms(SPHInclusiveSum, SPHExclusiveSum))[0]
        func = sph_sum.function
        is_even = True

        

        b = sph_sum.index
        parity_test_var = sympy.Symbol("f")
        if func.has(kernel_func):
            #deriv = tuple(func.atoms(sympy.Derivative))[0]
            #is_even = deriv.derivative_count % 2 == 0
            kernels = tuple(kernel for kernel in func.atoms(kernel_func) if b in kernel.args[1].free_symbols)  # most likely only 1 kernel
            assert len(kernels) <= 1
            if kernels:
                kernel = kernels[0]
                is_even = kernel.args[0] % 2 == 0
                kill_me, _ = func.as_independent(kernel)
            else:
                is_even = True
                kill_me = func
            kill_me = kill_me.subs(x_func(b, t), x_func(global_idx, t) - parity_test_var).subs(x_func(global_idx, t), 0)  # this is to get rid of x(b, t) - x(a, t)
            parity = parity_wrt(kill_me, parity_test_var)
            if parity == 0:
                is_even = False
            if parity == -1:
                is_even = not is_even
            #print(func, kill_me, is_even, parity)
            # tried is_constant here, but it kept saying a - b is constant wrt b
        #if func.find(sympy.sin(k * (x_func(global_idx, t) - x_func(b, t)))) or func.find(sympy.sin(k * (x_func(b, t) - x_func(global_idx, t)))) or \
        #   func.find(sympy.sin(k * x_func(global_idx, t) - k * x_func(b, t))) or func.find(sympy.sin(k * x_func(b, t) - k * x_func(global_idx, t))):
        #    is_even = not is_even

        if is_even:
            lin_eq4 += addend

    #print("fucking kill me", lin_eq4)

    return lin_eq4

def _prettify(lin_eq4):
    lin_eq5 = 0
    for addend in lin_eq4.as_ordered_terms():
        if not addend.has(SPHInclusiveSum, SPHExclusiveSum):
            lin_eq5 += addend
            continue

        sph_sum = tuple(addend.atoms(SPHInclusiveSum, SPHExclusiveSum))[0]
        func, idx = sph_sum.function, sph_sum.index

        indep, dep = func.as_independent(idx)
        lin_eq5 += addend.replace(func, dep) * indep

    return lin_eq5

def _linearise_one_singular_equation_man(normal_funcs, stat_funcs, pert, equation, x_func, kernel_func, variable, global_idx):
    lin_eq = _add_perturbations_and_expand(normal_funcs, stat_funcs, pert, equation, x_func, kernel_func, variable, global_idx)
    lin_eq2 = _linearise_perturbations(normal_funcs, pert, lin_eq, x_func, kernel_func, variable, global_idx)
    lin_eq3 = _reduce_to_global_idx(normal_funcs, pert, lin_eq2, x_func, kernel_func, variable, global_idx)
    lin_eq4 = _simplify_by_parity(lin_eq3, x_func, kernel_func, variable, global_idx)
    lin_eq5 = _prettify(lin_eq4)

    return lin_eq5

def _prettify_end(sph_dr):
    sph_dr = sph_dr.as_poly()
    sph_sums = sph_dr.atoms(SPHInclusiveSum, SPHExclusiveSum)
    coeff_list = sph_dr.all_coeffs()
    new_sph_dr = 0

    for idx, coeff in enumerate(reversed(coeff_list)):
        expr = coeff.collect(sph_sums)
        new_sum = 0
        for addend in expr.as_ordered_terms():
            if addend.has(SPHInclusiveSum, SPHExclusiveSum):
                sph_sum = tuple(addend.atoms(SPHInclusiveSum, SPHExclusiveSum))[0]
                felk = addend.as_independent(sph_sum)
                gcd = sympy.gcd(felk[0].as_ordered_terms())
                new_sum += felk[0].collect(gcd) * felk[1]
            else:
                gcd = sympy.gcd(addend.as_ordered_terms())
                new_sum += addend.collect(gcd)
        new_sph_dr += new_sum * w ** idx

    return new_sph_dr

def system_sph_dr(system, x_func, kernel_func, variable, global_idx, beautify=True):
    funcs = set().union(*[set(function.func for function in equation.atoms(sympy.Function)) for equation in system]) - {kernel_func}
    normal_funcs = set(funcs) - {x_func}
    
    pert = {func: sympy.Function(f"delta_{func.name}") for func in funcs}
    stat_funcs = {func: sympy.Symbol(f"{func.name}_c") for func in normal_funcs}

    lin_eqs = [_linearise_one_singular_equation_man(normal_funcs, stat_funcs, pert, equation, x_func, kernel_func, variable, global_idx) for equation in system]
    matrix = sympy.linear_eq_to_matrix(lin_eqs, [p(global_idx, variable) for p in pert.values()])[0]
    disp_rel = matrix.det().as_poly(w)

    #print(matrix)

    if beautify:
        disp_rel = _prettify_end(disp_rel)

    return stat_funcs, disp_rel

def system_sph_dr_dimless(system, x_func, density_func, kernel_func, mass_symbol, variable, global_idx, beautify=True):
    phi, K, H = sympy.symbols("phi K H", positive=True)
    W_dimless = sympy.Function("W_wave")

    local_idx_wild = sympy.Wild("b", properties=[lambda x: x != global_idx])
    
    stat_funcs, disp_rel = system_sph_dr(system, x_func, kernel_func, variable, global_idx, beautify=beautify)
    rho_c = stat_funcs[density_func]
    disp_rel = disp_rel.subs(mass_symbol, phi * H * rho_c)
    disp_rel = disp_rel.subs(k, 2 * sympy.pi / (K * H))
    disp_rel = disp_rel.replace(x_func(local_idx_wild, variable), x_func(global_idx, variable) - (global_idx - local_idx_wild) * phi * H)
    disp_rel = disp_rel.subs(x_func(global_idx, variable), 0)  # alledgedly they'll all be removed

    def dimless(expr):
        if not expr.args:
            return expr
        if expr.func == kernel_func:
            order, j = expr.args
            return W_dimless(order, j * phi) / H ** (order + 1)
        if expr.func == SPHInclusiveSum:
            idx = expr.index
            j = sympy.Symbol("j")
            func = expr.function.subs(idx, global_idx - j).subs(global_idx, 0)  # might cause problems
            func = dimless(func)
            zeroth = func.subs(j, 0)
            if zeroth == 0:
                return sympy.Sum(2 * func, (j, 0, sympy.floor(1 / phi)))
            else:
                return sympy.Sum(func, (j, -sympy.floor(1 / phi), sympy.floor(1 / phi)))
        if expr.func == SPHExclusiveSum:
            idx = expr.index
            j = sympy.Symbol("j")
            func = expr.function.subs(idx, global_idx - j).subs(global_idx, 0) # might cause problems
            return sympy.Sum(2 * dimless(func), (j, 1, sympy.floor(1 / phi)))
        return expr.func(*[dimless(arg) for arg in expr.args])

    disp_rel = dimless(disp_rel)

    return [W_dimless, phi, K, H], stat_funcs, disp_rel

#m = system_sph_dr_dimless(system_, x_, rho_, W_, m_, t_, a_)
