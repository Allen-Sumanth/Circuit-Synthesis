from sympy import symbols, apart, fraction, sympify, init_printing, div, Poly, Mul
from computation import cauer_division

s = symbols("s")

num_input = "s**4+5*s**2+4"
denom_input = "s*(s**2+2)"
cauer_num = "s*(s+2)*(s+6)"
cauer_denom = "2*(s+1)*(s+3)"
inv_cauer_num = "(1/s)*((1/s)+2)*((1/s)+6)"
inv_cauer_denom = "2*((1/s)+1)*((1/s)+3)"

cauer_num_simp, cauer_denom_simp = sympify(inv_cauer_num), sympify(inv_cauer_denom)
cau_num_poly, cau_denom_poly = cauer_num_simp.as_poly(1/s), cauer_denom_simp.as_poly(1/s)
print(cau_num_poly)
print(cau_denom_poly)
# num, denom = sympify(num_input), sympify(denom_input)
# numerator, denominator = num.as_poly(s), denom.as_poly(s)
# highest_degree = max(numerator.degree(), denominator.degree())
# print("Highest degree: ", highest_degree) 

# new_num = numerator.as_expr() * (s**(-highest_degree))
# new_denom = denominator.as_expr() * (s**(-highest_degree))

# print("New numerator: ", new_num)
# print("New denominator: ", new_denom)

# poly_num = new_num.as_poly(1/s)
# poly_denom = new_denom.as_poly(1/s)
# print("Poly numerator: ", poly_num)
# print("Poly denominator: ", poly_denom)

# print("Poly numerator degree: ", poly_num.degree())

# ans = div(poly_num, poly_denom, domain='QQ')
# quotient = ans[0]
# remainder = ans[1]
# print("Quotient: ", quotient)

# print(quotient.coeff_monomial(1/s))

numerator, denominator = sympify("x**2+2*x+1"), sympify("(x+1)*(x+4)")
zeroes = numerator.as_poly(s).all_roots()
poles = denominator.as_poly(s).all_roots()
# partial_fraction = apart(numerator / (denominator*s))
quo, rem = cauer_division(cau_num_poly, cau_denom_poly, True)
print(quo)
print(rem)
# for term in (partial_fraction).as_ordered_terms():
#     term = term*s
#     print(term)