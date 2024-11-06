from sympy import symbols, apart, fraction, sympify, init_printing, div, Poly, Mul

s = symbols("s")

num_input = "s**4+5*s**2+4"
denom_input = "s*(s**2+2)"
part_test_num = "(s+1)*(s+4)"
part_test_denom = "(s+3)*(s+5)**3*(s+1)"

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
print("Zeroes: ", zeroes)
# for term in (partial_fraction).as_ordered_terms():
#     term = term*s
#     print(term)