from sympy import symbols, apart, fraction, sympify, init_printing, div, Poly, Mul

s = symbols("s")

num_input = "s**4+5*s**2+4"
denom_input = "s*(s**2+2)"
part_test_num = "(s+1)*(s+4)"
part_test_denom = "(s+3)*(s+5)"

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

dividend = sympify("42*s")
divisor = sympify("6*s")
ans = div(dividend, divisor, domain='QQ')
quotient = ans[0]
remainder = ans[1]
print(remainder.degree())