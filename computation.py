from sympy import symbols, apart, fraction, sympify, init_printing, div, Poly, Mul
import schemdraw as schem
import schemdraw.elements as elm 
import userinterface as ui
from math import ceil

init_printing()


s = symbols("s")

optno = 0

def synthesize(num_raw, den_raw):
    numerator, denominator = sympify(num_raw), sympify(den_raw)
    partial_fraction = apart(numerator / denominator)
    #for debugging purposes
    print("\nPartial fraction:", partial_fraction,"\n")
    
    # Assign each option that we get to the specific function eg. foster1_computing_function
    if optno == 0:
        foster1_cf(partial_fraction)
    elif optno == 1:
        foster2_cf(partial_fraction)
    elif optno == 2:
        cauer1_cf(numerator, denominator)
    elif optno == 3:
        cauer2_cf(numerator, denominator)
    else:
        print("Invalid option")
        exit()

def foster1_cf(partial_fraction):
    c0 = 0
    l_inf = 0
    lc_pairs = []
    
    for term in partial_fraction.as_ordered_terms():
        num, denom = fraction(term)
        
        if denom.is_number and num.is_polynomial(s) and num.as_poly(s).degree() == 1:
            l_inf = num.as_poly(s).coeff_monomial(s)/denom
            print("l_inf: ", l_inf)
            
        if num.is_number and denom.is_polynomial(s):
            denom_poly = denom.as_poly(s)
            if denom_poly.degree() == 1:
                c0 = (denom.as_poly(s).all_coeffs()[0])/num
                print("C0: ", c0)
        
        if denom.is_polynomial(s) and denom != 1:
            denom_poly = denom.as_poly(s)
            if denom_poly and denom_poly.degree() == 2 and denom_poly.coeff_monomial(s) == 0:
                coeffs = denom_poly.all_coeffs()
                
                k = num.as_poly(s).coeff_monomial(s)/coeffs[0]
                sigma=  coeffs[2]/coeffs[0]
                
                l_value = k/sigma
                c_value = 1/k
                lc_pairs.append((l_value, c_value))
                print("LC pairs: ", lc_pairs)           
    foster1_df(c0, l_inf, lc_pairs)

def foster2_cf(partial_fraction):
    l0 = 0
    c_inf = 0
    lc_pairs = []
    
    for term in partial_fraction.as_ordered_terms():
        num, denom = fraction(term)
        
        if denom.is_number and num.is_polynomial(s) and num.as_poly(s).degree() == 1:
            c_inf = num.as_poly(s).coeff_monomial(s)/denom
            print("c_inf: ", c_inf)
            
        if num.is_number and denom.is_polynomial(s):
            denom_poly = denom.as_poly(s)
            if denom_poly.degree() == 1:
                l0 = (denom.as_poly(s).all_coeffs()[0])/num
                print("l0: ", l0)
        
        if denom.is_polynomial(s) and denom != 1:
            denom_poly = denom.as_poly(s)
            if denom_poly and denom_poly.degree() == 2 and denom_poly.coeff_monomial(s) == 0:
                coeffs = denom_poly.all_coeffs()
                
                k = num.as_poly(s).coeff_monomial(s)/coeffs[0]
                sigma=  coeffs[2]/coeffs[0]
                
                l_value = 1/k
                c_value = k/sigma
                lc_pairs.append((l_value, c_value))
                print("LC pairs: ", lc_pairs)
    return l0, c_inf, lc_pairs

def cauer1_cf(num, denom):
    lc_terms = []
    
    # recursive function to find the l value
    def find_l_c_value(dividend, divisor):
        print(dividend.degree(), divisor.degree())
        if dividend.degree() >= 1:
            if dividend.degree() > divisor.degree():
                ans = div(dividend, divisor, domain='QQ')
                quotient = ans[0]
                remainder = ans[1]
                print("Quotient: ", quotient)
                lc_terms.append(quotient.coeff_monomial(s))
                print("LC terms: ", lc_terms)
                find_l_c_value(dividend=divisor, divisor=remainder)
            else:
                lc_terms.append(0)
                print("LC terms: ", lc_terms)

                find_l_c_value(dividend=divisor, divisor=dividend)
        elif dividend.degree() == 0:
            raise PolynomialDegreeZeroException("Dividend degree is 0")
        
    numerator, denominator = num.as_poly(s), denom.as_poly(s)
    try:
        find_l_c_value(dividend=numerator, divisor=denominator)
    except PolynomialDegreeZeroException as e:
        print(e)   
    return lc_terms
    
def cauer2_cf(num, denom):
    cl_terms = []
    
    # recursive function to find the l value
    def find_c_l_value(dividend, divisor):
        print(dividend.degree(), divisor.degree())
        if dividend.degree() >= 1:
            if dividend.degree() > divisor.degree():
                ans = div(dividend, divisor, domain='QQ')
                quotient = ans[0]
                remainder = ans[1]
                print("Quotient: ", quotient)
                cl_terms.append(1/(quotient.coeff_monomial(1/s)))
                print("CL terms: ", cl_terms)
                find_c_l_value(dividend=divisor, divisor=remainder)
            else:
                cl_terms.append(0)
                print("CL terms: ", cl_terms)

                find_c_l_value(dividend=divisor, divisor=dividend)
        elif dividend.degree() == 0:
            raise PolynomialDegreeZeroException("Dividend degree is 0")
    
    
    numerator, denominator = num.as_poly(s), denom.as_poly(s)
    highest_degree = max(numerator.degree(), denominator.degree())
    print("Highest degree: ", highest_degree) 
    
    new_num = numerator.as_expr() * (s**(-highest_degree))
    new_denom = denominator.as_expr() * (s**(-highest_degree))

    poly_num = new_num.as_poly(1/s)
    poly_denom = new_denom.as_poly(1/s)

    try:
        find_c_l_value(dividend=poly_denom, divisor=poly_num)
    except PolynomialDegreeZeroException as e:
        print(e)
        
    return cl_terms

class PolynomialDegreeZeroException(Exception):
    pass

# # Add a resistor
# d.add(elm.Resistor().label('R'))

# # Add an inductor
# d.add(elm.Inductor().right().label('L'))



def foster1_df(c0, l_inf, lc_pairs):
    c0_label = f"C0 = {c0}F"
    l_inf_label = f"L_inf = {l_inf}H"    
    pair_count = ceil(len(lc_pairs)/2)

    with schem.Drawing() as d:
        elm.Dot()
        elm.Line().right()
        elm.Capacitor().label(c0_label).right()
        
        for i in range(pair_count):
            elm.Line().up()
            elm.Inductor().right().label(f"{lc_pairs[i][0]}H")
            elm.Line().down()
            d.push()
            elm.Line().down()
            elm.Capacitor().left().label(f"{lc_pairs[i][1]}F")
            elm.Line().up()
            d.pop()
            elm.Line().right()
        
        elm.Line().right()
        elm.Line().down()
        elm.Inductor().down().label(l_inf_label)
        elm.Line().down()
        for i in range(pair_count+4):
            elm.Line().left()
            
        elm.Dot()   
        d.draw()
    
