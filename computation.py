from sympy import symbols, apart, fraction, sympify, init_printing
import schemdraw as schem
import schemdraw.elements as elm 
import userinterface as ui
from math import ceil

init_printing()

s = symbols("s")

optno = 0

def synthesize(num_raw, den_raw):
    numerator, denominator = sympify(num_raw), sympify(den_raw)
    zeroes = numerator.as_poly(s).all_roots()
    print("Zeroes: ", zeroes)
    poles = denominator.as_poly(s).all_roots()
    print("Poles: ", poles)
    partial_fraction = apart(numerator / denominator)
    #for debugging purposes
    print("\nPartial fraction:", partial_fraction,"\n")
    
    # Assign each option that we get to the specific function eg. foster1_computing_function
    flag = False
    if optno == 0:
        for term in partial_fraction.as_ordered_terms():
            num, denom = fraction(term)
            if denom.as_poly(s).degree() >= 2:
                foster1_lc(partial_fraction)
                flag = True
        if not flag:
            if zeroes[-1] < poles[-1]:
                foster1_rc(partial_fraction)

            if zeroes[-1] >  poles[-1]:
                modded_f1_fracs = apart(numerator / (denominator*s))
                foster1_rl(modded_f1_fracs)
    elif optno == 1:
        for term in partial_fraction.as_ordered_terms():
            num, denom = fraction(term)
            if denom.as_poly(s).degree() >= 2:
                foster2_lc(partial_fraction)
                break
        else:
            if zeroes[-1] > poles[-1]:
                modded_f2_fracs = apart(numerator / (denominator*s))
                foster2_rc(modded_f2_fracs)
            if zeroes[-1] < poles[-1]:
                foster2_rl(partial_fraction)
    elif optno == 2:
        for term in partial_fraction.as_ordered_terms():
            num, denom = fraction(term)
            if denom.as_poly(s).degree() >= 2:
                cauer1_lc(numerator, denominator)
                break
        else:
            if zeroes[-1] < poles[-1]:
                for term in partial_fraction.as_ordered_terms():
                    if term.is_negative:
                        cauer1_rc(numerator, denominator*s)
                        break
                else: 
                    cauer1_rl(numerator, denominator)
            if zeroes[-1] >  poles[-1]:
                for term in partial_fraction.as_ordered_terms():
                    if term.is_negative:
                        cauer1_rl(numerator, denominator*s)
                        
                else:
                    cauer1_rc(numerator, denominator)
    elif optno == 3:
        for term in partial_fraction.as_ordered_terms():
            num, denom = fraction(term)
            if denom.as_poly(s).degree() == 2:
                cauer2_lc(numerator, denominator)
                break
        else:
            if zeroes[-1] > poles[-1]:
                for term in partial_fraction.as_ordered_terms():
                    if term.is_negative:
                        cauer2_rc(numerator, denominator*s)
                else: 
                    cauer2_rl(numerator, denominator)
            if zeroes[-1] < poles[-1]:
                for term in partial_fraction.as_ordered_terms():
                    if term.is_negative:
                        cauer2_rl(numerator, denominator*s)
                else:
                    cauer2_rc(numerator, denominator)    
    else:
        print("Invalid option")
        exit()

def foster1_lc(partial_fraction):
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
    print("calling drawing function")
    foster1_lc_df(c0, l_inf, lc_pairs)

def foster2_lc(partial_fraction):
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
    foster2_lc_df(l0, c_inf, lc_pairs)

def cauer1_lc(num, denom):
    lc_terms = []
    
    # recursive function to find the lc value
    def find_l_c_value(dividend, divisor):
        print(dividend.degree(), divisor.degree())
        if dividend.degree() >= 1:
            if dividend.degree() > divisor.degree():
                quotient, remainder = cauer_division(dividend, divisor)
                lc_terms.append(quotient.as_poly(s).coeff_monomial(s))
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
    cauer1_lc_df(lc_terms)
    
def cauer2_lc(num, denom):
    cl_terms = []
    
    # recursive function to find the l value
    def find_c_l_value(dividend, divisor):
        print(dividend.degree(), divisor.degree())
        if dividend.degree() >= 1:
            if dividend.degree() > divisor.degree():
                quotient, remainder = cauer_division(dividend, divisor, True)
                cl_terms.append(1/(quotient.as_poly(1/s).coeff_monomial(1/s)))
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
        find_c_l_value(dividend=poly_num, divisor=poly_denom)
    except PolynomialDegreeZeroException as e:
        print(e)
        
    cauer2_cl_df(cl_terms)

def foster1_rc(partial_fraction):
    c0 = 0 
    res = 0
    rc_pairs = []
    
    for term in partial_fraction.as_ordered_terms():
        num, denom = fraction(term)
        
        if denom.is_number and num.is_number:
            res = num/denom
            print("Res: ", res)
        if num.is_number and denom.is_polynomial(s):
            denom_poly = denom.as_poly(s)
            coeffs = denom_poly.all_coeffs()
    
            # The polynomial shd be of the form as + 0
            if len(coeffs) == 2 and coeffs[-1] == 0 and coeffs[0] != 0:
                c0 = coeffs[0]/num
                print("C0: ", c0)
                
        if num.is_number and denom.is_polynomial(s):
            denom_poly = denom.as_poly(s)
            coeffs = denom_poly.all_coeffs()
            
            # The polynomial shd be of the form as + b
            if len(coeffs) == 2 and coeffs[-1] != 0 and coeffs[0] != 0:
                
                k = num/coeffs[0]
                sigma = coeffs[1]/coeffs[0]

                rc_pairs.append((k/sigma, 1/k))
                print("RC pairs: ", rc_pairs)
    foster1_rc_df(c0, res, rc_pairs)
                
def foster1_rl(partial_fraction):
    print(f"foster partial_fraction: {partial_fraction}")
    l0 = 0 
    res = 0
    rl_pairs = []
    
    for term in partial_fraction.as_ordered_terms():
        term = term*s
        num, denom = fraction(term)
        
        if denom.is_number and num.is_number:
            res = num/denom
            print("Res: ", res)
        if denom.is_number and num.is_polynomial(s):
            num_poly = num.as_poly(s)
            coeffs = num_poly.all_coeffs()
    
            # The polynomial shd be of the form as + 0
            if len(coeffs) == 2 and coeffs[-1] == 0 and coeffs[0] != 0:
                l0 = coeffs[0]/denom
                print("L0: ", l0)
                
        if num.is_polynomial(s) and denom.is_polynomial(s):
            denom_poly = denom.as_poly(s)
            denom_coeffs = denom_poly.all_coeffs()
            num_poly = num.as_poly(s)
            num_coeffs = num_poly.all_coeffs()
            
            # The polynomial shd be of the form as / as+b
            if len(denom_coeffs) == 2 and denom_coeffs[-1] != 0 and denom_coeffs[0] != 0:
                if len(num_coeffs) == 2 and num_coeffs[-1] == 0 and num_coeffs[0] != 0:
                    k = num_coeffs[0]/denom_coeffs[0]
                    sigma = denom_coeffs[1]/denom_coeffs[0]

                    rl_pairs.append((k, k/sigma))
                    print("RL pairs: ", rl_pairs)
    foster1_rl_df(l0, res, rl_pairs)
    
def foster2_rc(partial_fraction):
    c0 = 0 
    res = 0
    rc_pairs = []
    
    for term in partial_fraction.as_ordered_terms():
        term = term*s
        num, denom = fraction(term)
        
        if denom.is_number and num.is_number:
            res = denom/num
            print("Res: ", res)
        if denom.is_number and num.is_polynomial(s):
            num_poly = num.as_poly(s)
            coeffs = num_poly.all_coeffs()
    
            # The polynomial shd be of the form as + 0
            if len(coeffs) == 2 and coeffs[-1] == 0 and coeffs[0] != 0:
                c0 = coeffs[0]/denom
                print("c0: ", c0)
                
        if num.is_polynomial(s) and denom.is_polynomial(s):
            denom_poly = denom.as_poly(s)
            denom_coeffs = denom_poly.all_coeffs()
            num_poly = num.as_poly(s)
            num_coeffs = num_poly.all_coeffs()
            
            # The polynomial shd be of the form as / as+b
            if len(denom_coeffs) == 2 and denom_coeffs[-1] != 0 and denom_coeffs[0] != 0:
                if len(num_coeffs) == 2 and num_coeffs[-1] == 0 and num_coeffs[0] != 0:
                    k = num_coeffs[0]/denom_coeffs[0]
                    sigma = denom_coeffs[1]/denom_coeffs[0]

                    rc_pairs.append((1/k, k/sigma))
                    print("RC pairs: ", rc_pairs)
    foster2_rc_df(c0, res, rc_pairs)

def foster2_rl(partial_fraction):
    l0 = 0 
    res = 0
    rl_pairs = []
    
    for term in partial_fraction.as_ordered_terms():
        num, denom = fraction(term)
        
        if denom.is_number and num.is_number:
            res = denom/num
            print("Res: ", res)
        if num.is_number and denom.is_polynomial(s):
            denom_poly = denom.as_poly(s)
            coeffs = denom_poly.all_coeffs()
    
            # The polynomial shd be of the form as + 0
            if len(coeffs) == 2 and coeffs[-1] == 0 and coeffs[0] != 0:
                l0 = coeffs[0]/num
                print("L0: ", l0)
                
        if num.is_number and denom.is_polynomial(s):
            denom_poly = denom.as_poly(s)
            coeffs = denom_poly.all_coeffs()
            
            # The polynomial shd be of the form as + b
            if len(coeffs) == 2 and coeffs[-1] != 0 and coeffs[0] != 0:
                
                k = num/coeffs[0]
                sigma = coeffs[1]/coeffs[0]

                rl_pairs.append((sigma/k, 1/k))
                print("RL pairs: ", rl_pairs)

def cauer1_rc(num, denom):
    rc_terms = []
    
    # recursive function to find the l value
    def find_r_c_value(dividend, divisor):  
        print(dividend.degree(), divisor.degree())
        if dividend.degree() >= 1:
            if dividend.degree() >= divisor.degree():
                quotient, remainder = cauer_division(dividend, divisor)
                
                if quotient.is_number:
                    print("Quotient: ", quotient)
                    rc_terms.append(quotient)
                else:
                    rc_terms.append(quotient.as_poly(s).coeff_monomial(s))
                    
                print("RC terms: ", rc_terms)
                find_r_c_value(dividend=divisor, divisor=remainder)
            else:
                rc_terms.append(0)
                print("RC terms: ", rc_terms)

                find_r_c_value(dividend=divisor, divisor=dividend)
        elif dividend.degree() == 0:
            rc_terms.append(dividend/divisor)
            raise PolynomialDegreeZeroException("Dividend degree is 0")
        
    numerator, denominator = num.as_poly(s), denom.as_poly(s)
    try:
        find_r_c_value(dividend=numerator, divisor=denominator)
    except PolynomialDegreeZeroException as e:
        print(e)   
    print("RC terms: ", rc_terms)
    cauer1_rc_df(rc_terms)

def cauer1_rl(num, denom):
    rl_terms = []
    
    # recursive function to find the l value
    def find_r_l_value(dividend, divisor):
        print(dividend.degree(), divisor.degree())
        print("n Dividend: ", dividend)
        print("n Divisor: ", divisor)
        if dividend.degree() > 0:
            if dividend.degree() >= divisor.degree():
                quotient, remainder = cauer_division(dividend, divisor)
                print(f"remainder: {remainder}")
                if quotient.is_number:
                    print("Quotient: ", quotient)
                    rl_terms.append(1/quotient)
                else:
                    rl_terms.append(quotient.as_poly(s).coeff_monomial(s))
                print("RL terms: ", rl_terms)
                find_r_l_value(dividend=divisor, divisor=remainder)
            else:
                rl_terms.append(0)
                print("RL terms: ", rl_terms)

                find_r_l_value(dividend=divisor, divisor=dividend)
        elif dividend.degree() == 0:
            rl_terms.append(divisor/dividend)
            print("RL terms: ", rl_terms)
            raise PolynomialDegreeZeroException("Dividend degree is 0")
        
    numerator, denominator = num.as_poly(s), denom.as_poly(s)
    try:
        find_r_l_value(dividend=numerator, divisor=denominator)
    except PolynomialDegreeZeroException as e:
        print(e)   
    cauer1_rl_df(rl_terms)

def cauer2_rc(num, denom):
    rc_terms = []
    
    # recursive function to find the l value
    def find_r_c_value(dividend, divisor):
        print(dividend.degree(), divisor.degree())
        if dividend.degree() >= 1:
            if dividend.degree() >= divisor.degree():
                quotient, remainder = cauer_division(dividend, divisor, True)
                if quotient.is_number:
                    print("Quotient: ", quotient)
                    rc_terms.append(1/quotient)
                else:
                    rc_terms.append(1/quotient.as_poly(1/s).coeff_monomial(1/s))
                print("RC terms: ", rc_terms)
                find_r_c_value(dividend=divisor, divisor=remainder)
            else:
                rc_terms.append(0)
                print("RC terms: ", rc_terms)
                find_r_c_value(dividend=divisor, divisor=dividend)
        elif dividend.degree() == 0:
            rc_terms.append(divisor/dividend)
            raise PolynomialDegreeZeroException("Dividend degree is 0")
    
    
    numerator, denominator = num.as_poly(s), denom.as_poly(s)
    highest_degree = max(numerator.degree(), denominator.degree())
    print("Highest degree: ", highest_degree) 
    
    new_num = numerator.as_expr() * (s**(-highest_degree))
    new_denom = denominator.as_expr() * (s**(-highest_degree))

    poly_num = new_num.as_poly(1/s)
    poly_denom = new_denom.as_poly(1/s)

    try:
        find_r_c_value(dividend=poly_num, divisor=poly_denom)
    except PolynomialDegreeZeroException as e:
        print(e)
    print("RC terms: ", rc_terms)
    cauer2_rc_df(rc_terms)

def cauer2_rl(num, denom):
    rl_terms = []
    
    # recursive function to find the l value
    def find_r_l_value(dividend, divisor):
        print(dividend.degree(), divisor.degree())
        if dividend.degree() >= 1:
            if dividend.degree() >= divisor.degree():
                quotient, remainder = cauer_division(dividend, divisor, True)
                if quotient.is_number:
                    print("Quotient: ", quotient)
                    rl_terms.append(quotient)
                else:
                    rl_terms.append(1/quotient.as_poly(1/s).coeff_monomial(1/s))
                print("RL terms: ", rl_terms)
                find_r_l_value(dividend=divisor, divisor=remainder)
            else:
                rl_terms.append(0)
                print("dividend degree less than divisor! RL terms: ", rl_terms)
                find_r_l_value(dividend=divisor, divisor=dividend)
        elif dividend.degree() == 0:
            rl_terms.append(dividend/divisor)
            raise PolynomialDegreeZeroException("Dividend degree is 0")
    
    
    numerator, denominator = num.as_poly(s), denom.as_poly(s)
    highest_degree = max(numerator.degree(), denominator.degree())
    print("Highest degree: ", highest_degree) 
    
    new_num = numerator.as_expr() * (s**(-highest_degree))
    new_denom = denominator.as_expr() * (s**(-highest_degree))

    poly_num = new_num.as_poly(1/s)
    poly_denom = new_denom.as_poly(1/s)

    try:
        find_r_l_value(dividend=poly_num, divisor=poly_denom)
    except PolynomialDegreeZeroException as e:
        print(e)
    print("RL terms: ", rl_terms)
    cauer2_rl_df(rl_terms)

class PolynomialDegreeZeroException(Exception):
    pass

def cauer_division(dividend, divisor, is2=False): #takes in two polynomials
    print(f"dividend: {dividend}")
    print(f"divisor: {divisor}")
    if is2:
        dividend_monom = dividend.as_expr().as_ordered_terms()[-1]
        divisor_monom = divisor.as_expr().as_ordered_terms()[-1] 
    else:
        dividend_monom = dividend.as_expr().as_ordered_terms()[0]
        divisor_monom = divisor.as_expr().as_ordered_terms()[0]
    print(f"dividend_monom: {dividend_monom} divisor_monom: {divisor_monom}")
    cou_quotient = dividend_monom/divisor_monom
    subtrahend = divisor*cou_quotient
    cou_remainder = dividend - subtrahend
    
    return cou_quotient, cou_remainder
    # quotient: int, remainder: polyn

def foster1_lc_df(c0, l_inf, lc_pairs):
    c0_label = f"C0 = {c0}F"
    l_inf_label = f"L_inf = {l_inf}H"    
    pair_count = len(lc_pairs)

    with schem.Drawing() as d:
        elm.Dot()
        elm.Line().right()
        if c0 != 0:
            elm.Capacitor().label(c0_label).right()
        else:
            elm.Line().right()
        
        for i in range(pair_count):
            elm.Line().up()
            if lc_pairs[i][0] != 0:
                elm.Inductor().right().label(f"{lc_pairs[i][0]}H")
            else:
                elm.Line().right()
            elm.Line().down()
            d.push()
            elm.Line().down()
            if lc_pairs[i][1] != 0:
                elm.Capacitor().left().label(f"{lc_pairs[i][1]}F")
            else: 
                elm.Line().left()
            elm.Line().up()
            d.pop()
            elm.Line().right()
        
        elm.Line().right()
        elm.Line().down()
        if l_inf != 0:
            elm.Inductor().down().label(l_inf_label)
        else:
            elm.Line().down()
        elm.Line().down()
        for i in range(pair_count+4):
            elm.Line().left()
            
        elm.Dot()
        d.draw()

def foster2_lc_df(l0, c_inf, lc_pairs):
    l0_label = f"{l0}H"
    c_inf_label = f"{c_inf}F"    
    pair_count = len(lc_pairs)

    
    with schem.Drawing() as d:
        elm.Dot()
        elm.Line().right()
        elm.Line().right()  
        d.push()
        elm.Line().down()
        if l0 != 0:
            elm.Inductor().down().label(l0_label)
        else:
            elm.Switch().down()
        elm.Line().down()
        elm.Line().left()
        elm.Line().left()
        elm.Dot()
        
        if c_inf != 0:
            d.pop()
            elm.Line().right()
            elm.Line().right()  
            d.push()
            elm.Line().down()
            elm.Capacitor().down().label(c_inf_label)
            elm.Line().down()
            elm.Line().left()
            elm.Line().left()        
        
        for i in range(pair_count):
            d.pop()
            elm.Line().right()
            elm.Line().right()  
            d.push()
            elm.Inductor().down().label(f"{lc_pairs[i][0]}H")
            elm.Line().down()
            elm.Capacitor().down().label(f"{lc_pairs[i][1]}F")
            elm.Line().left()
            elm.Line().left()       
        
        d.draw()
        
def cauer1_lc_df(lc_pairs):  
    pair_count = ceil(len(lc_pairs)/2)
    
    with schem.Drawing() as d:
        elm.Dot()
        for i in range(pair_count):
            if i != 0:
                d.pop()
            elm.Line().right()
            if lc_pairs[2*i] != 0:
                elm.Inductor().right().label(f"{lc_pairs[2*i]}H")
            else:
                elm.Line().right()
            elm.Line().right()  
            d.push()
            elm.Line().down()
            if (2*i+1)<len(lc_pairs) and lc_pairs[2*i+1] != 0:
                elm.Capacitor().down().label(f"{lc_pairs[2*i+1]}F")
            else:
                elm.Line().down()
            elm.Line().down()
            elm.Line().left()
            elm.Line().left()       
            elm.Line().left()
            if i == 0:
                elm.Dot()        
        d.draw()
                
def cauer2_cl_df(cl_pairs):  
    pair_count = ceil(len(cl_pairs)/2)
        
    with schem.Drawing() as d:
        elm.Dot()
        for i in range(pair_count):
            if i != 0:
                d.pop()
            elm.Line().right()
            if cl_pairs[2*i] != 0:
                elm.Capacitor().right().label(f"{cl_pairs[i]}F")
            else:
                elm.Line().right()
            elm.Line().right()  
            d.push()
            elm.Line().down()
            if (2*i+1)<len(cl_pairs) and cl_pairs[2*i+1] != 0:
                elm.Inductor().down().label(f"{cl_pairs[2*i+1]}H")
            else:
                elm.Line().down()
            elm.Line().down()
            elm.Line().left()
            elm.Line().left()       
            elm.Line().left()
            if i == 0:
                elm.Dot()        
        d.draw()

def foster1_rc_df(c0, res, rc_pairs):
    c0_label = f"{c0}F"
    res = f"{res}\N{GREEK CAPITAL LETTER OMEGA}"   
    pair_count = len(rc_pairs)

    with schem.Drawing() as d:
        elm.Dot()
        elm.Line().right()
        if c0 != 0:
            elm.Capacitor().label(c0_label).right()
        else:
            elm.Line().right()
        
        for i in range(pair_count):
            elm.Line().up()
            if rc_pairs[i] != 0:
                elm.ResistorIEEE().right().label(f"{rc_pairs[i][0]}\N{GREEK CAPITAL LETTER OMEGA}")
            else:
                elm.Line().right()
            elm.Line().down()
            d.push()
            elm.Line().down()
            if rc_pairs[i][1] != 0:
                elm.Capacitor().left().label(f"{rc_pairs[i][1]}F")
            else:
                elm.Line().left()
            elm.Line().up()
            d.pop()
            elm.Line().right()
        
        elm.Line().right()
        elm.Line().down()
        if res != 0:
            elm.ResistorIEEE().down().label(res)
        else:
            elm.Line().down()
        elm.Line().down()
        for i in range(pair_count+5):
            elm.Line().left()
            
        elm.Dot()
        d.draw()

def foster1_rl_df(l0, res, rl_pairs):
    l0_label = f"{l0}H"
    res_label = f"{res}\N{GREEK CAPITAL LETTER OMEGA}"    
    pair_count = len(rl_pairs)
    
    with schem.Drawing() as d:
        elm.Dot()
        elm.Line().right()
        elm.Line().right()  
        d.push()
        elm.Line().down()
        if l0 != 0:
            elm.Inductor().down().label(l0_label)
        else:
            elm.Switch().down()
        elm.Line().down()
        elm.Line().left()
        elm.Line().left()
        elm.Dot()
        
        if res != 0:
            d.pop()
            elm.Line().right()
            elm.Line().right()  
            d.push()
            elm.Line().down()
            elm.ResistorIEEE().down().label(res_label)
            elm.Line().down()
            elm.Line().left()
            elm.Line().left()        
        
        for i in range(pair_count):
            d.pop()
            elm.Line().right()
            elm.Line().right()  
            d.push()
            elm.Inductor().down().label(f"{rl_pairs[i][0]}H")
            elm.Line().down()
            elm.ResistorIEEE().down().label(f"{rl_pairs[i][1]}\N{GREEK CAPITAL LETTER OMEGA}")
            elm.Line().left()
            elm.Line().left()       
        
        d.draw()
        
def foster2_rc_df(c0, res, rc_pairs):
    c0_label = f"{c0}F"
    res_label = f"{res}\N{GREEK CAPITAL LETTER OMEGA}"    
    pair_count = len(rc_pairs)
    
    with schem.Drawing() as d:
        elm.Dot()
        elm.Line().right()
        elm.Line().right()  
        d.push()
        elm.Line().down()
        if c0 != 0:
            elm.Capacitor().down().label(c0_label)
        else:
            elm.Switch().down()
        elm.Line().down()
        elm.Line().left()
        elm.Line().left()
        elm.Dot()
        
        if res != 0:
            d.pop()
            elm.Line().right()
            elm.Line().right()  
            d.push()
            elm.Line().down()
            elm.ResistorIEEE().down().label(res_label)
            elm.Line().down()
            elm.Line().left()
            elm.Line().left()        
        
        for i in range(pair_count):
            d.pop()
            elm.Line().right()
            elm.Line().right()  
            d.push()
            elm.ResistorIEEE().down().label(f"{rc_pairs[i][0]}H")
            elm.Line().down()
            elm.Capacitor().down().label(f"{rc_pairs[i][1]}\N{GREEK CAPITAL LETTER OMEGA}")
            elm.Line().left()
            elm.Line().left()       
        
        d.draw()
            
def foster2_rl_df(l0, res, rl_pairs):
    {}
    
def cauer1_rc_df(rc_pairs):
    pair_count = ceil(len(rc_pairs)/2)
    
    with schem.Drawing() as d:
        elm.Dot()
        for i in range(3):
            elm.Line().right()
            d.push()
        for i in range(pair_count):
            
            elm.Line().down()
            if rc_pairs[2*i] != 0:
                elm.Capacitor().down().label(f"{rc_pairs[2*i]}F")
            else:
                elm.Line().down()
            elm.Line().down()
            elm.Line().left()
            elm.Line().left()       
            elm.Line().left()
            if i == 0:
                elm.Dot()        
                
            d.pop()
            elm.Line().right()
            if (2*i+1)<len(rc_pairs) and rc_pairs[2*i+1] != 0:
                elm.ResistorIEEE().right().label(f"{rc_pairs[2*i+1]}\N{GREEK CAPITAL LETTER OMEGA}")
            else:
                elm.Line().right()
            elm.Line().right()  
            d.push()
            
        d.draw()
    
def cauer1_rl_df(rl_pairs):
    pair_count = ceil(len(rl_pairs)/2)
    
    with schem.Drawing() as d:
        elm.Dot()
        for i in range(3):
            elm.Line().right()
            d.push()
        for i in range(pair_count):
            
            elm.Line().down()
            if rl_pairs[2*i] != 0:
                elm.ResistorIEEE().down().label(f"{rl_pairs[2*i]}\N{GREEK CAPITAL LETTER OMEGA}")
            else:
                elm.Line().down()
            elm.Line().down()
            elm.Line().left()
            elm.Line().left()       
            elm.Line().left()
            if i == 0:
                elm.Dot()        
                
            d.pop()
            elm.Line().right()
            if (2*i+1)<len(rl_pairs) and rl_pairs[2*i+1] != 0:
                elm.Inductor().right().label(f"{rl_pairs[2*i+1]}H")
            else:
                elm.Line().right()
            elm.Line().right()  
            d.push()
            
        d.draw()

def cauer2_rc_df(rc_pairs):
    pair_count = ceil(len(rc_pairs)/2)
    
    with schem.Drawing() as d:
        elm.Dot()
        for i in range(pair_count):
            if i != 0:
                d.pop()
            elm.Line().right()
            if rc_pairs[2*i] != 0:
                elm.Capacitor().right().label(f"{rc_pairs[2*i]}F")
            else:
                elm.Line().right()
            elm.Line().right()  
            d.push()
            elm.Line().down()
            if (2*i+1)<len(rc_pairs) and rc_pairs[2*i+1] != 0:
                elm.Resistor().down().label(f"{rc_pairs[2*i+1]}\N{GREEK CAPITAL LETTER OMEGA}")
            else:
                elm.Line().down()
            elm.Line().down()
            elm.Line().left()
            elm.Line().left()       
            elm.Line().left()
            if i == 0:
                elm.Dot()        
        d.draw()
        
def cauer2_rl_df(rl_pairs):
    pair_count = ceil(len(rl_pairs)/2)
    
    with schem.Drawing() as d:
        elm.Dot()
        for i in range(pair_count):
            if i != 0:
                d.pop()
            elm.Line().right()
            if rl_pairs[2*i] != 0:
                elm.ResistorIEEE().right().label(f"{rl_pairs[2*i]}\N{GREEK CAPITAL LETTER OMEGA}")
            else:
                elm.Line().right()
            elm.Line().right()  
            d.push()
            elm.Line().down()
            if (2*i+1)<len(rl_pairs) and rl_pairs[2*i+1] != 0:
                elm.Capacitor().down().label(f"{rl_pairs[2*i+1]}F")
            else:
                elm.Line().down()
            elm.Line().down()
            elm.Line().left()
            elm.Line().left()       
            elm.Line().left()
            if i == 0:
                elm.Dot()        
        d.draw()