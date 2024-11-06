from sympy import symbols, apart, fraction, sympify, init_printing, div, Poly, Mul
from computation import cauer_division
import schemdraw as schem
import schemdraw.elements as elm 

s = symbols("s")

num_input = "s**4+5*s**2+4"
denom_input = "s*(s**2+2)"
cauer_num = "s*(s+2)*(s+6)"
cauer_denom = "2*(s+1)*(s+3)"
inv_cauer_num = "(1/s)*((1/s)+2)*((1/s)+6)"
inv_cauer_denom = "2*((1/s)+1)*((1/s)+3)"

lc_pairs = [2, 1, 1/3, 9/2, 1/6]
pair_count = 3
with schem.Drawing() as d:
    elm.Dot()
    for i in range(2):
        if i != 0:
            d.pop()
        elm.Line().right()
        if lc_pairs[i] != 0:
            elm.Inductor().right().label(lc_pairs[i])
        else:
            elm.Line().right()
        elm.Line().right()  
        d.push()
        elm.Line().down()
        if lc_pairs[i+1] != 0:
            elm.Capacitor().down().label(lc_pairs[i+1])
        else:
            elm.Swtich().down()
        elm.Line().down()
        elm.Line().left()
        elm.Line().left()       
        elm.Line().left()
        if i == 0:
            elm.Dot()
    d.draw()