import numpy as np
from datetime import datetime
from time import time

"""
A function that receives a coeff list
Returns the deriavtive 
"""
def derivative(p):
    # Creates an array from n to 1 and multiples by each coeff execpt last (which will be removed due to derivation)
    return p[:-1] * np.arange(len(p) - 1, 0, -1)

"""
A function that creates a list of powers and a list of complex representations of coeff 
uses them to calculate p(x) = an*x^n+.....+a1*x+a0 (SLIDE 4).
"""
def polyval(p, x):
    n = len(p) 
    x_pow = np.ones(n)*x
    # x^0 is always 1
    x_pow[0] = 1
    # Get the powers for x (each element is the product of the prev by itself)
    x_pow = np.multiply.accumulate(x_pow)
    # Reversing pow to match coeff array (so the powers are from larger to smaller)
    # Then activating dot product to calc the value of poly
    return (x_pow[::-1] @ p).item()

"""
This function simply divides value of poly at x and poly' at x.
Two cases here to avoid overflow for very large numbers.
Receives two lists and a complex num
""" 
def divide(p, q, x):
    if np.abs(x) <=1:
        return (polyval(p, x)/polyval(q, x))
    else:
        # Avoiding overflow -> same equation that is shown in slide 8.
        return x*(polyval(p[::-1], 1/x) / polyval(q[::-1], 1/x)) 
    
"""
Helper function for init roots, receives R and theta returns the complex representation
(SLIDE 4)
"""
def euler_equation(R, theta):
    real_part = R * np.cos(theta)
    imag_part = R * np.sin(theta)
    return real_part + 1j * imag_part

"""
Calculates the initial roots. Entire implementation of equation from SLIDE 4 
"""
def init_roots(p):
    # Number of roots depend on the degree of the poly, - 1 because the free number.
    num_of_roots = len(p) - 1 
    # Adding a  very small number to the division to avoid accidentally dividing by zero.
    R = 1 + max(np.abs(p[1:])) / (np.abs(p[0]) + 1e-10)
    roots = np.zeros(num_of_roots, dtype=complex)
    for i in range(num_of_roots):
        theta = (2*np.pi*i)/len(p)
        root = euler_equation(R, theta)
        roots[i] = root
    return roots

"""
This function calculates 1 / the diffrence between roots[k] and all the other roots (except k) 
"""
def sum_except_at_index(roots, k):
    # second part of the equation is to create an array that excludes the root at index k
    total_sum = np.sum(1 / (roots[k] - roots[np.arange(len(roots)) != k]))
    return total_sum

"""
SLIDE 5, entire equation.
Receives coeff list, derivative list, roots list and returns offsets list
"""
def getOffsets(p, p_tag, roots):
    #for each root in roots, apply divide function with p and ptag lists
    numerator = np.array([divide(p, p_tag, root) for root in roots])
    # call sum_except_at_index function for each root (at index k) in roots
    except_sums = np.array([sum_except_at_index(roots, k) for k in range(len(roots))])
    denominator = 1 - np.multiply(numerator, except_sums)
    # Complete the equation
    w_list = np.divide(numerator, denominator)
    return w_list

"""
Combining all the functions to eventually implement the alg and return the final roots
"""
def aberth_erlich(p, p_tag, epsilon, max_tries):
    tries = 0
    roots = init_roots(p)
    while True:
        w = getOffsets(p, p_tag, roots)
        maxW = max(w)
        print(f"maxW: {maxW}")
        # the alg stops when each offset is smaller than the defined epsilon
        # to avoid checking each offset, we simply take the max one and checking if it's smaller than epsilon
        if maxW < epsilon or tries >= max_tries:
            break
        # Updating roots with the offset -> SLIDE 6
        roots -= w
        tries += 1
    return roots

"""
Helper function that returns the context of the file input as a list of doubles
"""
def get_coefficient_from_file(file_name):
    with open(file_name, 'r') as file:
        return np.array([float(line.strip()) for line in file], dtype=np.float64)

"""
A function to print the results
"""
def format_roots_print(roots):
    formatted_roots = []
    for root in roots:
        if np.iscomplex(root):
            real_part = root.real
            imag_part = root.imag
            formatted_real = f"{real_part:.3f}"
            formatted_imag = f"{abs(imag_part):.3f}j"
            if imag_part < 0:
                formatted_imag = "-" + formatted_imag
            formatted_roots.append(formatted_real + " + " + formatted_imag)
        else:
            formatted_roots.append(f"{root:.3f}")
    print("Found roots:")
    for formatted_root in formatted_roots:
        print(formatted_root)

"""
This main function constructs the general flow of Abert Erlich algorithm 
Also manages parameters that help control when the algorithm reaches an end
"""
def main():
    max_tries = 850
    epsilon = 1e-5
    p = get_coefficient_from_file("poly_coeff(997).txt")
    p_tag = derivative(p)
    start_time = time()
    print(f"Starting the operation at: {str(datetime.now())}")
    roots = aberth_erlich(p, p_tag, epsilon, max_tries)
    end_time = time()
    elapsed_time = end_time - start_time
    format_roots_print(roots)
    print("The operation took:", elapsed_time, "seconds.")
    print(f"\nEnding the operation at: {str(datetime.now())}")

if __name__ == "__main__":
    main()