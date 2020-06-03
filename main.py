import math
from sympy import *
from sympy.parsing.sympy_parser import *
import matplotlib.pyplot as plt
import numpy as np

# _____________________________________
num = [1]
den = [1, 125, 5100, 65000, 0]
poles_real = [0, -25, -50, -50]
poles_complex = [0, 0, 10, -10]
# _____________________________________
EPS = 0.00000000001
n = len(den) - 1
m = len(num) - 1


class Function:
    def __init__(self, function):
        """the function is a valid string representing the function  """
        """the function is assumed to be a function of x """

        transformation = standard_transformations + (implicit_application,)
        self.f = parse_expr(function.replace('-0.0', '0'), transformations=transformation)
        self.__x = symbols('x')

    def get_root(self):
        return solve(self.f, self.__x)

    def get_derivative_root(self):
        #  roots = solve(self.f.diff(self.__x), self.__x)
        # print(self.f.diff(self.__x))
        return solve(self.f.diff(self.__x), self.__x)

    def get_value_at(self, x):
        """ get the value of the function at x """
        ans = self.f.subs(self.__x, x)
        return ans.evalf()

    def get_derivative_value_at(self, x):
        derivative = self.f.diff(self.__x)
        ans = derivative.subs(self.__x, x)
        return ans.evalf()

    def get_expr(self):
        return self.f


def get_funtion_string():
    function_string = "("
    for i in range(0, len(den)):
        function_string += str(den[i]) + "* x**" + str(len(den) - i - 1)
        if i < len(den) - 1:
            function_string += "+"
    function_string += ")"
    return function_string


def get_angles():
    res = set(())
    for i in range(0, n - m):
        res.add((((int)(180 * (2 * i + 1) / (n - m))) % 360 + 360) % 360)
        res.add((((int)(-180 * (2 * i + 1) / (n - m))) % 360 + 360) % 360)

    return res


def get_centroid():
    sum_poles = sum(poles_real)
    return sum_poles / (n - m)


def get_breaking_down_s():
    function_string = get_funtion_string()
    function = Function(function_string)
    roots_of_derivatives = function.get_derivative_root()
    res = list()
    function = Function("-1 *" + function_string)
    for root in roots_of_derivatives:
        try:
            val = eval(str(function.get_value_at(root)))
            if (val > 0):
                # print(root.evalf())
                res.append(root.evalf())
        except:
            pass

    return res


def get_angle_of_depature():
    res = []
    sz = len(den)
    for i in range(0, n):
        sum = 0
        if poles_complex[i] != 0:
            for j in range(0, n):
                if j == i:
                    continue
                else:
                    denominator = poles_real[i] - poles_real[j]
                    if abs(denominator) < EPS:
                        if (poles_complex[i] - poles_complex[j] > 0):
                            sum += 90
                        else:
                            sum += 270
                    else:
                        sum += math.atan(
                            (poles_complex[i] - poles_complex[j]) / (poles_real[i] - poles_real[j])) * 180 / math.pi
            res.append(((180 - sum) % 360 + 360) % 360)
        else:
            res.append(None)
    return res


def get_intersection_with_img_axis():
    sz = len(den)

    raws = n
    cols = 1 + (int)(sz / 2)
    table = [["0" for i in range(cols)] for j in range(raws)]
    k = 0
    for i in range(0, sz, 2):
        table[0][k] = "(" + str(den[i]) + ")"
        k += 1
    k = 0
    for i in range(1, sz, 2):
        table[1][k] = "(" + str(den[i]) + ")"
        k += 1
    if (sz % 2 == 1):
        if eval(table[0][-1]) != 0:
            table[0][-1] += 'x'
        else:
            table[0][-1] = 'x'
    else:
        if eval(table[1][-1]) != 0:
            table[1][-1] += 'x'
        else:
            table[1][-1] = 'x'
    # steps of routh
    cols = len(table[0])
    for i in range(2, n):
        for j in range(0, cols - 1):
            current = "("
            current += "(" + str(table[i - 1][0]) + "*" + table[i - 2][j + 1] + "-" + table[i - 2][0] + "*" + \
                       table[i - 1][j + 1] + ")"
            current += "/" + str(table[i - 1][0])
            current += ")"
            try:
                val = eval(current)
                table[i][j] = "(" + str(val) + ")"
            except:
                table[i][j] = current
            # table[i].append( (table[i-1][0] "*" table[i-2][j+1] "+" table[i-2][0] "*" table[i-1][j+1])/table[i-1][0])

    # print(table)
    equation = table[-1][0]
    function = Function(equation)
    x = function.get_root()[0]
    table[-2][1] = eval(table[-2][1])
    w = math.sqrt(eval(str(table[-2][1]) + "/" + str(table[-2][0])))
    return [w, -w]


def plot():
    # Draw the real locus by changing the values of K
    xs = []
    ys = []
    s = get_funtion_string()
    for k in range(0, 5000000, 50000):
        f = Function(s + "+" + str(k))
        l = f.get_root()
        for root in l:
            c = complex(root)
            xs.append(c.real)
            ys.append(c.imag)
    plt.plot(xs, ys, 'r*', mew=0.05)

    # Drawing axis and asymptotes
    x = get_centroid()
    val = 50
    plt.plot([val, -val], [0, 0], '--', lw=1.1 , color='#BFBFBF')
    plt.plot([0, 0], [val, -val], '--', lw=1.1 , color='#BFBFBF')
    plt.plot([val, -val], [val + abs(x), -(val - abs(x))], '--', lw=1.3, color='#BFBFBF')
    plt.plot([val, -val], [-(val + abs(x)), (val - abs(x))], '--', lw=1.3, color='#BFBFBF')


    # draw manually the root locus
    xs = [-50, -70]
    ys = [10, 30]
    plt.plot(xs, ys , color = 'black')

    xs = [-50, -70]
    ys = [-10, -30]
    plt.plot(xs, ys ,  color = 'black')

    xs = [-25, 0]
    ys = [0, 0]
    plt.plot(xs, ys, color='black')

    x = complex(get_breaking_down_s()[0]).real
    i = x
    xs = []
    ys = []
    #x = 0.018275 y ^ 2 - 9.1503
    while(i < 35):
        xs.append(i)
        ys.append(sqrt(i - x ) * 6.3)
        i+= 0.05

    plt.plot(xs , ys , 'k-')

    i = x
    xs = []
    ys = []
    # x = 0.018275 y ^ 2 - 9.1503
    while (i < 35):
        xs.append(i)
        ys.append(-sqrt(i - x) * 6.3)

        i += 0.05

    plt.plot(xs, ys, 'k-')


    plt.grid()
    plt.show()

print("---------The results  of the program-----")
print("The equation you enter is ", get_funtion_string())
print("angles of asymptotes are ", get_angles())
print("Astoroid is ", get_centroid())
print("S which leads to a break down is ", get_breaking_down_s()[0].args[0])
print("angles of departure for the given poles are ", get_angle_of_depature())
print("The Intersection with the img axis occurs at ", get_intersection_with_img_axis())
plot()
