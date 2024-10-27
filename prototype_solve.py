'''
Caiman Moreno-Earle

The following is my implementation of a numerical solver for ordinary differential equations

Features Added
- Solve a single first order differential equation
- Solve a system of first order differential equations
- Plot the solution of the differential equation
- Plot the solution of the system of differential equations

Features to Add
    - Add support for higher order differential equations by converting them to a system of first order differential equations
    - create interpolation function to approximate the solution of the differential equation at any point
    - Add support for boundary value problems

Last Edit: 10-27-2024
'''

import re
import matplotlib.pyplot as plt
import math
import random

#miscellaneous mathematicla methods used b
class Misc_Methods:
        # Generator for non-integer range
        def frange(start, stop, step):
            while start < stop:
                yield start
                start += step

class PrototypeSolve:
    def __init__(self, equation, x, y , x0, y0, domain):
        
        #number of points used in approximation
        self.num_points = 100

        self.equation = equation
        self.x = x
        self.y = y
        self.x0 = x0
        self.y0 = y0
        self.domain = domain

        #time step for the runge-kutta method
        self.h = self.calculate_timestep()

        substring, f, dependent = self.parse_input(equation)
        self.function_definition = f
        self.dependent = dependent

        #define the function and the variables x and y
        self.function = self.string_to_function(substring[1])

        #points that solver has discretly solved for the integrated function at
        self.points = []

    def parse_input(self, equation):
        # Remove all spaces from the equation
        equation = equation.replace(' ', '')

        sub_strings = equation.split('=')
        derivative_pattern = re.compile(r"([a-zA-Z]+)'\(([^)]+)\)")
        
        for i, part in enumerate(sub_strings):
            match = derivative_pattern.search(part)
            if match:
                f = match.group(1)
                x = match.group(2)
                return sub_strings, f, x
        
        return sub_strings, None, None

    def string_to_function(self, func_str):
        # Define a function from a string
        exec(f"def func({self.x}, {self.y}): return {func_str}", globals())
        return eval('func') # type: ignore

    #select the time step for the runge-kutta method such that 1000 points are generated
    def calculate_timestep(self):
        return (self.domain[1] - self.domain[0]) / self.num_points

    #advance the equation by one step h using the runge-kutta method
    def RK4(self, x, y):

        k1 = self.function(x, y)
        k2 = self.function(x + self.h/2, y + self.h/2*k1)
        k3 = self.function(x + self.h/2, y + self.h/2*k2)
        k4 = self.function(x + self.h, y + self.h*k3)

        y1 = y + self.h/6*(k1 + 2*k2 + 2*k3 + k4)
        return y1


    def solve(self):
        #reset points to be empty
        self.points = []

        x = self.x0
        y = self.y0
        # Solve the equation for the given interval by solving step by step using runge-kutta method
        for _ in Misc_Methods.frange(self.domain[0], self.domain[1], self.h):
            self.points.append((x, y))
            y = self.RK4(x, y)
            x += self.h


    
    def plot_points(self):
        # Extract x and y coordinates from points
        x_coords = [point[0] for point in self.points]
        y_coords = [point[1] for point in self.points]

        # Plot the points
        plt.plot(x_coords, y_coords, linestyle='-', color='blue')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.title("Shitty Solve For " + self.function_definition+f"({self.x})")
        plt.show()




#Solve a system of differential equations and return the points of the solution for each function
class System_Solver:
    def __init__(self, equations, x, funcs , x0, initial_cons, domain):
        
        #number of points used in approximation
        self.num_points = 100

        self.equations = equations
        self.x = x
        self.funcs = funcs
        self.x0 = x0
        self.initial_cons = initial_cons
        self.domain = domain

        #time step for the runge-kutta method
        self.h = self.calculate_timestep()

        #array of functions representing the system of equations
        self.derivatives = []
        for i in range(len(equations)):
            substring, f, dependent = self.parse_input(equations[i])
            #substring 1 contains the equation fro derivative in terms of the variables
            self.derivatives.append(self.string_to_function(substring[1]))
            self.dependent = dependent
    

        #points that solver has discretly solved for the integrated function at
        self.points = []

    def parse_input(self, equation):
        # Remove all spaces from the equation
        equation = equation.replace(' ', '')

        sub_strings = equation.split('=')
        derivative_pattern = re.compile(r"([a-zA-Z]+)'\(([^)]+)\)")
        
        for i, part in enumerate(sub_strings):
            match = derivative_pattern.search(part)
            if match:
                f = match.group(1)
                x = match.group(2)
                return sub_strings, f, x
        
        
        return sub_strings, None, None

    def string_to_function(self, func_str):
        '''
        in a system the funtion is dependent on multiple variables, up to all the varibalbes in the system, 
        so include them all
        as potential arguments
        '''
        var_list = ', '.join(self.funcs)
        # Define a function from a string
        func_def = f"def func({self.x}, vars):\n"
        func_def += f"    {var_list} = vars\n"
        func_def += f"    return {func_str}\n"
        print(func_def)
        exec(func_def, globals())
        return eval('func')  # type: ignore

    #select the time step for the runge-kutta method such that 1000 points are generated
    def calculate_timestep(self):
        return (self.domain[1] - self.domain[0]) / self.num_points

    #increment the y variables based on step increase and the slope of derivative
    def increment_vars(self, vars, step, k):
        new_vars = []
        for i in range(len(vars)):
            new_vars.append(vars[i] + step*k[i])
        return new_vars

    #DON'T USE THIS FUNCTION UNLESS ABSOLUTELY NECCESARY AS YOU WILL MESS UP THE ALLIGNMENT OF Y and X
    #advance the equation by one step h using the runge-kutta method
    def RK4(self, derivative, var_index, x, vars):
        k1 = derivative(x, vars)
        k2 = derivative(x + self.h/2, self.increment_vars(vars, self.h/2*k1))
        k3 = derivative(x + self.h/2, self.increment_vars(vars,  self.h/2*k2))
        k4 = derivative(x + self.h, self.increment_vars(vars, self.h*k3))

        new_var = vars[var_index] + self.h/6*(k1 + 2*k2 + 2*k3 + k4)
        return new_var

    #performs rk4 to integrate a system of coupled first order diffiqs one step
    def system_RK4(self, x, vars):
        #array of slopes evaluated at start point
        k1 = [self.derivatives[i](x,vars) for i in range(len(vars))]
        #array of slopes evaluated at midpoint
        k2 = [self.derivatives[i](x+self.h/2, self.increment_vars(vars, self.h/2, k1)) for i in range(len(vars))]
        k3 = [self.derivatives[i](x+self.h/2, self.increment_vars(vars, self.h/2, k2)) for i in range(len(vars))]
        #array of slopes evaluated at end piint
        k4 = [self.derivatives[i](x+self.h, self.increment_vars(vars, self.h, k3)) for i in range(len(vars))]
        new_vars = [vars[i]+self.h/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(len(vars))]
        return new_vars

    def solve(self):
        #points for each solution are stored in diciotnary where the key is the function and the value is the list of points
        self.func_points = [[] for _ in range(len(self.funcs))]

        #starting time and initial conditions
        x_n = self.x0
        func_n = self.initial_cons
        func_n1 = [0 for _ in range(len(func_n))]

        # Solve the equation for the given interval by solving step by step using runge-kutta method
        for _ in Misc_Methods.frange(self.domain[0], self.domain[1], self.h):
            #integrate each of the equations at this timestep
            x_n += self.h
            func_n = self.system_RK4(x_n, func_n)
            #add ppints to the the plot to interpolate function later
            for i in range(len(self.funcs)):
                self.func_points[i].append((x_n, func_n[i]))

    #plot the points for one function
    def plot_points(self, points, func):
        # Extract x and y coordinates from points
        x_coords = [point[0] for point in points]
        y_coords = [point[1] for point in points]

        # Plot the points with a different color for each function
        plt.plot(x_coords, y_coords, linestyle='-', label=func)

    #plot all functions overlayed upon each other
    def plot_all_points(self):
        print(self.func_points)
        for i in range(len(self.funcs)):
            self.plot_points(self.func_points[i], self.funcs[i])

        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.title('Overlayed Plot of Functions')
        plt.legend()
        plt.grid()
        plt.show()




equation = "f'(x) = f*(0.2*math.cos(3.2*x)-0.2*math.sin(x))"
equations = [
    "y1'(x) = y2",
    "y2'(x) = -y1"
    ]

#solver = PrototypeSolve(equation, x='x', y='f', x0=0, y0=0.01, domain=(0, 10))
system_solver = System_Solver(equations, x='x', funcs=['y1', 'y2'], x0=0, initial_cons=[1,0], domain=(0, 10))
#solver.solve()
system_solver.solve()
#solver.plot_points()
system_solver.plot_all_points()

# result = solver.parse_input(equation)
# solver.solve(equation)
# print(result)  # Output: (['f\'(x) ', ' 2x + 3'], 'f', 'x')