import numpy as np
import matplotlib.pyplot as plt

def generate_exact_functions(u1, u2, plot = False):
    N = 10000 #can be made higher for higher precision
    grid = np.linspace(0, 3, 10000)
    u1_ex = u1(grid)
    u2_ex = u2(grid)
    if plot:
        plt.plot(grid, u1_ex, label="u1_ex", color='green')
        plt.plot(grid, u2_ex, label="u2_ex", color='red')
        plt.xlabel("x")
        plt.ylabel("u_ex(x)")
        plt.title("u(x)")
        plt.legend()
        f = plt.gcf()
        plt.show()
    return u1_ex, u2_ex

def solveP(N, plotting, u1ex, u2ex, f1, f2, u1, u2):
    #define grid
    h = 3 / N
    xgrid = np.linspace(0, 3, N + 1)

    #generate BC vector
    g = np.zeros_like(xgrid[1:-1])
    g[0] = 1
    g[-1] = 1

    #generate Laplace FD matrix
    A = np.zeros((N-1, N-1))
    np.fill_diagonal(A[:-1, 1:], -1)
    np.fill_diagonal(A[1:, :-1], -1)
    np.fill_diagonal(A, 2)
    if plotting:
        #output the numerical eigenvalues
        print("Eigenvalues are:")
        print(np.linalg.eig(A * 1/ h **2)[0])
        plt.plot(xgrid, f1(xgrid), label = "f1")
        plt.plot(xgrid, f2(xgrid), label = "f2")
        plt.xlabel("x")
        plt.ylabel("f(x)")
        plt.title("f(x)")
        plt.legend()
        f  = plt.gcf()
        plt.show()
        plt.cla()
        #imshow the matrix
        plt.imshow(A)
        plt.colorbar()
        plt.show()
    #compute forcing terms:
    f_1 = f1(xgrid[1:-1]) * h ** 2 + g
    f_2 = f2(xgrid[1:-1]) * h ** 2 + g

    #now solving system Au = f_1
    u_1 = np.linalg.solve(A, f_1)
    u1n = np.ones_like(xgrid)
    u1n[1:-1] = u_1

    u_2 = np.linalg.solve(A, f_2)
    u2n = np.ones_like(xgrid)
    u2n[1:-1] = u_2

    if plotting:

        plt.plot(xgrid, u1n, linestyle="--", color="blue", label = 'u1_FD')
        plt.plot(xgrid, u1(xgrid),  label = 'u1_ex', color = 'green')
        plt.plot(xgrid, u2n, linestyle="--", color="orange", label = 'u2_FD')
        plt.plot(xgrid, u2(xgrid),  label = 'u2_ex', color = 'red')
        plt.xlabel("X")
        plt.ylabel("u(x)")
        plt.title("u(x) for N = 5   ")
        plt.legend()
        plt.show()

    #now calculate the maximum error:
    output1 = np.zeros_like(u1ex)
    output2 = np.zeros_like(u2ex)
    err1 = np.sum(np.abs(u1n - u1(xgrid))**2)
    err1 = np.sqrt(err1 / (N-1))
    err2 = np.sum(np.abs(u2n - u2(xgrid)) ** 2)
    err2 = np.sqrt(err2 / (N - 1))
    return err1, err2

def convergence_study():
    output1, output2 = [], []
    for x in range(5, 1000):
        err1, err2 = solveP(x, False, u1_ex, u2_ex, f1, f2, u1, u2)
        output1.append(np.log(err1))
        output2.append(np.log(err2))
    #plt.plot(np.log(np.linspace(5, 1000, 995)), output1, label = "err1")
    plt.plot(np.log(np.linspace(5, 1000, 995)), output2, label = "err2", color = "orange")
    plt.legend()
    plt.xlabel("log(N)")
    plt.ylabel("log(err)")
    plt.show()

#define forcing and exact functions as obtained previously
f1 = lambda x : 3 * x - 2
f2 = lambda x: x ** 2 + 3 * x - 2
u1 = lambda u: u **2 - 0.5 * (u ** 3) + 9/6 * u + 1
u2 = lambda u: u **2 - 1/12 * (u ** 4) - 0.5 * (u ** 3) + 45/12 * u + 1


#calculate exact solutions:
u1_ex, u2_ex = generate_exact_functions(u1, u2, plot = False)

#solve the problem for N = 5
err1, err2 = solveP(5, True, u1_ex, u2_ex, f1, f2, u1, u2)

#if required, print the two global error terms:
print(err1, err2)

#perform convergence study
convergence_study()





