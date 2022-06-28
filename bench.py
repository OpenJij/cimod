import dimod
import cimod
import time

fil = open("benchmark", "w")
fil.write("N t_dimod t_cimod\n")

def benchmark(N, test_fw):
    linear = {}
    quadratic = {}

    spin = {}

    # interactions

    for i in range(N):
        spin[i] = 1

    for elem in range(N):
        linear[elem] = 2.0*elem;

    for i in range(N):
        for j in range(i+1, N):
            if i != j:
                quadratic[(i,j)] = (i+j)/(N)

    t1 = time.time()

    # initialize
    a = test_fw.BinaryQuadraticModel(linear, quadratic, 0, test_fw.BINARY)
    a.change_vartype(test_fw.SPIN)

    # calculate energy for 50 times.
    for _ in range(50):
        print(a.energy(spin))

    t2 = time.time()

    return t2-t1

d_arr = []
c_arr = []

for N in [25, 50, 100, 200, 300, 400, 600, 800,1000, 1600, 2000, 3200, 5000]:
    print("N {}".format(N))
    d = benchmark(N, dimod)
    c = benchmark(N, cimod)
    print("{} {} {}".format(N, d, c))
    fil.write("{} {} {}\n".format(N, d, c))
