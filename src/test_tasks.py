import re
import subprocess

MAX_PROC = 8


def pairs_of_factors(val):
    return [
        *[(int(i), int(val / i)) for i in range(1, int(val**0.5)+1) if val % i == 0],
        *[(int(val / i), int(i)) for i in range(1, int(val**0.5)+1) if val % i == 0],
    ]


def main():
    print("Test correct parrallel")

    for Nx, Ny, K1, K2, maxiter, tol in [
        (1, 1, 1, 1, 10, 0),
        (5, 5, 0, 10, 100, 0),
        (5, 5, 10, 0, 100, 0),
        (10, 10, 1, 1, 10, 0),
        (10, 10, 2, 2, 1000, 0.001),
    ]:
        results = []

        for n in range(1, MAX_PROC+1):
            result = subprocess.run([
                "./bin/task1",
                "--Nx", str(Nx),
                "--Ny", str(Ny),
                "--K1", str(K1),
                "--K2", str(K2),
                "--maxiter", str(maxiter),
                "--tol", str(tol),
                "-n", str(n),
                ], capture_output=True,
            )
            print(' '.join(result.args))
            assert result.returncode == 0
            results.append(re.findall(r"Number of iteration: .*?\n", result.stdout.decode()))

        for n in range(1, MAX_PROC+1):
            for Px, Py in pairs_of_factors(n):
                print(f"mpirun -np {n} ./bin/task2 --Nx {Nx} --Ny {Ny} --K1 {K1} --K2 {K2} --Px {Px} --Py {Py} --maxiter {maxiter} --tol {tol}")
                result = subprocess.run(
                    [f"mpirun -np {n} ./bin/task2 --Nx {Nx} --Ny {Ny} --K1 {K1} --K2 {K2} --Px {Px} --Py {Py} --maxiter {maxiter} --tol {tol}"],
                    capture_output=True,
                    shell=True,
                )

                assert result.returncode == 0
                results.append(re.findall(r"Number of iteration: .*?\n", result.stdout.decode()))

        for row1, row2 in zip(results[:-1], results[1:]):
            assert set(row1) == set(row2)


if __name__ == "__main__":
    main()
