import os
import re
import sys
import subprocess

import pandas as pd
import matplotlib.pyplot as plt


MAX_N = 12


def min_pair_of_factors(val):
    pairs = [(int(i), int(val / i)) for i in range(1, int(val**0.5)+1) if val % i == 0]
    min_distance = min([abs(x-y) for x, y in pairs])
    return [(x, y) for x, y in pairs if abs(x-y) == min_distance][0]


class GflopsResult:
    def __init__(self, *, dot=0, axpby=0, SpMV=0, VVbe=0, solver=0):
        self.dot = dot
        self.axpby = axpby
        self.SpMV = SpMV
        self.VVbe = VVbe
        self.solver = solver

    @classmethod
    def parse(cls, stdout):
        return cls(
            dot=float(re.search(r"Time of dot: .*? seconds, gflops: (.*?)\n", stdout).group(1)),
            axpby=float(re.search(r"Time of axpby: .*? seconds, gflops: (.*?)\n", stdout).group(1)),
            SpMV=float(re.search(r"Time of SpMV: .*? seconds, gflops: (.*?)\n", stdout).group(1)),
            VVbe=float(re.search(r"Time of VVbe: .*? seconds, gflops: (.*?)\n", stdout).group(1)),
            solver=float(re.search(r"Time of solver: .*? seconds, gflops: (.*?)\n", stdout).group(1)),
        )

    def __add__(self, other):
        return self.__class__(
            dot=self.dot+other.dot,
            axpby=self.axpby+other.axpby,
            SpMV=self.SpMV+other.SpMV,
            VVbe=self.VVbe+other.VVbe,
            solver=self.solver+other.solver,
        )

    def __truediv__(self, other):
        return self.__class__(
            dot=self.dot/other,
            axpby=self.axpby/other,
            SpMV=self.SpMV/other,
            VVbe=self.VVbe/other,
            solver=self.solver/other,
        )

    def to_dict(self):
        return {
            "dot": self.dot,
            "axpby": self.axpby,
            "SpMV": self.SpMV,
            "VVbe": self.VVbe,
            "solver": self.solver,
        }


def report1():
    gflops = []

    for i, (Nx, Ny, K1, K2, maxiter, tol, n, count_calls) in enumerate([
        (100, 100, 1, 0, 10, 0.0, 1, 1000),
        (1000, 100, 1, 0, 10, 0.0, 1, 100),
        (1000, 1000, 1, 0, 10, 0.0, 1, 2),
        (10000, 1000, 1, 0, 10, 0.0, 1, 2),
    ]):
        gflops.append(GflopsResult())
        for _ in range(count_calls):
            result = subprocess.run([
                "./bin/task1",
                "--Nx", str(Nx),
                "--Ny", str(Ny),
                "--K1", str(K1),
                "--K2", str(K2),
                "--maxiter", str(maxiter),
                "--tol", str(tol),
                "-n", str(n),
            ], capture_output=True)

            assert result.returncode == 0

            stdout = result.stdout.decode()

            gflops[-1] += GflopsResult.parse(stdout)

        gflops[-1] /= count_calls


    df = (
        pd.DataFrame([g.to_dict() for g in gflops])
        .assign(N=[10**4, 10**5, 10**6, 10**7])
    )

    df.to_excel("./reports/report_1.xlsx", index=False)

    plt.plot(df["N"], df['dot'])
    plt.plot(df["N"], df['axpby'])
    plt.plot(df["N"], df['SpMV'])
    plt.plot(df["N"], df['VVbe'])
    plt.plot(df["N"], df['solver'])
    plt.legend(["dot", "axpby", "SpMV", "VVbe", "solver"])
    plt.title("GFLOPS (N)")
    plt.xlabel("N")
    plt.ylabel("GFLOPS")
    plt.xscale("log")
    plt.savefig('./reports/report1.jpg')
    plt.close()


def report2():
    gflops = []

    for n in range(1, MAX_N + 1):
        print(n)
        gflops.append([])
        for i, (Nx, Ny, K1, K2, maxiter, tol, n, count_calls) in enumerate([
            (100, 100, 1, 0, 10, 0.0, n, 1000),
            (1000, 100, 1, 0, 10, 0.0, n, 100),
            (1000, 1000, 1, 0, 10, 0.0, n, 2),
            (10000, 1000, 1, 0, 10, 0.0, n, 2),
        ]):
            gflops[-1].append(GflopsResult())
            for _ in range(count_calls):
                result = subprocess.run([
                    "./bin/task1",
                    "--Nx", str(Nx),
                    "--Ny", str(Ny),
                    "--K1", str(K1),
                    "--K2", str(K2),
                    "--maxiter", str(maxiter),
                    "--tol", str(tol),
                    "-n", str(n),
                ], capture_output=True)

                assert result.returncode == 0

                stdout = result.stdout.decode()

                gflops[-1][-1] += GflopsResult.parse(stdout)

            gflops[-1][-1] /= count_calls


    for i, n in enumerate([10**4, 10**5, 10**6, 10**7]):
        df = (
            pd.DataFrame([row[i].to_dict() for row in gflops])
            .assign(T=range(1, MAX_N+1))
        )

        df.to_excel(f"./reports/report_2_{n}.xlsx", index=False)

        plt.plot(df["T"], df['dot'])
        plt.plot(df["T"], df['axpby'])
        plt.plot(df["T"], df['SpMV'])
        plt.plot(df["T"], df['VVbe'])
        plt.plot(df["T"], df['solver'])
        plt.legend(["dot", "axpby", "SpMV", "VVbe", "solver"])
        plt.title(f"GFLOPS (T) N={n}")
        plt.xlabel("T")
        plt.ylabel("GFLOPS")
        plt.savefig(f'./reports/report_2_gflops_{n}.jpg')
        plt.close()

        plt.plot(df["T"], df['dot']/df.loc[0, 'dot'])
        plt.plot(df["T"], df['axpby']/df.loc[0, 'axpby'])
        plt.plot(df["T"], df['SpMV']/df.loc[0, 'SpMV'])
        plt.plot(df["T"], df['VVbe']/df.loc[0, 'VVbe'])
        plt.plot(df["T"], df['solver']/df.loc[0, 'solver'])
        plt.legend(["dot", "axpby", "SpMV", "VVbe", "solver"])
        plt.title(f"Acceleration (T) N={n}")
        plt.xlabel("T")
        plt.ylabel("acceleration")
        plt.savefig(f'./reports/report_2_acceleration_{n}.jpg')
        plt.close()

        plt.plot(df["T"], df['dot']/df.loc[0, 'dot']/df["T"])
        plt.plot(df["T"], df['axpby']/df.loc[0, 'axpby']/df["T"])
        plt.plot(df["T"], df['SpMV']/df.loc[0, 'SpMV']/df["T"])
        plt.plot(df["T"], df['VVbe']/df.loc[0, 'VVbe']/df["T"])
        plt.plot(df["T"], df['solver']/df.loc[0, 'solver']/df["T"])
        plt.legend(["dot", "axpby", "SpMV", "VVbe", "solver"])
        plt.title(f"Efficiency (T) N={n}")
        plt.xlabel("T")
        plt.ylabel("efficiency")
        plt.savefig(f'./reports/report_2_efficiency_{n}.jpg')
        plt.close()


def report3():
    gflops = []

    for n in range(1, MAX_N + 1):
        print(n)
        gflops.append([])
        for i, (Nx, Ny, K1, K2, maxiter, tol, n, count_calls) in enumerate([
            (100, 100, 1, 0, 10, 0.0, n, 10),
            (1000, 100, 1, 0, 10, 0.0, n, 10),
            (1000, 1000, 1, 0, 10, 0.0, n, 2),
            (10000, 1000, 1, 0, 10, 0.0, n, 2),
        ]):
            Px, Py = min_pair_of_factors(n)

            gflops[-1].append(GflopsResult())
            for _ in range(count_calls):
                # print(_, f"mpirun -np {n} ./bin/task2 --Nx {Nx} --Ny {Ny} --K1 {K1} --K2 {K2} --Px {Px} --Py {Py} --maxiter {maxiter} --tol {tol}")
                result = subprocess.run(
                    [f"mpirun -np {n} ./bin/task2 --Nx {Nx} --Ny {Ny} --K1 {K1} --K2 {K2} --Px {Px} --Py {Py} --maxiter {maxiter} --tol {tol}"],
                    capture_output=True,
                    shell=True,
                )

                assert result.returncode == 0

                stdout = result.stdout.decode()

                gflops[-1][-1] += GflopsResult.parse(stdout)

            gflops[-1][-1] /= count_calls


    for i, n in enumerate([10**4, 10**5, 10**6, 10**7]):
        df = (
            pd.DataFrame([row[i].to_dict() for row in gflops])
            .assign(T=range(1, MAX_N+1))
        )

        df.to_excel(f"./reports/report_3_{n}.xlsx", index=False)

        plt.plot(df["T"], df['dot'])
        plt.plot(df["T"], df['axpby'])
        plt.plot(df["T"], df['SpMV'])
        plt.plot(df["T"], df['VVbe'])
        plt.plot(df["T"], df['solver'])
        plt.legend(["dot", "axpby", "SpMV", "VVbe", "solver"])
        plt.title(f"GFLOPS (T) N={n}")
        plt.xlabel("T")
        plt.ylabel("GFLOPS")
        plt.savefig(f'./reports/report_3_gflops_{n}.jpg')
        plt.close()

        plt.plot(df["T"], df['dot']/df.loc[0, 'dot'])
        plt.plot(df["T"], df['axpby']/df.loc[0, 'axpby'])
        plt.plot(df["T"], df['SpMV']/df.loc[0, 'SpMV'])
        plt.plot(df["T"], df['VVbe']/df.loc[0, 'VVbe'])
        plt.plot(df["T"], df['solver']/df.loc[0, 'solver'])
        plt.legend(["dot", "axpby", "SpMV", "VVbe", "solver"])
        plt.title(f"Acceleration (T) N={n}")
        plt.xlabel("T")
        plt.ylabel("acceleration")
        plt.savefig(f'./reports/report_3_acceleration_{n}.jpg')
        plt.close()

        plt.plot(df["T"], df['dot']/df.loc[0, 'dot']/df["T"])
        plt.plot(df["T"], df['axpby']/df.loc[0, 'axpby']/df["T"])
        plt.plot(df["T"], df['SpMV']/df.loc[0, 'SpMV']/df["T"])
        plt.plot(df["T"], df['VVbe']/df.loc[0, 'VVbe']/df["T"])
        plt.plot(df["T"], df['solver']/df.loc[0, 'solver']/df["T"])
        plt.legend(["dot", "axpby", "SpMV", "VVbe", "solver"])
        plt.title(f"Efficiency (T) N={n}")
        plt.xlabel("T")
        plt.ylabel("efficiency")
        plt.savefig(f'./reports/report_3_efficiency_{n}.jpg')
        plt.close()


if __name__ == "__main__":
    report1()
    report2()
    # report3()
