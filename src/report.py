import os
import re
import subprocess

import pandas as pd
import matplotlib.pyplot as plt


MAX_N = 12


def report1():
    time_of_dot = [0] * 4
    time_of_axpby = [0] * 4
    time_of_SpMV = [0] * 4
    time_of_VVbe = [0] * 4

    gflops_of_dot = [0] * 4
    gflops_of_axpby = [0] * 4
    gflops_of_SpMV = [0] * 4
    gflops_of_VVbe = [0] * 4


    for i, (Nx, Ny, K1, K2, maxiter, tol, n, count_calls) in enumerate([
        (100, 100, 1, 0, 10, 0.0, 1, 1000),
        (1000, 100, 1, 0, 10, 0.0, 1, 100),
        (1000, 1000, 1, 0, 10, 0.0, 1, 1),
        (10000, 1000, 1, 0, 10, 0.0, 1, 1),
    ]):
        for _ in range(count_calls):
            result = subprocess.run([
                "./bin/main",
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

            t = re.search(r"Time of dot: (.*?) seconds, gflops: (.*?)\n", stdout)
            time_of_dot[i] += float(t.group(1))
            gflops_of_dot[i] += float(t.group(2))
            t = re.search(r"Time of axpby: (.*?) seconds, gflops: (.*?)\n", stdout)
            time_of_axpby[i] += float(t.group(1))
            gflops_of_axpby[i] += float(t.group(2))
            t = re.search(r"Time of SpMV: (.*?) seconds, gflops: (.*?)\n", stdout)
            time_of_SpMV[i] += float(t.group(1))
            gflops_of_SpMV[i] += float(t.group(2))
            t = re.search(r"Time of VVbe: (.*?) seconds, gflops: (.*?)\n", stdout)
            time_of_VVbe[i] += float(t.group(1))
            gflops_of_VVbe[i] += float(t.group(2))

    df = pd.DataFrame(
        [
            [g/t for t, g in zip(time_of_dot, gflops_of_dot)],
            [g/t for t, g in zip(time_of_axpby, gflops_of_axpby)],
            [g/t for t, g in zip(time_of_SpMV, gflops_of_SpMV)],
            [g/t for t, g in zip(time_of_VVbe, gflops_of_VVbe)],
            [
                g/t for g, t in zip(
                    [sum(args) for args in zip(
                    gflops_of_dot,
                    gflops_of_axpby,
                    gflops_of_SpMV,
                    gflops_of_VVbe,
                    )],
                    [sum(args) for args in zip(
                    time_of_dot,
                    time_of_axpby,
                    time_of_SpMV,
                    time_of_VVbe,
                    )]
                )
            ],
        ],
        index=["dot", "axpby", "SpMV", "VVbe", "solve"],
    ).T

    df["N"] = [10**4, 10**5, 10**6, 10**7]
    
    df.to_excel("report_1.xlsx", index=False)

    plt.plot(df["N"], df['dot'])
    plt.plot(df["N"], df['axpby'])
    plt.plot(df["N"], df['SpMV'])
    plt.plot(df["N"], df['VVbe'])
    plt.plot(df["N"], df['solve'])
    plt.legend(["dot", "axpby", "SpMV", "VVbe", "solve"])
    plt.title("GFLOPS (N)")
    plt.xlabel("N")
    plt.ylabel("GFLOPS")
    plt.xscale("log")
    plt.savefig('report1.jpg')
    plt.close()


def report2():
    time_of_dot = [[0]*MAX_N for i in range(4)]
    time_of_axpby = [[0]*MAX_N for i in range(4)]
    time_of_SpMV = [[0]*MAX_N for i in range(4)]
    time_of_VVbe = [[0]*MAX_N for i in range(4)]
    time_of_solving = [[0]*MAX_N for i in range(4)]

    for n in range(1, MAX_N + 1):
        print(n)
        for i, (Nx, Ny, K1, K2, maxiter, tol, n, count_calls) in enumerate([
            (100, 100, 1, 0, 10, 0.0, n, 1000),
            (1000, 100, 1, 0, 10, 0.0, n, 100),
            (1000, 1000, 1, 0, 10, 0.0, n, 1),
            (10000, 1000, 1, 0, 10, 0.0, n, 1),
        ]):
            for _ in range(count_calls):
                result = subprocess.run([
                    "./bin/main",
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

                time_of_dot[i][n-1] += float(re.search(r"Time of dot: (.*?) seconds", stdout).group(1))
                time_of_axpby[i][n-1] += float(re.search(r"Time of axpby: (.*?) seconds", stdout).group(1))
                time_of_SpMV[i][n-1] += float(re.search(r"Time of SpMV: (.*?) seconds", stdout).group(1))
                time_of_VVbe[i][n-1] += float(re.search(r"Time of VVbe: (.*?) seconds", stdout).group(1))

    for i, n in enumerate([10**4, 10**5, 10**6, 10**7]):
        df = pd.DataFrame(
            [
                [time_of_dot[i][0]/t for t in time_of_dot[i]],
                [time_of_axpby[i][0]/t for t in time_of_axpby[i]],
                [time_of_SpMV[i][0]/t for t in time_of_SpMV[i]],
                [time_of_VVbe[i][0]/t for t in time_of_VVbe[i]],
                [
                    sum([time_of_dot[i][0], time_of_axpby[i][0], time_of_SpMV[i][0], time_of_VVbe[i][0]]) 
                    / sum(args) for args in zip(
                        time_of_dot[i],
                        time_of_axpby[i],
                        time_of_SpMV[i],
                        time_of_VVbe[i],
                    )
                ]
            ],
            index=["dot", "axpby", "SpMV", "VVbe", "solve"],
        ).T

        df["T"] = range(1, MAX_N+1)

        df.to_excel(f"report_{n}.xlsx", index=False)

        plt.plot(df["T"], df['dot'])
        plt.plot(df["T"], df['axpby'])
        plt.plot(df["T"], df['SpMV'])
        plt.plot(df["T"], df['VVbe'])
        plt.plot(df["T"], df['solve'])
        plt.legend(["dot", "axpby", "SpMV", "VVbe", "solve"])
        plt.title(f"Acceleration (T) N={n}")
        plt.xlabel("T")
        plt.ylabel("acceleration")
        plt.savefig(f'report_{n}.jpg')
        plt.close()


if __name__ == "__main__":
    # report1()
    report2()
