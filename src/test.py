import os
import subprocess

def main():
    print("Test correct parrallel")

    for Nx, Ny, K1, K2, maxiter, tol in [
        (1, 1, 1, 1, 10, 0),
        (5, 5, 0, 10, 100, 0),
        (5, 5, 10, 0, 100, 0),
        (10, 10, 1, 1, 10, 0),
        (10, 10, 2, 2, 1000, 0.001),
    ]:
        for n in range(1, 8):
            result = subprocess.run([
                "./bin/main.exe",
                "--Nx", str(Nx),
                "--Ny", str(Ny),
                "--K1", str(K1),
                "--K2", str(K2),
                "--maxiter", str(maxiter),
                "--tol", str(tol),
                "-n", str(n),
                "--filename", "tmp"],
                capture_output=True,
            )

            assert result.returncode == 0

            if n == 1:
                with open("tmp") as f:
                    s = f.read()
            else:
                with open("tmp") as f:
                    assert s == f.read()

    os.remove("tmp")


if __name__ == "__main__":
    main()
