import sys, time
sys.path.append("../src")
from kie import KIE_Calculation

nreps = 10
start = time.time()

for _ in range(nreps):
    calc = KIE_Calculation(
        "kie.config",
        "SM.out",
        "grid_230_270.out",
        style="g09"
    )

end = time.time()
elapsed = end - start

print(f"\nelapsed: {elapsed:.3f} s")
print(f"per run: {elapsed/nreps*1000:.3f} ms")
