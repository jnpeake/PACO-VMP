import subprocess as sp

VMRMin = "1"
VMCMin = "1"
VMRMax = "8"
VMCMax = "4"
serverCPU = "16"
serverRAM = "32"
numArray = ["100","200","300","400","500","1000"]


for i in range(0,100):
    for n in numArray:
        filename = "VMP_B" + str(int(n)+i)
        args = ["./VMPProblemGen", n, VMCMin, VMRMin, VMCMax, VMRMax, serverCPU, serverRAM, filename]
        sp.run(args)
