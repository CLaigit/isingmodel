# isingmodel

Execute 'make' to obtain the executable file 'andromeda'
  - make with parameters:
    - bls:    block size in CUDA code (default 256)
    - score:  initial conditions using host functions (default using device)
              This parameter can be set to any value

How to access Titan box:
ssh username@ieng6.ucsd.edu
ssh username@igpu6-210.ucsd.edu

How to run the cuda code:
make cuda
./2dising_cuda.sh
