An nonlinear programming test for [CasADi](https://web.casadi.org/)

Follow this [link](https://github.com/casadi/casadi/wiki/InstallationLinux) to install CasADi and Ipopt.

Only thing need to notice is when build CasADi, instead of using 
```
cmake -DWITH_PYTHON=ON ..
``` 
in the instruction,
the ipopt compatible should be included, like
```
cmake -DWITH_PYTHON=ON -DWITH_IPOPT=ON ..
```