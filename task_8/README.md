## Overview
    Approximate system of 2 linear DE
    My first attempt to do something on GPU,
    so I'm authorizing myself to build methods
    that may require huge attention to realization
    details when applied to general problem

## Status
	Pre-rel

## Contains
    src/methods/rk2.py      -- Runge–Kutta 2nd order with constant step
    src/methods/rk2wc.py    -- Runge–Kutta 2nd order with auto step
	src/methods/rk4.py      -- Runge–Kutta 4th order with constant step
    src/methods/rk4wc.py    -- Runge–Kutta 4th order with auto step
	src/methods/aux         -- Cauchy problem, analytical solution
	src/methods/kernels     -- CUDA kernels
	src/main.py             -- Base logic

## Dependencies
    python      >= 3.12
    cuda        >= 12.6.85
    numpy       >= 2.1
	numba       >= 0.60.0
    pycuda      >= 2024.1.2 (optional) [will be available soon]
	nvidia gpu  >= Turing (optional)   [will be available soon]

## Init
    git clone https://github.com/vshulcz/vman
    mv vman/src/vman.sh .
    chmod +x vman.sh
    ./vman.sh
    rm -r vman vman.sh
    ./venv/bin/pip install numpy numba
	./venv/bin/pip install pycuda (optional)

## Start
    ./venv/bin/python main.py
