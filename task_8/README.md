## Overview
    Approximate system of 2 linear DE

    My first attempt to do something on GPU,
    so I'm authorizing myself to make methods  
    that may require huge attention to realization
    details when applied to general problem

## Status
    In development, WIP

## Contains
    src/methods/rk2.py      -- Runge–Kutta 2nd order with constant step
    src/methods/rk2wc.py    -- Runge–Kutta 2nd order with auto step

## Dependencies
    python      >= 3.12
    cuda        >= 12.6.85
    numpy       >= 2.1
    pycuda      >= 2024.1.2
    nvidia gpu  >= Turing

## Init
    git clone https://github.com/vshulcz/vman
    mv vman/src/vman.sh .
    chmod +x vman.sh
    ./vman.sh
    rm -r vman vman.sh
    ./venv/bin/pip install pycuda 
    ./venv/bin/pip install numpy 

## Start
    ./venv/bin/python main.py
