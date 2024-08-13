****************************************************************
To compile the abert_erlich_haskell.hs use the follwing command:
1. ghc -O2 -prof -fprof-auto -rtsopts -o aberth_erlich_haskell aberth_erlich_haskell.hs
2. aberth_erlich_haskell.exe +RTS -p
OR - Could simply run
1. ghc -o aberth_erlich_haskell aberth_erlich_haskell.hs
2. aberth_erlich_haskell.exe
**First option is better for later analysis.
-> In this implementation, there's a limition of 850 iterations for demonstration purposes, if needed accurate results simply remove it (best accuarcy can be achieved with epsilon = 1e-4).
-> .hi, .o files are produced while compiling with ghc. .prof (open it with notepad or any editor) is produced when compiling with +RTS -p and is used for analysis.
****************************************************************

****************************************************************
To run aberth_erlich_efficient.py use:
Any IDE that works with python
OR
Run in terminal: python aberth_erlich_efficient.py
-> In this implementation, there's a limition of 850 iterations for demonstration purposes, if needed accurate results simply remove it (best accuarcy can be achieved with epsilon = 1e-4).
-> This implementation uses as many vectorized operations as possible, mostly with numpy library.
****************************************************************

****************************************************************
To run aberth_erlich_less_efficient.py use:
Any IDE that works with python
OR
Run in terminal: python aberth_erlich_less_efficient.py
-> In this implementation, there's a limition of 850 iterations for demonstration purposes, if needed accurate results simply remove it (best accuarcy can be achieved with epsilon = 1e-4).
****************************************************************

Wiki: https://en.wikipedia.org/wiki/Aberth_method


