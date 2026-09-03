"""
    startRMTool()

Prints a short welcome message and returns `wishartpol(c)` (the Marčenko–Pastur
polynomial with symbolic parameter `c`), mirroring the MATLAB start-up script
`startRMTool.m` (`syms m c z r g s y myu; b = wishartpol(c)`).

The symbolic variables are available as `m`, `z`, `c` and through [`RMV`](@ref).
"""
function startRMTool()
    println("RMTool - A Random Matrix and Free Probability Calculator (Julia port)")
    println("Symbolic variables: m, z, c exported; all variables in RMV: ", join(VARNAMES, ", "))
    b = wishartpol(c)
    println("b = wishartpol(c) = ", b)
    b
end
