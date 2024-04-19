module FittingFunction


export PL
export SBPL



"""
    PL(E,N,B;E0=1.)

Computes a power-law with spectral index 'β' and normalization 'N' at input 'E'. You can normalize the power-law at 'E₀'.


# Examples
```jldoctest
PL(3.,1.,-1.)

# output

0.3333333333333333
```
"""
function PL(E,N,B; E0=1.)
    spec = .^(E./E0,B)
    return N .* spec
end





"""
    SBPL(E,N,A1,A2,Eb)

Computes a smoothly joint broken power-law with spectral index 'A1' and 'A2', a break at 'Eb' and normalization 'N' at input 'E'.


# Examples
```jldoctest
SBPL(6.,1.,-1.,-1.5,5.)

# output

0.15214515486254615
```
"""
function SBPL(E,N,A1,A2,Eb)
    f1 = PL(E,N,A1)
    n = PL(Eb,N,A1) ./ PL(Eb,1.0,A2)
    f2 = PL(E,n,A2)
    return ifelse.(E .<= Eb, f1, f2)
end



end
