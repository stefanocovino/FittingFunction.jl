module FittingFunction


export PL



"""
    PL(E,N,β;E₀ = 1.)

Computes a power-law with spectral index 'β' and normalization 'N' at input 'E'. You can normalize the power-law at 'E₀'.


# Examples
```jldoctest
PL(3.,1.,-1.)

# output

0.3333333333333333
```
"""
function PL(E,N,β;E₀ = 1.)
    spec = .^(E./E₀,β)
    return N .* spec
end


end
