var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = FittingFunction","category":"page"},{"location":"#FittingFunction","page":"Home","title":"FittingFunction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FittingFunction.","category":"page"},{"location":"","page":"Home","title":"Home","text":"A set of functions I am often using in carrying out regressions.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [FittingFunction]","category":"page"},{"location":"#FittingFunction.Band-NTuple{5, Any}","page":"Home","title":"FittingFunction.Band","text":"Band(E,α,β,E₀,A)\n\nCompute the so-called 'Band function', at input energy 'E'. 'α' is the low-energy index, 'β' the high-energy index and 'E₀' the break energy. 'A' is the normalization. The peak-energy is 'Ep = (2+α)E₀'.\n\nExamples\n\nBand(300.,-2.,-1.,200.,1.)\n\n# output\n\n-0.45304697140984085\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.CPL-NTuple{4, Any}","page":"Home","title":"FittingFunction.CPL","text":"CPL(E,N,α,Ec;E0=1.)\n\nComputes a power-law with spectral index 'α' and normalization 'N' at input 'E'. You can normalize the power-law at 'E₀'.\n\nExamples\n\nCPL(3.,1.,-1.,2.)\n\n# output\n\n0.0743767200494766\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.Extinction-Tuple{Any, Any}","page":"Home","title":"FittingFunction.Extinction","text":"Extinction(wave,EBV;gal=\"SMC\",Rv=FFGals[\"SMC\"],z=0.)\n\nCompute the UV/optical/extinction at the input wavelength 'wave', in Angstrom. 'EBV' is E(B-V) in magnitudes, 'gal' is one of the galaxy extinction recipes listed in the 'FFGals' (exported) dictionary. 'Rv' is the seelctive extinction and 'z' is the redshift of the absorpber. References about the adopted extinction curves are discussed in the documentation of the DustExtinction package.\n\nExamples\n\nExtinction(5500.,1,gal=\"MW\",Rv=FFGals[\"MW\"],z=0.)\n\n# output\n\n0.06017150537641227\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.PL-Tuple{Any, Any, Any}","page":"Home","title":"FittingFunction.PL","text":"PL(E,N,α;E0=1.)\n\nComputes a power-law with spectral index 'α' and normalization 'N' at input 'E'. You can normalize the power-law at 'E₀'.\n\nExamples\n\nPL(3.,1.,-1.)\n\n# output\n\n0.3333333333333333\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.SBPL-NTuple{5, Any}","page":"Home","title":"FittingFunction.SBPL","text":"SBPL(E,N,α,β,Eb)\n\nComputes a smoothly joint broken power-law with spectral index 'α' and 'β', a break at 'Eb', normalization 'N' and input energy'E'.\n\nExamples\n\nSBPL(6.,1.,-1.,-1.5,5.)\n\n# output\n\n0.15214515486254615\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.SBPL2-NTuple{7, Any}","page":"Home","title":"FittingFunction.SBPL2","text":"SBPL2(E,N,α,β,γ,Eb1,Eb2)\n\nComputes a smoothly joint double broken power-law with spectral index 'α', 'β' and 'γ', and two breaks at 'Eb1' and 'Eb2'. The normalization is 'N' and the input energy 'E'.\n\nExamples\n\nSBPL2([4.,6.,8.],1.,-1.,-1.5,-2.,5.,7.)\n\n# output\n\n3-element Vector{Float64}:\n 0.25\n 0.06804138174397717\n 0.04133986423538423\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.XAbs-Tuple{Any}","page":"Home","title":"FittingFunction.XAbs","text":"XAbs(E; NH=1e20, z=0)\n\nComputes the effective absorption cross section per hydrogen atom following Morrison & McCammon (1983). 'E' is the energy (KeV, in the 0.03-10 KeV range), 'NH' is the column density of hydrogen atom (particle per square cm) and 'z' is the redshift.\n\nExamples\n\nXAbs([0.5,1.25,2.], NH=1e20, z=0)\n\n# output\n\n3-element Vector{Float64}:\n 0.9290803992951834\n 0.9868927382232576\n 0.9957079871273042\n\n\n\n\n\n","category":"method"}]
}
