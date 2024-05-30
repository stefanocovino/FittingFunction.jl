var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = FittingFunction","category":"page"},{"location":"#FittingFunction","page":"Home","title":"FittingFunction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FittingFunction.","category":"page"},{"location":"","page":"Home","title":"Home","text":"A set of functions I am often using in carrying out regressions.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [FittingFunction]","category":"page"},{"location":"#FittingFunction.Band-NTuple{5, Any}","page":"Home","title":"FittingFunction.Band","text":"Band(E,α,β,E₀,A)\n\nCompute the so-called 'Band function', at input energy 'E'. 'α' is the low-energy index, 'β' the high-energy index and 'E₀' the break energy. 'A' is the normalization. The peak-energy is 'Ep = (2+α)E₀'.\n\nExamples\n\nBand(300.,-2.,-1.,200.,1.)\n\n# output\n\n-0.45304697140984085\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.CPL-NTuple{4, Any}","page":"Home","title":"FittingFunction.CPL","text":"CPL(E,N,α,Ec;E0=1.)\n\nComputes a cut-off power-law with spectral index 'α' and normalization 'N' at input 'E' and cut-off energy 'Ec'. You can normalize the power-law at 'E0'.\n\nExamples\n\nCPL(3.,1.,-1.,2.)\n\n# output\n\n0.0743767200494766\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.Counts2Mag-Tuple{Any, Any}","page":"Home","title":"FittingFunction.Counts2Mag","text":"Counts2Mag(cts,ects;zp=25.0,ezp=0.0)\n\nConvert counts (or flux) 'cts', with uncertainty 'ects', to magnitude scale with zero-point 'zp' and relative uncertainty.\n\nExamples\n\nCounts2Mag(100,10)\n\n# output\n\n(20.0, 0.10857362047581294)\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.Extinction-Tuple{Any, Any}","page":"Home","title":"FittingFunction.Extinction","text":"Extinction(wave,EBV;gal=\"SMC\",Rv=FFGals[\"SMC\"],z=0.)\n\nCompute the UV/optical/extinction at the input wavelength 'wave', in Angstrom. 'EBV' is E(B-V) in magnitudes, 'gal' is one of the galaxy extinction recipes listed in the 'FFGals' (exported) dictionary. 'Rv' is the seelctive extinction and 'z' is the redshift of the absorpber. References about the adopted extinction curves are discussed in the documentation of the DustExtinction package.\n\nExamples\n\nExtinction(5500.,1,gal=\"MW\",Rv=FFGals[\"MW\"],z=0.)\n\n# output\n\n0.06017150537641227\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.GaussAbs-NTuple{4, Any}","page":"Home","title":"FittingFunction.GaussAbs","text":"GaussAbs(E,Ed,σ,El)\n\nComputes a Gaussian absorption with depth 'Ed' at energy 'E' with center in 'El' and width (FWHM/2.35) 'σ'. 'A.\n\nExamples\n\nGaussAbs(2.,1.,0.1,2.)\n\n# output\n\n0.018510395162776326\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.Gaussian-NTuple{4, Any}","page":"Home","title":"FittingFunction.Gaussian","text":"Gaussian(E,A,σ,El)\n\nComputes a Gaussian at energy 'E' with center in 'El' and width (FWHM/2.35) 'σ'. 'A is the normalization.\n\nExamples\n\nGaussian(3.,1.,0.1,2.)\n\n# output\n\n4.009419869036418\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.GetAtomicData","page":"Home","title":"FittingFunction.GetAtomicData","text":"GetAtomicData(table::String=\"\")::DataFrame\n\nProvides a table with atomic data among those available. Calling the function with no parameter shows the available tables. At present we have:\n\n\"FitLyman\". This is the original table generated for the \"FitLyman\" ESO-MIDAS package.\n\"JitrikBunge04\". This is a collection of atomic data prepared by Jitrik & Bunge (2004).\n\"VALD03\". This is a collection of atomic data extracted from the VALD3 database.\n\nExamples\n\nfly = GetAtomicData(\"FitLyman\")\n\n\n\n\n\n","category":"function"},{"location":"#FittingFunction.Mag2Counts-Tuple{Any, Any}","page":"Home","title":"FittingFunction.Mag2Counts","text":"Mag2Counts(mag,emag;zp=25.0)\n\nConvert magnitudes 'mag' with uncertainty 'emag' to counts (or flux) with zero-point 'zp'.\n\nExamples\n\nMag2Counts(20,0.1)\n\n# output\n\n(100.0, 9.210340371976184)\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.PL-Tuple{Any, Any, Any}","page":"Home","title":"FittingFunction.PL","text":"PL(E,N,α;E0=1.)\n\nComputes a power-law with spectral index 'α' and normalization 'N' at input 'E'. You can normalize the power-law at 'E0'.\n\nExamples\n\nPL(3.,1.,-1.)\n\n# output\n\n0.3333333333333333\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.Pol2Stokes-NTuple{4, Any}","page":"Home","title":"FittingFunction.Pol2Stokes","text":"Pol2Stokes(i, p, t, c; ep=0.0, et=0.0, ec=0.0)\n\nComputes the Stokes parameters for a given polarisation degree and position angle. 'i' is total intensity, 'p' the polarisation degree, 't' the position angle (randians) and 'c' the circular polarisation. The errors 'ep', 'ew, and 'e' are optional. The ouput is: intensity 'q', 'u', 'v' Stokes parameters and respective errors.\n\nExamples\n\nPol2Stokes(1.,0.1,0.1,0.)\n\n# output\n\n(1.0, 0.09800665778412417, 0.019866933079506124, 0.0, 0.0, 0.0, 0.0)\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.SBPL-NTuple{5, Any}","page":"Home","title":"FittingFunction.SBPL","text":"SBPL(E,N,α,β,Eb)\n\nComputes a smoothly joint broken power-law with spectral index 'α' and 'β', a break at 'Eb', normalization 'N' and input energy'E'.\n\nExamples\n\nSBPL(6.,1.,-1.,-1.5,5.)\n\n# output\n\n0.15214515486254615\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.SBPL2-NTuple{7, Any}","page":"Home","title":"FittingFunction.SBPL2","text":"SBPL2(E,N,α,β,γ,Eb1,Eb2)\n\nComputes a smoothly joint double broken power-law with spectral index 'α', 'β' and 'γ', and two breaks at 'Eb1' and 'Eb2'. The normalization is 'N' and the input energy 'E'.\n\nExamples\n\nSBPL2([4.,6.,8.],1.,-1.,-1.5,-2.,5.,7.)\n\n# output\n\n3-element Vector{Float64}:\n 0.25\n 0.06804138174397717\n 0.04133986423538423\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.SigmaClip","page":"Home","title":"FittingFunction.SigmaClip","text":"SigmaClip(x, ex=ones(size(x)); sigmacutlevel=2)\n\nFilters an input array 'x' with optional uncertainties 'ex' with one-iteration sigma clipping at the level 'sigmacutlevel'. It reports a mask to select the surviving elements in the input arrays or other relwted arrays.\n\nExamples\n\nx = [4.,6.,8.,1.,3.,5.,20.]\nmask = SigmaClip(x)\nx[mask]\n\n# output\n\n6-element Vector{Float64}:\n 4.0\n 6.0\n 8.0\n 1.0\n 3.0\n 5.0\n\n\n\n\n\n","category":"function"},{"location":"#FittingFunction.Stokes2Pol-NTuple{4, Any}","page":"Home","title":"FittingFunction.Stokes2Pol","text":"Stokes2Pol(i, q, u, v; eq=0.0, eu=0.0, ev=0.0)\n\nComputes the polarisation degree and position angle given the Stokes parameters: 'i', total intensity, and 'q', 'u' and 'v'. The errors 'eq', 'eu, and 'ev' are optional. The ouput is: intensity 'i', polarisation degree, 'p', position angle (randians) 'theta' and circular polarisation 'chi', with respective errors.\n\nExamples\n\nStokes2Pol(1.,0.1,0.1,0.)\n\n# output\n\n(1.0, 0.14142135623730953, 0.39269908169872414, 0.0, 0.0, 0.0, 0.0)\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.TauVoigt-NTuple{5, Any}","page":"Home","title":"FittingFunction.TauVoigt","text":"TauVoigt(λ,NHI,DopplerBroadening,z,transition)\n\nComputes the opacity due to a given absorption transition. lambda are the wavelenghths (cm), N the column density (cm^-2), 'DopplerBroadening' is expressed in (cms^-1), 'z' is the source redshift and 'transition' is a tuple formed by the central wawelength (cm), the oscillator strength and the damping coefficient (s^-1) for the given transition. The Doppler broadening factor is typicaly due to either an intrinsic thermal broadening (b_K = sqrt2KTm) or to a turbulent motion of the gas (b_T = sqrtsigma_T, where sigma_T is the inner velocity dispersion). These two factors can be summed in quadrature, i.e. b = sqrtb_K^2 + b_T^2.\n\nExamples\n\nTauVoigt(range(start=1494.5,stop=1495.5,step=0.1) .* 1e-8,1e18,1e5,0.,[1495.05*1e-8,0.54,8.106e8])\n\n# output\n\n11-element Vector{Float64}:\n  0.5408331361420329\n  0.8079604230648749\n  1.3357686994197793\n  2.618871809580839\n  7.282358478732964\n 66.40911467182524\n 66.40911467182524\n  7.282358478732964\n  2.618871809580839\n  1.3357686994184856\n  0.8079604230648749\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.VoigtFunctTG-Tuple{Any, Any}","page":"Home","title":"FittingFunction.VoigtFunctTG","text":"VoigtFuncTG(a,v)\n\nApproximates the Voigt (or line broadening function) function. It was introduced by Tepper-Garcia (2006). It is of high accuracy provided that a le 10^-4. The a and v parameters are defined as: a = fracGamma lambda_c4pi b 10^13 and v = frac(lambda_c - lambda)cb lambda_c sqrt2 log 2, where lambda_c is the central wavelength of the transition, c the speed of light, b the Doppler broadening factor and Gamma the damping coefficient.\n\nExamples\n\nFittingFunction.VoigtFunctTG(1e-5,1.)\n\n# output\n\n0.36788094737611143\n\n\n\n\n\n","category":"method"},{"location":"#FittingFunction.XAbs-Tuple{Any}","page":"Home","title":"FittingFunction.XAbs","text":"XAbs(E; NH=1e20, z=0)\n\nComputes the effective absorption cross section per hydrogen atom following Morrison & McCammon (1983). 'E' is the energy (KeV, in the 0.03-10 KeV range), 'NH' is the column density of hydrogen atom (particle per square cm) and 'z' is the redshift.\n\nExamples\n\nXAbs([0.5,1.25,2.], NH=1e20, z=0)\n\n# output\n\n3-element Vector{Float64}:\n 0.9290803992951834\n 0.9868927382232576\n 0.9957079871273042\n\n\n\n\n\n","category":"method"}]
}