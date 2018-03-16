#pragma TextEncoding = "MacRoman"		// For details execute DisplayHelpTopic "The TextEncoding Pragma"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "smallwood"

//Similar to FitRayleighBkgd, but includes an option to fit the righ-hand side as well. 
//Example of how to use:
//Run DTfileconvert_single_cls to produce pump-probe arrays.
//Run phased version of FTSA to produce the file "pumpprobe.dat" in the output folder.
//Run the function below, locating the appropriate path for "PumpProbe" and choosing appropriate limits.
//e.g., FitScatteringBkgd(-inf,1455,0.2)
function FitScatteringBkgd(ene_min,ene_max,delay,[E2min,E2max])
variable ene_min, ene_max, delay
variable E2min, E2max
LoadMatLabDat("PumpProbe",OriginalFilename="PumpProbe.dat")
LoadMatLabDat("Ref",OriginalFilename="Ref_Filtered.dat")
wave PumpProbe, Ref
wavestats/Q Ref; Ref -= V_min
ene_min = ene_min < dimoffset(PumpProbe,0) ? dimoffset(PumpProbe,0) : ene_min
ene_max = ene_max > dimfinal(PumpProbe,0) ? dimfinal(PumpProbe,0) : ene_max
E2min = E2min < dimoffset(PumpProbe,0) ? dimoffset(PumpProbe,0) : E2min
E2max = paramisdefault(E2max) ? inf : E2max
E2max = E2max > dimfinal(PumpProbe,0) ? dimfinal(PumpProbe,0) : E2max
if(!paramisdefault(E2min))
	duplicate/O PumpProbe PumpProbe2
	PumpProbe2 = x < ene_min ? nan : ( x < ene_max ? PumpProbe2(x) : ( x < E2min ? nan : ( x < E2max ? PumpProbe2(x) : nan ) ) )
else
	duplicate/O/R=(ene_min,ene_max) PumpProbe PumpProbe2
endif
display PumpProbe PumpProbe2
make/N=4/O W_coef = {0.5,0.5,delay,0.1}
FuncFit/H="0010"/NTHR=0 ScatteringBkgd W_coef PumpProbe2 /D 
wave Fit_PumpProbe2
duplicate/O PumpProbe Fit_PumpProbe2b
Fit_PumpProbe2b = ScatteringBkgd(W_coef,x)
appendtograph Fit_PumpProbe2b
ModifyGraph rgb(PumpProbe)=(0,0,0),rgb(fit_PumpProbe2b)=(0,0,65535)
end

Function ScatteringBkgd(w0,w) : FitFunc
	Wave w0
	Variable w
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(w) = A+B*cos( t0*x - phi0 )
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w0[0] = A
	//CurveFitDialog/ w0[1] = B
	//CurveFitDialog/ w0[2] = t0
	//CurveFitDialog/ w0[3] = phi0
	wave Ref, PumpProbe
	variable ma
	variable planck = 4.13566766 // meV/THz
	wavestats/M=1/Q PumpProbe
	ma = max(abs(V_max),abs(V_min))
	return Ref(w)*(ma*w0[0]+ma*w0[1]*cos( 2*pi/planck*w0[2]*w - pi*w0[3] ))
End
