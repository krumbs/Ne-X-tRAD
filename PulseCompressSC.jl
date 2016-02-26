## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to perform serial Pulse Compression

include("NR_ReadBinFile.jl") 
include("NR_ReadRefSig.jl")
export readRawData
export readRefSig

function SCPC(nPulses_::Int, rawDataFile_::AbstractString, refSigFile_::AbstractString)

## Read in refernce signal from file
dataRaw = readRawData(nPulses_, rawDataFile_)	#"/home/stephaniejonkers/NetRAD/e10_10_21_1927_55_P1_1_130000_S0_1_2047_node3.bin"

## Read the chirp from the txt file
refSig = readRefSig(refSigFile_);		#"/home/stephaniejonkers/NetRAD/refSigN3N3Pl5000.txt"
	
# Array and Variable Initializations
#dataPostiveHalf		= zeros(Complex64, nPulses_, 1024);
dataMatched 		= zeros(Complex64, nPulses_, 2048);
dataIFFT		= zeros(Complex64, nPulses_, 2048);
#x 			= zeros(Complex64, nPulses_, 2048);
# Pulse Compression
P = plan_fft(dataRaw[1,:]);
for n = 1:nPulses_
	dataMatched[n,:] = refSig.*(P*dataRaw[n,:]);
end

# Clear memory
dataRaw		= 0;
refSigRaw 	= 0;
refSig 		= 0;
gc()

P = plan_ifft(dataMatched[1,:]);
for n = 1:nPulses_
	dataIFFT[n,:] = P*dataMatched[n,:];
end

# Clear memory
dataMatched	= 0;
gc()

#=
# Create the analytic signal transform
dataMatched[:,2:2048] 	= 2.*dataMatched[:,2:2048];

# Flatten the negative frequencies
dataPostiveHalf 	= hcat(dataPostiveHalf,dataMatched[:,1025:end]);

# Transform back into time domain
for n = 1:nPulses_
	x[n,:] 		= ifft(dataPostiveHalf[n,:]);
end
=#

# Clear memory
dataMatched		= 0;
#dataPostiveHalf	= 0;
gc()

return dataIFFT
end
tic()
data = SCPC(8000, "/home/stephaniejonkers/NetRAD/e10_10_21_1927_55_P1_1_130000_S0_1_2047_node3.bin", "/home/stephaniejonkers/NetRAD/refSigN3N3Pl5000.txt")
toc()























