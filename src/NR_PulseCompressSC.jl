## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to perform serial Pulse Compression

include("NR_ReadBinFile.jl") 
include("NR_ReadRefSig.jl")
export readRawData
export readRefSig

function SCPC(nPulses_::Int, rawDataFile_::AbstractString, refSigFile_::AbstractString, analytic_::Bool=false)

## Read in refernce signal from file
d = readRawData(nPulses_, rawDataFile_)	# Function that processes the radar data

## Read the chirp from the txt file
r = readRefSig(refSigFile_);		# Function that processes the reference signal 
	
## Array and Variable Initializations
dPostiveHalf = zeros(Complex64, nPulses_, 1024);
dMatched = zeros(Complex64, nPulses_, 2048);
dIFFT = zeros(Complex64, nPulses_, 2048);
dOut = zeros(Complex64, nPulses_, 2048);

## Pulse Compression
P = plan_fft(d[1,:]);
for n = 1:nPulses_
	dMatched[n,:] = r.*(P*d[n,:]);
end

## Clear memory
d = 0;
r = 0;
gc()

P = plan_ifft(dMatched[1,:]);

## Branch for analytic signal
if analytic_==false
	for n = 1:nPulses_
		dIFFT[n,:] = P*dMatched[n,:];
	end

	## Clear memory
	dMatched = 0;
	gc()
	return dIFFT
else
	## Create the analytic signal transform
	dMatched[:,2:2048] = 2.*dMatched[:,2:2048];

	## Flatten the negative frequencies
	dPostiveHalf = hcat(dPostiveHalf,dMatched[:,1025:end]);

	## Transform back into time domain
	for n = 1:nPulses_
		dOut[n,:] = ifft(dPostiveHalf[n,:]);
	end

	## Clear memory
	dMatched = 0;
	dPositiveHalf = 0;
	gc()
	return dOut
end

end
tic()
data = SCPC(10000, "/home/jnkste003/thesis_code/Ne-X-tRAD/data/e10_10_21_1927_55_P1_1_130000_S0_1_2047_node3.bin", "/home/jnkste003/thesis_code/Ne-X-tRAD/data/refSigN3N3Pl5000.txt", true)
toc()


