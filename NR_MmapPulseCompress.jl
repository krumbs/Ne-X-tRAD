## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to use memory mapping for pulse compression

function createDataMmap(N_::Int, fileName_::AbstractString)
local RadarData

# Open the bin file
binStream = open(fileName_, "r");

# Create a Matrix whose values are linked to a file saved on disk rather than memory in cache
radarData = Mmap.mmap(binStream, Array{UInt16,1}, 130000*2048);

# Remove bias introduced by ADC
rawData = radarData - mean(radarData[1:2048]);
return rawData 
end
#-----------------------------------------------------------------------------------------
using DSP

function readRefSig(file_name_)
local refSigRaw_
refSigRaw_  	= readcsv(file_name_, Float32);	

## Process resig
refSigRaw_ 	= refSigRaw_/sqrt(mean(abs(refSigRaw_).^2));
refSigRaw_	= (hanning(size(refSigRaw_,2)))'.*refSigRaw_; 
refSig_	 	= [refSigRaw_ zeros(1,2048-size(refSigRaw_,2))];
refSig_	 	= fft(refSig_);
refSig_		= conj(refSig_);
return refSig_
end
#-----------------------------------------------------------------------------------------
## Read in refernce signal from file 
dataRaw = createDataMmap(1, "/home/stephanie/Desktop/e10_10_21_1927_55_P1_1_130000_S0_1_2047_node3.bin"); 

## Read the chirp from the txt file
refSig = readRefSig("/home/stephanie/Desktop/refSigN3N3Pl5000.txt");

function SCPC(nPulses_, radarData_, refSig_)
local w
local P
local specData = zeros(2048*nPulses_,1)*im;

# Transform to frequency domain and convolve
P = plan_fft(radarData_[1:2048]);
for w = 1:2048:2048*nPulses_
	specData[w:w+2047] = refSig_[1:2048].*(P*radarData_[w:w+2047]);
end

# Create a file for Mmapping
outFile = open("/home/stephanie/Desktop/02-03/SCPC.bin", "w+");
write(outFile, nPulses_);

# Transform back to time domain and write to file
Plan = plan_ifft(specData[1:2048]);
for w = 1:2048:2048*nPulses_
	write(outFile, P*specData[w:w+2047]);
end
close(outFile);
end
gc()
