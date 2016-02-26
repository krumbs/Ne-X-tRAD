## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to use memory mapping for pulse compression

function createDataMmap(N_::Int, fileName_::AbstractString)
local radarData
local binStream

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
function readRefSig(fileName_)
local refSigRaw_
refSigRaw_  	= readcsv(fileName_, Float32);	

## Process resig
refSigRaw_ 	= refSigRaw_/sqrt(mean(abs(refSigRaw_).^2));
refSigRaw_	= (hanning(size(refSigRaw_,2)))'.*refSigRaw_; 
refSig_	 	= [refSigRaw_ zeros(1,2048-size(refSigRaw_,2))];
refSig_	 	= fft(refSig_);
refSig_		= conj(refSig_);
return refSig_
end 
#-----------------------------------------------------------------------------------------
function SCPC(nPulses_::Int, radarData_, refSig_)

local w
local P
local specData = zeros(Complex128, 2048,1);

# Create a file for Mmapping
outFile = open("PC.bin", "w+");

# Plan FFTW
PlanRtoC = plan_fft(radarData_[1:2048]);
PlanCtoC = plan_ifft(specData[1:2048]);

# Clear Memory
specData = 0; 

for w = 1:2048:2048*nPulses_
# Transform to frequency domain and convolve then transform back to time domain and write to file
	write(outFile, PlanCtoC*(refSig_[1:2048].*(PlanRtoC*radarData_[w:w+2047])));
end
close(outFile);
end
gc()
#----------------------------------------------------------------------------------------------
@everywhere function SCPC1(nPulses_::Int, radarData_, refSig_)

local w
local P
local specData = zeros(Complex128, 2048*nPulses_,1);

# Plan FFTW
PlanRtoC = plan_fft(radarData_[1:2048]);
PlanCtoC = plan_ifft(specData[1:2048]);

for w = 1:2048:2048*nPulses_
# Transform to frequency domain and convolve then transform back to time domain
	specData[w:w+2047] = PlanCtoC*(refSig_[1:2048].*(PlanRtoC*radarData_[w:w+2047]));
end
return specData
end
#----------------------------------------------------------------------------------------------
## Altered distribute "DistributedArrays.jl" distribute() function to distribute the raw data matrix into specified distribution
using DistributedArrays
function mydistribute(A::AbstractArray, mydist::AbstractArray; myprocs=workers()[1:min(nworkers(), maximum(size(A)))])
    owner = myid()
    rr = RemoteRef()
    put!(rr, A)
    d = DArray(size(A), myprocs, mydist) do I
        remotecall_fetch(()->fetch(rr)[I...], owner)
    end
    # Ensure that all workers have fetched their localparts.
    for chunk in d.chunks
        wait(chunk)
    end
    return d
end
#-----------------------------------------------------------------------------------------------
## Parallel Pulse Compression function
function MCPC(nPulses_::Int, radarDataFileName_::AbstractString, refSigFileName_::AbstractString)

## Array and Variable Initializations 
local const P = nprocs()-1;  
local const PCpart = convert(Int, floor(nPulses_/P));

local Ddata = zeros(Complex128, 2048*nPulses_);
local radarData
local binStream1
local DdataOut
local DdataRefs

## Read in raw data from file 
radarData = createDataMmap(nPulses_, radarDataFileName_); 

## Read the chirp from the txt file
refSig = readRefSig(refSigFileName_); #"/home/stephanie/Desktop/refSigN3N3Pl5000.txt"

# Distribute the raw data split in the row dimension
Ddata = mydistribute(radarData[1:2048*nPulses_], [P,]);

# Move all mycorrelate() to all Ddata local to each worker and store it to a list, note that @spawnat output is a remote reference
Drefs = [(@spawnat (procs(Ddata))[w] SCPC1(PCpart, localpart(Ddata), refSig)) for w=1:P];

# Clear distr. data memory
#Ddata = 0;
#gc()	

# Fetch all the remote references
DdataArray = Array{Any,1};
DdataArray = pmap(fetch, Drefs);


# Move the array of processed chunks back into local memory
DdataArray = convert(Array, DdataArray);

DdataRefs = DdataArray[1];
for w = 2:P
DdataOut = vcat(DdataRefs, DdataArray[w]);
DdataRefs = DdataOut;
end

# Clear memory
DdataOut = 0;
gc()

# Write to file
binStream1 = open("PC1.bin", "w");
write(binStream1, DdataRefs);

end

tic()
MCPC(130000, "/home/stephanie/Desktop/e10_10_21_1927_55_P1_1_130000_S0_1_2047_node3.bin", "/home/stephanie/Desktop/refSigN3N3Pl5000.txt")
toc()
