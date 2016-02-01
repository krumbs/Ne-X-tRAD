## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to perform parallel Pulse Compression on NetRAD data

Pkg.add("DistributedArrays")
using DistributedArrays

include("NR_ReadBinFile.jl") 
include("NR_ReadRefSig.jl")
export readRawData
export readRefSig

#---------------------------------------------------------------------------------------------------------------------------------
## Altered distribute "DistributedArrays.jl" distribute() function to distribute the raw data matrix into specified distribution
function mydistribute(A::AbstractArray, mydist::AbstractArray;
                    myprocs=workers()[1:min(nworkers(), maximum(size(A)))])
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
#---------------------------------------------------------------------------------------------------------------------------------	
## Serial Pulse Compression function to run on multiple cores
@everywhere function mycorrelate(n_pulses::Int, ref_sig::AbstractArray, raw_data::AbstractArray)
# Array Initialization
d = zeros(n_pulses, 2048)*im;
pl = plan_fft(raw_data[1,:]);

# Loop to correlate frequency domain refSig and data array, and convert back to time domain
for n = 1:n_pulses
	d[n,:] = (ref_sig).*(pl*raw_data[n,:]);
end
pl = plan_ifft(d[1,:]);
for n = 1:n_pulses
	d[n,:] = (pl*d[n,:]);
end
#println(typeof(d))
return d
end
#---------------------------------------------------------------------------------------------------------------------------------	
## Parallel Pulse Compression function
function MCPC(nPulses_::Int, rawDataFile_::AbstractString, refSigFile_::AbstractString)

# Array and Variable Initializations
P = nprocs()-1;  
PCpart = convert(Int, floor(nPulses_/P));
Ddata = zeros(nPulses_,2048);


# Read the data from the bin file
dataRaw = readRawData(nPulses_,rawDataFile_);	#"/home/stephaniejonkers/NetRAD/e10_10_21_1927_55_P1_1_130000_S0_1_2047_node3.bin"
# Read the chirp from the txt file
refSig = readRefSig(refSigFile_);		#"/home/stephaniejonkers/NetRAD/refSigN3N3Pl5000.txt"

# Distribute the raw data split in the row dimension
Ddata = mydistribute(dataRaw, [P,1]);

# Clear local raw data memory
dataRaw = 0;
gc()

# Move all mycorrelate() to all Ddata local to each worker and store it to a list, note that @spawnat output is a remote reference
refs = [(@spawnat (procs(Ddata))[w] mycorrelate(PCpart, refSig, localpart(Ddata))) for w=1:P];

# Clear distributed data memory
Ddata = 0;
gc()	

# Fetch all the remote references
DdataArray = Array{Any,1};
DdataArray = pmap(fetch,refs);

# Move the array of processed chunks back into local memory
DdataArray = convert(Array, DdataArray);
DdataRefsTemp = DdataArray[1];

# reconstruct the array from list of arrays
for w = 2:P
DdataOut = vcat(DdataRefsTemp, DdataArray[w]);
DdataRefsTemp = DdataOut;
end
DdataArray = 0;
DdataOut = 0;
gc()

## Create a file for Memory Mapping
myfile = open("NR_PCdata.bin", "w+");

# Setting the dimensions of the Array in binary file
write(myfile, size(DdataRefsTemp,1)::Int64);
write(myfile, size(DdataRefsTemp,2)::Int64);

write(myfile, DdataRefsTemp);
close(myfile);
end

## Run the Pulse Compression function
DdataIn = MCPC(80000, "/home/stephanie/Desktop/e10_10_21_1927_55_P1_1_130000_S0_1_2047_node3.bin", "/home/stephanie/Desktop/refSigN3N3Pl5000.txt")
gc()



