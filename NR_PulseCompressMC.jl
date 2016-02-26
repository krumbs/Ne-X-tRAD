## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to perform parallel Pulse Compression on NetRAD data

Pkg.add("DistributedArrays")
using DistributedArrays

include("/home/stephanie/Desktop/NR_ReadBinFile.jl") 
include("/home/stephanie/Desktop/NR_ReadRefSig.jl")
export readRawData
export readRefSig

#---------------------------------------------------------------------------------------------------------------------------------
## Altered distribute "DistributedArrays.jl" distribute() function to distribute the raw data matrix into specified distribution
function NR_Distribute(A::AbstractArray, mydist::AbstractArray; myprocs=workers()[1:min(nworkers(), maximum(size(A)))])
    owner = myid()
    rr = RemoteRef();
    put!(rr, A);
    d = DArray(size(A), myprocs, mydist) do I
        remotecall_fetch(()->fetch(rr)[I...], owner);
    end
    # Ensure that all workers have fetched their localparts.
    for chunk in d.chunks
        wait(chunk);
    end
    return d
end
#---------------------------------------------------------------------------------------------------------------------------------	
## Serial Pulse Compression function to run on multiple cores
@everywhere function NR_Correlate(nPulses_::Int, refSig_::AbstractArray, rawData_::AbstractArray; analytic_::Bool = false)
	# Array Initialization
	local d;
	local pl;
	d = zeros(nPulses_, 2048)*im;
	pl = plan_fft(rawData_[1,:]);

	# Loop to correlate frequency domain refSig and data array, and convert back to time domain
	for n = 1:nPulses_
		d[n,:] = (refSig_).*(pl*rawData_[n,:]);
	end
	# Flattens the negative half of the spectrum to create the analytic signal
	if (analytic_ == 1)
		local s;
		s = convert(Int64, nPulses_/2);
		for i = 1:n_pulses
			d[i,1:s] = zeros(s,1);
		end
	end
	pl = plan_ifft(d[1,:]);
	for n = 1:nPulses_
		d[n,:] = (pl*d[n,:]);
	end
	return d
end

#---------------------------------------------------------------------------------------------------------------------------------	
## Parallel Pulse Compression function

function MCPC(nPulses_::Int, rawDataFile_::AbstractString, refSigFile_::AbstractString, analytic_::Bool)
	# Array and Variable Initializations
	local PCchunk;
	local Ddata;
	local dataRaw;
	local refSig;
	local refs;
	local DdataArray;
	#local DdataOut;

	P = nprocs()-1;  
	PCchunk = convert(Int, floor(nPulses_/P));
	Ddata = zeros(nPulses_,2048);

	# Read the data from the bin file
	dataRaw = readRawData(nPulses_,rawDataFile_);
	
	# Read the chirp from the txt file
	refSig = readRefSig(refSigFile_);

	# Distribute the raw data split in the row dimension
	Ddata = NR_Distribute(dataRaw, [P,1]);

	# Clear local raw data memory
	dataRaw = 0;
	gc()

	refs = [(@spawnat (procs(Ddata))[w] NR_Correlate(PCchunk, refSig, localpart(Ddata)), analytic_) for w=1:P];

	# Fetch all the remote references
	DdataArray = Array{Any,1};
	DdataArray = pmap(fetch,refs);

	# Clear distributed data memory
	Ddata = 0;
	refs = 0;
	gc()

	## Create a file for Memory Mapping

	myFile = open("NR_PCdata.bin", "w+");

	for w = 1:P
		write(myFile, DdataArray[w][1]);

	end
	DdataArray = 0;
	DdataOut = 0;
	gc()
	close(myFile);

#=
	return DdataArray
=#
end
## Run the Pulse Compression function
data = MCPC(2048*35, "/home/stephanie/Desktop/e10_10_21_1927_55_P1_1_130000_S0_1_2047_node3.bin", "/home/stephanie/Desktop/refSigN3N3Pl5000.txt", true)






