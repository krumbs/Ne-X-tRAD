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
function NR_Distribute(A::AbstractArray, mydist::AbstractArray; myprocs=workers()[1:min(nworkers(), maximum(size(A)))])
	owner = myid()
	rr = RemoteRef();
	put!(rr, A);
	d = DArray(size(A), myprocs, mydist) do I
		remotecall_fetch(()->fetch(rr)[I...], owner);
	end
	return d
end
#---------------------------------------------------------------------------------------------------------------------------------	
## Serial Pulse Compression function to run on multiple cores
@everywhere function NR_Correlate(np_::Int, refSig_::AbstractArray, rawData_::AbstractArray; analytic_::Bool = false)
	## Array Initialization
	local d;
	local pl;
	d = zeros(np_, 2048)*im;
	pl = plan_fft(rawData_[1,:]);

	## Loop to correlate frequency domain refSig and data array, and convert back to time domain
	for n = 1:np_
		d[n,:] = (refSig_).*(pl*rawData_[n,:]);
	end
	## Flatten the negative half of the spectrum to create the analytic signal
	if (analytic_ == 1)
		local s;
		s = convert(Int64, np_/2);
		for i = 1:n_pulses
			d[i,1:s] = zeros(s,1);
		end
	end
	pl = plan_ifft(d[1,:]);
	for n = 1:np_
		d[n,:] = (pl*d[n,:]);
	end
	return d
end

#---------------------------------------------------------------------------------------------------------------------------------	
## Parallel Pulse Compression function

function MCPC(np_::Int, rawDataFile_::AbstractString, refSigFile_::AbstractString, analytic_::Bool)
	## Array and Variable Initializations
	local PCchunk;
	local Ddata;
	local dRaw;
	local rSig;
	local refs;
	local DdataArray;

	P = nprocs()-1;  
	PCchunk = convert(Int, floor(np_/P));
	Ddata = zeros(np_,2048);

	## Read the data from the bin file
	dRaw = readRawData(np_,rawDataFile_);
	
	## Read the chirp from the txt file
	rSig = readRefSig(refSigFile_);

	## Distribute the raw data split in the row dimension
	Ddata = NR_Distribute(dRaw, [P,1]);

	## Clear local raw data memory
	dRaw = 0;
	gc()

	refs = [(@spawnat (procs(Ddata))[w] NR_Correlate(PCchunk, rSig, localpart(Ddata)), analytic_) for w=1:P];

	## Fetch all the remote references
	DdataArray = Array{Any,1};
	DdataArray = pmap(fetch,refs);

	# Clear distributed data memory
	Ddata = 0;
	refs = 0;
	gc()
#=
	## Create a file for Memory Mapping

	myFile = open("NR_PCdata.bin", "w+");

	for w = 1:P
		write(myFile, DdataArray[w][1]);

	end
	DdataArray = 0;
	DdataOut = 0;
	gc()
	close(myFile);

=#
	return DdataArray

end
## Run the Pulse Compression function
#tic()
#data = MCPC(130000, "/home/jnkste003/thesis_code/Ne-X-tRAD/data/e10_10_21_1927_55_P1_1_130000_S0_1_2047_node3.bin", "/home/jnkste003/thesis_code/Ne-X-tRAD/data/refSigN3N3Pl5000.txt", true)
#toc()
