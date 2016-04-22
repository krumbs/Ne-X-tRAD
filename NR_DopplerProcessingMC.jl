## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to perform Doppler processing on match filtered NetRad data

include("NR_PulseCompressMC.jl")
export MCPC
export NR_Distribute

## Function that creates a Doppler filterbank 
@everywhere function NR_DopplerFilterBank(np_::Int64, chunk_::Int64, ds_::Int64, inputD_::AbstractArray)
	local sp;
	local pl;
	sp = zeros(Complex128, chunk_, np_);
	pl = plan_fft(inputD_[1,1:ds_]);
	for n = 1:chunk_
		for m = 1:ds_:np_
			sp[n,m:m+(ds_-1)] = fftshift(pl*(inputD_[n,m:m+(ds_-1)]));
		end
	end
	return sp
end

#using PyPlot
## Function to perform multi core doppler processing
function MCDP(np_::Int64, dftSize_::Int64, dFile_::AbstractString, rFile_::AbstractString)
	local data; 
	local PCdata;
	local p;
	local DPchunk;	
	
	## Perform multi-core pulse compression
	PCdata = MCPC(np_,dFile_,rFile_,true);
	
	## Convert back to local array
	Dd = convert(Array, PCdata);
	p = nprocs()-1;
	temp = Dd[1][1];
	for q = 2:p
		data = vcat(temp,Dd[q][1]);
	end

	## Distribute transformed array
	data = data';
	Ddata = zeros(2048,np_);
	Ddata = NR_Distribute(data, [p,1]);
	DPchunk = convert(Int, floor(2048/p));

	refs = [(@spawnat (procs(Ddata))[w] NR_DopplerFilterBank(np_, DPchunk, dftSize_, localpart(Ddata))) for w=1:p];
	DdataArray = Array{Any,1};
	DdataArray = pmap(fetch,refs);

	## Draw
	#=
	for k = 1:4
		for m = 1:dftSize_:np_
			imshow(log10(abs(DdataArray[k][1][:,m:m+(dftSize_-1)])), aspect = "auto");
			draw();
		end
	end
	=#
	return DdataArray
end
tic()
output = MCDP(256*80, 256, "/home/jnkste003/thesis_code/Ne-X-tRAD/data/e10_10_21_1927_55_P1_1_130000_S0_1_2047_node3.bin","/home/jnkste003/thesis_code/Ne-X-tRAD/data/refSigN3N3Pl5000.txt")
toc()
