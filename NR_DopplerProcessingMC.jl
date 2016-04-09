## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to perform Doppler processing on match filtered NetRad data

## Function that creates a Doppler filterbank 
@everywhere function NR_DopplerFilterBank(np::Int64, ds::Int64, rawData_::AbstractArray)
	local sp;
	local pl;
	sp = zeros(Complex128, 2048, np);
	pl = fft_plan(rawData_[1,:]);
	for n = 1:2048
		for m = 1:256:nPulses_
			sp[n,m:m+(ds-1)] = fftshift(fft(rawData_[n,m:m+(ds-1)]));
		end
	end
end
#=
## Function to parallelize the transposition across different CPUs
function NR_CornerTurn()
end
=#
using PyPlot
## Function to perform multi core doppler processing
function MCDP(nPulses_::Int64, dftSize_::Int64, FileName_::AbstractString)
	local d;
	local p;
	local DPchunk;	

	myFile = open(FileName_);
	d = read(myFile, Complex128, 2048, nPulses_);
	d = d';
	
	p = nprocs()-1;  
	DPchunk = convert(Int, floor(nPulses_/p));
	Ddata = zeros(nPulses_,2048);
	Ddata = NR_Distribute(d, [p,1]);
	refs = [(@spawnat (procs(Ddata))[w] NR_DopplerFilterBank(DPchunk, dftSize_, localpart(Ddata))) for w=1:p];
	
	DdataArray = Array{Any,1};
	DdataArray = pmap(fetch,refs);


	## Draw
	for k = 1:4
		for m = 1:dftSize_:nPulses_
			imshow(log10(abs(DdataArray[k][1][:,m:m+(dftSize_-1)])), aspect = "auto");
			draw();
		end
	end
end

MCDP(2048*35, 256, "NR_PCdata.bin")
