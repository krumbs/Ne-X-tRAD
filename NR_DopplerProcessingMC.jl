## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to perform Doppler processing on match filtered NetRad data

## Function that creates a Doppler filterbank 
function NR_DopplerFilterBank(nPulses_::Int64, dftSize_::Int64, rawData_::AbstractArray)
	local sp;
	local pl;
	sp = zeros(Complex128, 2048, nPulses_);
	pl = fft_plan(rawData_[1,:])
	for n = 1:2048
		for m = 1:256:nPulses_
			sp[n,m:m+(dftSize_-1)] = fftshift(fft(rawData_[n,m:m+(dftSize_-1)]));
		end
	end
end
#=
## Function to parallelize the transposition across different CPUs
function NR_CornerTurn()
end
## Function to perform multi core doppler processing
function MCDP()
=#

	
