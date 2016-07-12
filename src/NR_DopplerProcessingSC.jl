## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to perform serial pulse Doppler processing on match filtered NetRad data


## Import pulse compression script
include("NR_PulseCompressSC.jl")					

## Pulse-Doppler Processing
function NR_DopplerFilter(np_::Int64, ds_::Int64, inputData_::AbstractString, inputR_::AbstractString)
	## Perform pulse compression
	data = SCPC(np_,inputData_,inputR_,true)
	
	## Define local variables
	local a;
	local s = zeros(Complex64, 2048, np_)

	## Transpose the matrix for FFT along constant range
	a = data';

	## Clear memory
	data = 0;
	gc()
	
	P = plan_fft(s[1,1:ds_]);

	## FFT the data along the slow-time to obtain the range-Doppler matrix
	for n = 1:2048
		for m = 1:ds_:np_
			s[n,m:m+ds_-1] = fftshift(P*(a[n,m:m+ds_-1]));
		end
	end
return s
end

#----------------------------------------------------------------------------------------------------
# For Skippy:
# parameter 1 of NR_DopplerFilter() is the number of range lines (pulses) you consider. It's best if you choose the number of range lines to be a multiple of the DFT size.
# parameter 2 of NR_DopplerFilter() is the DFT size. Best to choose from the following (128, 256, 512)
# it might take a while to save to file
#-----------------------------------------------------------------------------------------------------

output = NR_DopplerFilter(256*200, 256,"/home/jnkste003/thesis_code/Ne-X-tRAD/data/e10_10_21_1927_55_P1_1_130000_S0_1_2047_node3.bin", "/home/jnkste003/thesis_code/Ne-X-tRAD/data/refSigN3N3Pl5000.txt")

## Create and save to binary file
myfile = open("DopplerProcessedData.bin", "w+");

# Setting the dimensions of the Array in binary file
write(myfile, size(output,1)::Int64);
write(myfile, size(output,2)::Int64);

write(myfile, output);
close(myfile);

## Plot one cut
#image = imshow(log10(abs(s[:,1:1+255])), aspect = "auto");
#using PyCall
#@pyimport matplotlib.pyplot as plt
#plt[:savefig]("range_doppler.png");	
