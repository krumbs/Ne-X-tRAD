## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to read and process the chirp from the txt file 

# Add DSP package 
Pkg.add("DSP")
using DSP

function readRefSig(file_name_::AbstractString)

refSigRaw_  	= readcsv(file_name_, Float32);	
#---------------------------------------------------------------
## Process refernece signal

refSigRaw_ 	= refSigRaw_/sqrt(mean(abs(refSigRaw_).^2));

# Apply Hanning window
refSigRaw_	= (hanning(size(refSigRaw_,2)))'.*refSigRaw_; 

# Zero fill to size of range line
refSig_	 	= [refSigRaw_ zeros(1,2048-size(refSigRaw_,2))];

# Clear memory
refSigRaw_ = 0;

# Fast Fourier tranform
refSig_	 	= fft(refSig_);

# Complex Conjugate
refSig_		= conj(refSig_);

return refSig_
end
