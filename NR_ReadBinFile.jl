## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to read and process the data from the bin file

function readRawData(N_, file_name_)
binFile 	= open(file_name_);

# Julia is column major so the read process reads into coloums before rows 	
dataArray_ 	= read(binFile, UInt16, 2048, N_);	
close(binFile);

# Transpose matrix to get (pulse x rangebins)
dataArray_	= dataArray_';	

# Remove bias introduced by ADC			
data_		= dataArray_ - mean(dataArray_[1,:]);   
data_ 		= convert(Array{Float32}, data_)

# Free memory
dataArray_	= 0;
gc()

return data_
end
