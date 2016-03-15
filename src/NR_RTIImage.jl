## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to create the radar image from Pulse Compressed 
## Note make sure you have Python and Matplotlib for Python installed before adding PyPlot package

using PyPlot

function RTIplot(binFile_::AbstractString)
file = open(binFile_);

# Get dimensions of stored matrix 
dim1 = read(file, Int64);
dim2 = read(file, Int64);

println("File data dimensions: (", dim1, ",", dim2, ")");

# Read in binary data as 2 Complex{Float64}
imageData = read(file, Complex128, (dim1,dim2));

# Create figure 
fig = imshow(log(abs(imageData)), aspect="auto");
xlabel("Range bins");
ylabel("Pulse number");
title("Ne(X)tRAD RTI Plot");

#Save figure to file
savefig("MCPC_RTIplot.png");

close(fig);
end
