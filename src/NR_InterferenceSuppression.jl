## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to suppress WiFi interference from NetRAD data files using adaptive LMS filter

# Load necessary packages
Pkg.add("MAT")
using MAT


# Open MAT file with contaminated data
file = matopen("e14_07_30_1330_37_P1_1_10000_S0_1_2047_node3_MF_refsig.mat");

# Read the data into a julia matrix
data = read(file, "Data_matched"); # Note you need to know the name of the matrix variable in this file

# Close the file
close(file)

# Variable definitions: these variables can be altered depending on the dataset, to give 
# better results
p = 10000;		# number of pulses
n = 60;			# number of filter weights
d = 3;			# size of the delay
mu = 1.8e-6;

# Variable definitions: these chould not be changed
inp = data[1:p,:];	# input data
Z = zeros(p, 1024)*im;
W = zeros(n)*im;

# Loop through constant range lines
for k = 1:1024
	# Loop through j-th range line
	for j = d:p-n
		X = inp[j+1:j+n,k];
		er = inp[j-d+1,k] - dot(X,W);
		W = W+2*mu*er*conj(X);
		Z[j,k] = er;
	end
end	

# Import matplotlib.pyplot
using PyCall
@pyimport maplotlib.pyplot as plt

println("Printing contaminated_data.png to disk...")
# Save an image of the contaminated data
plt[:savefig]("contaminated_data.png")
println("Printing interference_suppressed_data.pngor a Python file-like to disk...")
# Save an image of the interference suppressed data
plt[:savefig]("interference_suppressed_data.png")
