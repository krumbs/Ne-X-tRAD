## Author: Stephanie Jonkers
## Institution: University of Cape Town MSc(Eng)
## Code to read NetRAD radar data and radar parameters from HDF5 file 

Pkg.add("HDF5")
Pkg.add("Compat")
using HDF5
using Compat

function readHDF5(fileName_::AbstractString)
local dat
local a 
local d
# Open HDF5 file
file = h5open(userInput, "r+");

println("Displaying file structure:");
dump(file);

# g is the root group, g1 and g2 are primary group objects
g = root(file);
g1 = g["Data"];
g2 = g["Temp"];

# g1 contains the radar data
d = g1["SamplePerPulse_x_PulseNo"];
data  = read(d);

# g2 contains the radar parameters which are stored as attributes 
a = attrs(g2);

# Stores the attributes (name, value) in a dictionary
aArray = @compat[x=>read(a[x]) for x in names];

close(file);

# Prints the attributes
for (k,v) in aArray
  println("$(k) = $(v)");
end

return data
end
