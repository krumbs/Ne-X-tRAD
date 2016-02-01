# Ne-X-tRAD
Repository for all Julia infrastructure software to acquire, prepare radar data for Sea Clutter Anaylsis of NeXtRAD project.  

# Overview
Project NetRAD is an open source, multistatic radar system started to provide bistatic sea clutter data to aid the research and undestanding. This repository stores the scripts to acquire (read data from HDF5 files), process (match filter, Doppler filter, suppress interference signals) of the captured radar data. Development is done using Julia-lang, a language designed for scientific computing.

# Notes
This repository forms part of a greater project and is still heavily under construction.

# Pulse Compression 
Pulse compression, often also called match filtering is a method of increasing the 'signal-to-noise' ratio of a dataset using a matched filter. In radar we use the transmitted (or reference) signal as the matched filter. This reference is then convolved with each sampled raw input data array.

Once you have found the dataset containing the raw radar data. If it is in a binary (.bin) file then start by reading in the data and converting it to voltages using the NR_ReadBinFile.jl script. The same must be done for the reference signal which is typically stored in a text (.txt) file by using the NR_ReadRefSig.jl file.

Now you are ready to perform pulse compression on the data. PulseCompressMC.jl does this using multiple cores. Start julia in the terminal using (julia -p 4) to start it on an internal cluster of 4 cores. Then include NR_PulseCompressMC.jl to store the matched filtered data to a binary file.

Finally, NR_RTIImage.jl creates a range-time intensity plot of the filtered data, to get a "quick-look" at the data that you have chosen to analyse. 
