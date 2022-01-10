## Files

### compression_data.py

This file collects data into pandas data frames and throws it into csv files. It has a method for collecting data from the ucsc gtf files and from the encode api. From the gtf data file it generates csv recording the effect of set size, delimiter size, and method on the space taken by the array of start and end values for each genomic tuple.

### analysis.py

analysis.py uses the data generated by compression_data.py to create a several different graphs displaying different ascpects of the data.

### BITS_base.py

A reconstruction of the BITS algorithm.

### BITS_succinct1.py

The beginning of trying to use the different methods to change the array storing the start and end values into a succint data structure and using the BITS algorithm with this. The file is not finished, however, and was abandoned after finding the problems associated with method 3, and 4.

### check_prob_bistring.py

This file is used for getting the probability that a sub string of size "x" will occur in a string. It does this empirically, by generating all possible strings of a set size and then counting which strings have the substring. This can take too long if the string sizes are too large, though.

### gen_rnd_data.py

This file was meant for the BITS_base.py and BITS_succinct1.py files. It generates random start and stop values that the BITS algorithm can work on.