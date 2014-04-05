import numpy as np
# files I want to read 
data_block_size = 100



the_array = np.loadtxt('kin_energy.dat')
cor_len = np.size(the_array)-np.size(the_array)%data_block_size
right_size_array = the_array[:cor_len].reshape(data_block_size,-1)
averages = np.mean(right_size_array, axis = 0)
variances = np.std(right_size_array,axis = 0)
cv = 1.0/(2.0/(3.0*864.0)-(variances/averages)**2)

print np.mean(cv)
print np.std(cv)/np.sqrt(len(averages))


