To make parallel version type "make DelayNFE_Example1_PAR".  **ONLY WORKS IF OPENMP IS INSTALLED**


DelayNeuralField problem assumes delays of the form:

tau_ij = tau_c + |x_i-x_j|/c

where tau_c is a constant delay, |x-xp| is the distance between two points, and c is transmission delay speed.


To minimise the number of delays and prevent duplications, we take advantage of the symmetry property and only store the delay if i<j. The delay vector is populated as

delays = [tau_c, tau_01, tau_02, ... , tau_0n, tau_12, tau_13, ... , tau_1n, ... , tau_(n-1)n]


To access the delays when looping over i and j indices,

if (i==j) then the distance is 0 and therefore we access tau_c
if (i<j) then we access tau_ij
if (i>i) then we access tau_ji


To map from i and j to the correct entry in the delays vector we use the mapping

if (i==j) access delays(0)
if (i<j) access delays((2*i*n - i^2 + 2*j - 3*i - 2)/2 + 1)
if (i>j) access delays((2*j*n - j^2 + 2*i - 3*j - 2)/2 + 1)
