# Analysis of prestim window with frequency sliding

The problem of investigating pre-stimulus EEG or MEG activity close to t=0 without altering the data in some way or mixing it with post-stimulus activity is a very common one. Potential solutions are to cut the data at t=0 and zero-pad afterwards; or to include post-stimulus data, but ignore (or interpret with caution) anything found in the period where post-stimulus activity might play a role. Another interesting method that has been tried is to mirror-pad the data at t=0. However, this also introduces spurious frequencies, which depend on the phase, as well as frequency at the point of mirroring. This presentation will illustrate these issues using simulated data, as well as speculate about possible improvements or alternatives.

