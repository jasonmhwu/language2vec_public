This is the instruction on how to replicate the main results of our preprint paper.

1. Download the entire folder from OSF (TODO: specify website?)
2. open code/setup.m, specify your spm12 path, and execute code/setup.m
3. run main.m

Dataset Structure:
1. aal: stores the mask for the 8 semantic ROIs
2. behavioral: stimulus and behavioral response per run per subject
3. code: stores all codes
4. group_analysis: stores intermediate results for group_level analysis
5. harvard_oxford: contains mask for the 48 Harvard-Oxford ROIs
6. public_repositories: contains the fdr_bh function written by David Groppe and raizadalab-utilities folder by Dave Kleinschmidt.
7. subjects: preprocessed beta images (from SPM12) are stored under SPM_analysis_motionReg

References:
1. David Groppe (2020). fdr_bh (https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh), MATLAB Central File Exchange. Retrieved February 27, 2020.
