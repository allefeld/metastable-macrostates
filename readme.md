# Identification of metastable macrostates

Matlab code for C. Allefeld, H. Atmanspacher, and J. Wackermann, Mental states as macrostates emerging from brain electrical dynamics.  *Chaos*, 19(1):015102, 2009.

This code attempts to automatically identify metastable states (or "almost invariant sets") as macrostates from data describing a microstate dynamics. The main function is `almostInvariant`; if the data are not already discretized, `discretize` can be used to do this.  `transitionMatrix` and `pccap` are helper functions for `almostInvariant`.  The `example` script demonstrates how to use the functions, and also uses `scatterbox` to plot the sample data.

