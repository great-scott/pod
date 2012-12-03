pod~
===
v0.1 - Scott McCoid, Jay Clark, Greg Tronel

Overview
--------
pod~ is a perceptually influenced real-time onset detection external Pure Data (Pd) object.
The system does not assume an underlying pattern for onset occurrences and tries to be
robust across musical signals. The object differs from previous attempts by applying a
more comprehensive representation of the peripheral auditory system by using filter
models of the outer and middle ear acoustics.

We use band-wise processing to simulate the response of the Basilar Membrane.
Frequency values are weighted and scaled according to Bark frequency scale and additionally
weighed according to the equal-loudness contour. We then use a Spectral Flux based method
for determining where an onset occurs.

Setup
-----
Mac: An Xcode 4 project is provided to build the object.
- Specific path settings will need to be changed to automate copying the pod~.darwin
  target to the appropriate location.
- Additional settings will need to change to automatically launch the Pd binary and
  attach to that process (in order to debug)

Windows and Linux: It's definitely possible to build using another system. It should just
be a matter of downloading the source files. This link might potentially have a makefile 
template to use: http://puredata.info/docs/developer/MakefileTemplate
