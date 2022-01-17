# Trackers-Alignment

The directory contains:
-A file for the setup simulation
-A file for function minimization

////////////////Simulation Setup///////////////////
a)The centre of trackers is set in function Initialize()

b)The intersection between the trackers' plane and the muon line is decided in Projection

c)The intersection points are called:
-"true" if they are given directly in output
-"gaussian" if a smearing is applied to simulate tracker resolution. The function for this task is "tracker_resolution"

d)In order to save the data, launch function "simulation". The number of tracks need to be adjusted inside this function. Now it is set to 8000



////////////////Residuals Minimization///////////////////


