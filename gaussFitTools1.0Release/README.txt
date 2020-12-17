101115 - S Holden

Install instructions for GaussFitTools
For detailed instructions please see Twotone3.1_UserManual.pdf.

1. Ensure previous versions of GaussFitTools are uninstalled and removed from the Matlab search path.

2. (Linux only) For Ubuntu 32 bit install, the following packages are required from the repositories: liblapack-dev libblas-dev libcfitsio3 libcfitsio30dev.For Ubuntu 64 bit install, the following are required: libcfitsio3 libcfitsio3-dev. Users of other distributions will need to find equivalent packages.

3. Load MATLAB, and navigate to the `gaussFitTools1.0Release' directory within the MATLAB environment. At the MATLAB command prompt, type
>> installGaussFitTools
