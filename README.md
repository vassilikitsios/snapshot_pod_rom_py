-----------------------------------------------------------------------
1) Code overview
-----------------------------------------------------------------------

Author:       Vassili Kitsios

The suit of codes calculates the proper orthogonal decomposition (POD) modes, (also known as principal component analysis, or empirical orthgonal functions, or singular value decomposition), and Koopman modes (also known as dynamic mode decomposiion, or principal oscillation patterns).
The POD modes are later to be used to build a reduced order model using these modes as basis. 
The codes are written in python, and uses readily available libraries.

If you are to use this code in your own projects please cite the following documents:

Kitsios, V., Cordier, L., Bonnet, J.-P., Ooi, A. & Soria, J., 2011, On the coherent structures and stability properties of a leading edge separated aerofoil with turbulent recirculation, Journal of Fluid Mechanics, Vol. 683, pp 395-416. http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=8378263

Kitsios, V., 2010, Recovery of fluid mechanical modes in unsteady separated flows, PhD Thesis, The University of Melbourne https://minerva-access.unimelb.edu.au/handle/11343/35705

-----------------------------------------------------------------------
2) List of files and directories with brief explanations
-----------------------------------------------------------------------

Each of the following directories contain a results directory (results), a directory containing the python source (src), and an images directory (images) with a gnuplot script.

1modes:

	Contains the temporal and spatial proper orthogonal decomposition (POD) modes and their eigenvalues.

	The temporal correlation between an temporal POD mode and itself is the variance (the eigenvalue), and between any other mode the temporal correlation is 0 - see equation 5.5 of Kitsios (2010).

	The inner product between a spatial POD modes and itself is 1, and between any other mode the inner product is 0 - see equation 5.7 of Kitsios (2010). 

	Calculation of POD modes is now a standard process. Determining these modes is not included in the present examples as this would require the individual snapshots to also be provided, which is not practical.


2rom_coefficients.no_calibration:

	The coefficients required to run a POD ROM are calculated from the spatial POD modes, the spatial mean field, and their spatial derivatives. Here no calibration is used.


3rom.no_calibration:

	A POD ROM is run with no calibration.

	The images directory compares the POD ROM temporal integration results to the temporal POD modes.


4rom_coefficients.constant_linear_calibration:

	The coefficients required to run a POD ROM are calculated from the spatial POD modes, the spatial mean field, and their spatial derivatives.

	Here the constant and linear coefficients have been calibrated using additional information from the temporal modes.

	Note the soursce is identical to "2rom_coefficients.no_calibration", it is only the input parameters that have changed.


5rom.constant_linear_calibration:

	A POD ROM is run with calibration of the constant and linear coefficients. The images directory compares the POD ROM temporal integration results to the temporal POD modes.

	Note the source is identical to that in "3rom.no_calibration", it is only the input parameters that have changed.

	Notice the improvement in agreement within the snapshot period, and improved stability after the snapshot period.

-----------------------------------------------------------------------
