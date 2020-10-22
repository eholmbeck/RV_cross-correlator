<h2>RV_cross-correlator</h2>

<h4>Owner: Erika M. Holmbeck <br>
Last updated: 22 Oct 2020</h4>


This is a gui to calulate stellar radial velocities by cross-correlation. It relies on tk and has been thoroughly tested with Python 2.7.

Call `python gui.py` to run the code. Use Command-S to select the **S**cience spectrum to be RV-corrected, and use Command-T to select the **T**emplate spectrum to be used for the cross-correlation. The Template spectrum must have enough information in the header to extract the BCV (observatory/site lat, long, and alt; UT-DATE; etc.).

