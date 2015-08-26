# force-match-wgrouping
New version of forceMatch.py that supports grouping of atoms into single sites.
Consists of 5 necessary .py files and a .config file:
  - forceMatch.py
  - forceMatchInit.py
  - forceMatchAtomGroup.py
  - forceMatchProcess.py
  - forceMatchPlot.py
  - forceMatch.config

Required Python dependencies:
  - numpy
  - MDAnalysis
  - matplotlib
  
<b>Instructions</b>

1) Copy any .dcd, .force.dcd, and .psf into the program directory.

2) Edit 'forceMatch.config' to reflect the file names and analysis parameters you wish to use.

3) You may change the name of 'forceMatch.config' in the event you have data from multiple simulations, and therefore, multiple .config files. (ie: 'NaCl.config, KCl.config, rada.config')

4) Run:
      
      user@local:~$  python forceMatch.py <configFileName>
