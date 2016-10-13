# Analysis: Accurately Assessing The Global Quality of NMR-derived RNA Structures Using Chemical Shifts

NMR spectroscopy has emerged as a powerful tool for determining the complex architecture of non-coding RNAs, but assessing the quality of NMR-derived structures remains problematic in the field of RNA structural biology. We demonstrate that chemical shifts can serve as conformational “fingerprints” that can accurately assess the global quality of NMR structures of RNA. This was accomplished by carefully comparing observed and computed chemical shifts for pairs of NMR structures of the same RNA, one of which is known to be a better set of atomic coordinates. We discovered that in 7 out of 8 cases, we were able to correctly identify the more “accurate” NMR structure of the studied pairs. The promising results described in this report should pave the way for the development of tools that utilize chemical shifts to assess the global quality of NMR structures of RNA. As an initial step in this direction, we developed a graphical analysis and quality assessment tool, PyShifts, which we make freely available at https://github.com/atfrank/PyShifts.

- analysis.R -- R script that carryouts out all the analysis represent in the manuscript
- library.R -- R script that defines custom functions used in the analysis
- coors -- folder containing coordinate files from which chemical shifts are computed
- data -- folder containing observed and predicted chemical shifts used in the analysis
- errors -- folder containing computed errors
- figures -- folder containing generated figures