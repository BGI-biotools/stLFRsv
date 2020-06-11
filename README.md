[TOC]

# stLFRsv
## Introuction:
Structure variation(SV) pipeline for stLFR linked reads data.
This tool is applicable to stLFR technology and similar linked read data. Currently running on stLFR data, theoretically it is also applicable to other linked read data. You can analyze and test the barcode by converting it to the `read_id#XXX_XXX_XXX` format on read ID.

**Share barcode** information is used to detect breakpoint signals of structural variation (SV), such as: equilibrium translocation, INV, partial missing repetitions, and more complex structural breakpoints, which can be combined with CNV results and phase results. The accuracy of deriving the true structural variation of chromosomes is limited by the distribution density and length of linked reads on DNA molecules. For example, stLFR kits built with 1.5ng starting volume(about 30x significant depth) can produce data that can guarantee SV detection accuracy above 20K.
