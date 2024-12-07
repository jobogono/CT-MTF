# CT-MTF
Author: Joshua Bognar

Two ImageJ macros demonstrate the differences in a generated MTF from a PSF using different approaches. One method uses integration along a single axis to generate an LSF. The other method uses radial averaging around a PSF to generate an LSF.
Example scans are included for ease of demonstration.

Macro Files:
MTF - Integration.ijm
MTF - Radial Averaging.ijm

To run a macro in ImageJ, select "Plugins" -> "Macros" -> "Run..." and select the desired macro file. The user will then be prompted to select a folder.
All DICOM images in this folder will be opened for assessment. This permits averaging an LSF over multiple slices. 
The user will first be prompted to select slices for analysis (type "1" if only using a single slice). 
The user will then be prompted to draw a ROI around the wire object. This is only used for locating the peak value. A different ROI will be generated for analysis around the peak value.

Image folders:
Standard Resolution
High Resolution

Each image folder contains a single slice of a standard or high resolution scan.
