// Author: Joshua Bognar
// Generation of an MTF from a PSF. An LSF is generated usial radial averaging method.
// Permits averaging over multiple slices
// run macro, then select folder. All CT slices in folder will be opened. Slices can then be selected.

// discrete fourier transform for real input. Outputs complex absolute to provide real output only.
function dft(x_array)
{
	n = x_array.length;
	X_array = newArray(n);
	for (k=0; k<n-1; k++)
	{
		X_k = 0;
		sum_real = 0;
		sum_imag = 0;
		for (m=0; m<n-1; m++)
		{
			sum_real += x_array[m]*Math.cos(2*PI*k*m/n);
			sum_imag += x_array[m]*Math.sin(2*PI*k*m/n);
		}
		complex_abs = Math.sqrt(Math.sqr(sum_real) + Math.sqr(sum_imag));
		X_array[k] = complex_abs;
	}
	return X_array;
}

// takes an array and gives the position in the array which contains the minimum value
function minpos(array)
{
	j=1024;
	t=0;
	for (i=0; i<array.length; i++)
	{
		if (array[i] < j)
		{
			j = array[i];
			t = i;
		}
	}
	return t;
}

// takes an array and a value and identifies the array position of the closest two values in the array either side of the value
function nearest(xn, array)
{
	negative_array = newArray(array.length);
	positive_array = newArray(array.length);
	for (i=0; i<array.length; i++)
	{
		diff = array[i]-xn;
		if (diff < 0)
		{
			negative_array[i] = Math.abs(diff);
			positive_array[i] = 1024;
		}
		else if (diff > 0)
		{
			negative_array[i] = 1024;
			positive_array[i] = diff;
		}
		else
		{
			negative_array[i] = 0;
			positive_array[i] = 0;
		}
	}
	x0 = minpos(negative_array);
	x1 = minpos(positive_array);
	
	return newArray(x0, x1);
}

// linear interpolation function. Find yn for desired xn with known x0, x1, y0, y1
function interpolate(xn, x_array, y_array)
{
	array_pos = nearest(xn, x_array);
	x0 = x_array[array_pos[0]];
	x1 = x_array[array_pos[1]];
	y0 = y_array[array_pos[0]];
	y1 = y_array[array_pos[1]];
	yn = y0 + (xn-x0)*(y1-y0)/(x1-x0);
	// handle situation where xn == x0 == x1 (i.e., interpolation not required)
	if (isNaN(yn))
	{
		return y0;
		
	}
	return yn;
}

// Function takes a value n and returns true if n is odd or false if n is even
function isOdd(n)
{
	if (n % 2 == 0)
	{
		return false;
	}
	else
	{
		return true;
	}
}

// Prompt user to select a directory and open all files as a stack
inputDir = getDirectory("Choose input directory");
run("Image Sequence...", "dir=["+inputDir+"]");

// Prompt user to identify the first and last slice to use for calculation of MTF
waitForUser("Slice selection", "Identify the first and last slice number and then select OK");
run("Make Substack...");

// Extract pixel size / pitch - assumes size is the same in the x- and y- axes. Currently assumes unit is mm
getPixelSize(unit, pw, ph);
nyquist = 10/(2*pw); //nyquist in cm-1

// Define a physical size in mm for the ROI selection to generate the PSF.
// Value of 7 is currently set for use with CATPHAN so that objects other than wire are not in the ROIs
physical_width = 7;
physical_height = 7;

// Convert physical dimensions to pixel dimensions.
// These values are forced to be even to ensure a later division by 2 results in a whole number
roi_w = Math.round(physical_width/pw);
roi_h = Math.round(physical_height/pw);
if (isOdd(roi_w))
{
	roi_w++;
}
if (isOdd(roi_h))
{
	roi_h++;
}

// Define roi size in px for background calculation 
bg_w = roi_w*2;
bg_h = roi_h*2;

// Define number of samples and pixel resampling frequency
samples = 1024;
resample_freq = 0.25;

// Define top level arrays for PSF
multislice_psf = newArray(samples);

// Prompt user to draw a rectangle around the MTF object and get the bounds. These bounds will be used for each slice
// to identify the max value location
waitForUser("Draw a rectangle around the feature to be analysed and then select OK");
getSelectionBounds(rect_x, rect_y, rect_w, rect_h);

// Iterate through each slice to generate a PSF and add to top level arrays
for (i=1; i<=nSlices; i++)
{
	// Define arrays extracting distances and values within a roi
	distance_array = newArray((roi_w+1)*(roi_h+1));
	value_array = newArray((roi_w+1)*(roi_h+1));
	// Define arrays for reducing the distance and value arrays by taking the average of duplicate distances
	reduced_dist_array = newArray((roi_w+1)*(roi_h+1));
	reduced_val_array = newArray((roi_w+1)*(roi_h+1));
	
	
	// set the current slice for calculation
	Stack.setSlice(i);
	
	// Identify the position of the max pixel value in each slice
	makeRectangle(rect_x, rect_y, rect_w, rect_h);
	getRawStatistics(nPixels, mean, min, max);
	noise = max - (mean + 10);
	run("Find Maxima...", "noise=" + noise + " output=[Point Selection]");
	getSelectionBounds(x_max, y_max, w_max, h_max);
	
	// Create a large roi for the calculation of mean background px value
	bg_x = x_max - (bg_w/2);
	bg_y = y_max - (bg_h/2);
	makeRectangle(bg_x, bg_y, bg_w, bg_h);
	getStatistics(bg_nPixels, bg_mean, bg_min, bg_max);
	bg_roi_total = bg_mean * bg_nPixels;
	
	// Create small roi and get statistics for mean bg value
	roi_x = x_max - (roi_w/2);
	roi_y = y_max - (roi_h/2);
	makeRectangle(roi_x, roi_y, roi_w, roi_h);
	getStatistics(roi_nPixels, roi_mean, roi_min, roi_max);
	roi_total = roi_mean * roi_nPixels;
	
	// Calculate average background pixel value
	avg_bg = (bg_roi_total - roi_total) / (bg_nPixels - roi_nPixels);
	
	// within the small roi, iterate through all pixels and store the pixel value and distance from max pixel value
	roi_position = 0;
	for (x=0; x<=roi_w; x++)
	{
		for (y=0; y<=roi_h; y++)
		{
			distance = Math.sqrt(Math.sqr(roi_x+x-x_max) + Math.sqr(roi_y+y-y_max));
			value = getValue(roi_x+x, roi_y+y) - avg_bg;
			distance_array[roi_position] = distance;
			value_array[roi_position] = value;
			roi_position++;
		}
	}
	Array.sort(distance_array, value_array);
	
	// output unique distance values and mean values at duplicate distance values
	j=0;
	k=0;
	while (j < distance_array.length - 1)
	{
		slice1 = j;
		
		// handle situation for last set of distances to avoid array out of bounds error
		if (distance_array[j] == distance_array[distance_array.length - 1])
		{
			reduced_dist_array[k] = distance_array[j];
			Array.getStatistics(Array.slice(value_array,slice1,distance_array.length-1), smin, smax, smean);
		}
		else
		{
			while (distance_array[j] == distance_array[j+1])
			{
				j++;
			}
			reduced_dist_array[k] = distance_array[j];
			Array.getStatistics(Array.slice(value_array,slice1,j), smin, smax, smean);
		}
		
		if (isNaN(smean)) // mean returns NaN if taken for a single value
		{
			reduced_val_array[k] = value_array[j];
		}
		else
		{
			reduced_val_array[k] = smean;
		}
		j++;
		k++;
	}
	
	// interpolate reduced_val_array to have constant distance between values
	Array.getStatistics(reduced_dist_array, dmin, dmax);
	num_samples = (dmax-1) / resample_freq;
	resampled_vals = newArray(num_samples);
	for (j=0; j<num_samples; j++)
	{
		resample_position = j*resample_freq;
		resampled_vals[j] = interpolate(resample_position, reduced_dist_array, reduced_val_array);
	}
		
	// mirror interpolated data and average data to mutlislice array
	mirror_array = Array.reverse(Array.slice(resampled_vals,1));
	full_psf_array = Array.concat(mirror_array, resampled_vals);
	for (j=0; j < full_psf_array.length; j++)
	{
		multislice_psf[j] = multislice_psf[j] + (full_psf_array[j] / nSlices);
	}
	
}
mtf_array = dft(multislice_psf);
mtf_max = mtf_array[0];

pixpercm = 10/(pw*resample_freq);
freq = pixpercm/samples;

out_array_size = Math.ceil((1.1*nyquist)/(freq));

f_xArray = newArray(out_array_size);
f_yArray = newArray(out_array_size);
for (j=0; j<out_array_size; j++)
{
	f_xArray[j] = j*freq;
	f_yArray[j] = mtf_array[j]/mtf_max;
}

// find 50% and 10% MTF
fiftyMTF = interpolate(0.5, f_yArray, f_xArray);
tenMTF = interpolate(0.1, f_yArray, f_xArray);
fifty_label = "50%: " + toString(fiftyMTF);
ten_label = "10%: " + toString(tenMTF);

Table.create("MTF");
Table.setColumn("Frequency",f_xArray);
Table.setColumn("MTF",f_yArray);
Table.setLocationAndSize(100, 100, 200, 400);

Table.create("Metrics");
Table.set("50%",0,fiftyMTF);
Table.set("10%",0,tenMTF);
Table.setLocationAndSize(300, 100, 200, 200);

Plot.create("MTF", "Frequency (per cm)", "MTF", f_xArray, f_yArray);
Plot.setColor("red");
Plot.drawLine(0,0.5,fiftyMTF,0.5);
Plot.drawLine(fiftyMTF,0,fiftyMTF,0.5);
Plot.addText(fifty_label, 0.1, 0.5);
Plot.setColor("blue");
Plot.drawLine(0,0.1,tenMTF,0.1);
Plot.drawLine(tenMTF,0,tenMTF,0.1);
Plot.addText(ten_label, 0.1, 0.9);
Plot.setColor("black")

// Create array for neat output of LSF data
short_lsf = newArray(roi_w/resample_freq);
spatial_lsf = newArray(roi_w/resample_freq);
for (j=0; j<roi_w*Math.sqrt(2)/resample_freq; j++)
{
	short_lsf[j] = multislice_psf[j];
	spatial_lsf[j] = j*resample_freq;
}
Table.create("LSF");
Table.setColumn("Pixel", spatial_lsf);
Table.setColumn("LSF", short_lsf);
