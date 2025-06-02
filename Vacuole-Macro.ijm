//Close all images that are currently open
close("*");

//Prevent display of images between calls
setBatchMode(true);

//Start the bio-formats extensions function
run("Bio-Formats Macro Extensions");

//Get list of all files and set output directory
inputDirectory = getDirectory("Choose the input directory");
dataset = getString("What is the name of the data set being analyzed", "Replace me with Date & Expt");
fileList = getFileList(inputDirectory);
	File.makeDirectory(inputDirectory + dataset + "_" + "Macro Results");
	resultsDirectory = inputDirectory + File.separator + dataset + "_" + "Macro Results" + File.separator;

CellsData = resultsDirectory + dataset + "_" + "Infected Cell Data.csv";
	File.append("Image" + "," + "Slice" + "," + "%Area" + "," + "IntDen" + "," + "Mean", CellsData);
ResultsPath = resultsDirectory + dataset + "_" + "Bacterial Count Results.csv";
	File.append("Image" + "," + "Type" + "," + "Slice" + "," + "Circ" + "," + "Area" + "," + "AR" + "," + "Skew" + "," + "Kurt" + "," + "X" + "," + "Y" , ResultsPath);
CellsPath = resultsDirectory + dataset + "_" + "Cell Count Results.csv";

//Find the maximum number of frames in all the files, and the total number of ND2 files
nd2Count = 0;
maxFrames = 1;
for(a=0; a<fileList.length; a++){
	if(endsWith(fileList[a], ".nd2")){
		
		//Set the current file
		//Since there are no series in an nd2 file, do NOT set series - causes an error
		Ext.setId(inputDirectory + fileList[a]);
		//Ext.setSeries(1);

		//Count the number of frames from the metadata, and see if it is the largest number of frames
		Ext.getSizeT(frameCount);
		Ext.getSeriesCount(seriesCount);
		if(frameCount > maxFrames) maxFrames = frameCount;
		else if(seriesCount > maxFrames) maxFrames = seriesCount;

		//Check to make sure that all the files have the same number of channels
		if(nd2Count == 0){
			Ext.getSizeC(channelCount);
			nChannels = channelCount;
		}
		else{
			Ext.getSizeC(nChannels);	
		}

		//If the number of channels is inconsistent, exit with error
		if(nChannels != channelCount){
			reportChannels = getBoolean("Error: Inconsistent number of channels.  Would you like to see the number of channels per file?");
			
			//If the # of channels needs to be reported, output into results table
			run("Clear Results");
			if(reportChannels){
				for(b=0; b<fileList.length; b++){
					if(endsWith(fileList[b], ".nd2")){
						Ext.setId(inputDirectory + fileList[b]);
						Ext.getSizeC(nChannels);
						currentRow = nResults;
						setResult("File name", currentRow, fileList[b]);
						setResult("# Channels", currentRow, nChannels);
					}
				}
				updateResults();
			}
			exit("Inconsistent channel count in sample: " + fileList[a]);
		}

		//Add one to the nd2 count
		nd2Count += 1;
	}
}

//Ask the user to specify which channels are the DAPI and PI channels
//Reject any invalid entries
dapiChannel = channelCount + 1;
deadChannel = channelCount + 1;
gfpChannel = channelCount + 1;

while(dapiChannel > channelCount || dapiChannel < 1 || dapiChannel != round(dapiChannel)) dapiChannel = getNumber("Which channel is the DAPI channel: 1-" + channelCount, 1);
while(deadChannel > channelCount || deadChannel == dapiChannel || deadChannel < 1 || deadChannel != round(deadChannel) ) deadChannel = getNumber("Which channel is the Cell Death channel: 1-" + channelCount, 3);
while(gfpChannel > channelCount || gfpChannel == dapiChannel || gfpChannel == deadChannel || gfpChannel < 1 || gfpChannel != round(gfpChannel)) gfpChannel = getNumber("Which channel is the Bacterial channel: 1-" + channelCount, 2);


//Create an array for saving each sample name
sampleNameArray = newArray(nd2Count);

//Open each nd2 file, count the DAPI and PI nuclei, and output the counts to the final image
currentND2 = 0;
for(a=0; a<fileList.length; a++){
	//If the file is an nd2 file, then open it
	if(endsWith(fileList[a], ".nd2")){
		Ext.setId(inputDirectory + fileList[a]);
		Ext.getSizeC(channelCount);

        run("Bio-Formats Importer", "open=[" + inputDirectory + fileList[a] + "] autoscale color_mode=Default concatenate_series open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

        //If opened as series, this fixes the name
        rename(fileList[a]);
		
		//Show current progress
        showProgress(currentND2/nd2Count);
        showStatus("Processing sample " + currentND2+1 + " of " + nd2Count + ": " + fileList[a]);
        print("Processing sample " + currentND2+1 + " of " + nd2Count + ": " + fileList[a]);


//HIGHLY RECOMMENDED: AVOID COMMAS IN IMAGE NAMES AS THIS COMPLICATES THE SAVE FILES
//Save the file names in the input directory as a simple number (e.g. 1-1 or 4-6). 


        //Macro Here:
imageName = getTitle();
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Split Channels");
	run("Set Measurements...", "area mean integrated redirect=None decimal=3");
	run("Options...", "iterations=1 count=1 black do=Nothing");
	run("Clear Results");
	roiManager("reset");


//Set threshold for gfp channel
selectWindow("C" + gfpChannel + "-" + imageName); 
rename("gfp");
	target=round(nSlices*0.2);
	setSlice(target);
	resetMinAndMax();
	run("8-bit");
	run("Measure");
	gfpBackground = getResult("Mean");
	run("Subtract...", "value=gfpBackground stack");
	run("Duplicate...", "title=blur");
		run("Gaussian Blur...", "sigma=100 stack");
		imageCalculator("Subtract stack", "gfp","blur");
		close("blur");
	setSlice(target);
	setOption("BlackBackground", true);	
//Double check which threshold for GFP works best at 8h for your dataset.
		setAutoThreshold("MaxEntropy dark");	//worked well for exsA/pBAD
		run("Clear Results");


//Normalize background intensity of Nuclei using high pass filter	
selectWindow("C" + dapiChannel + "-" + imageName);
rename("nuc");
	setSlice(target);
	resetMinAndMax();
	run("Duplicate...", "title=blur");
		run("Gaussian Blur...", "sigma=100 stack");
		imageCalculator("Subtract stack", "nuc","blur");
		close("blur");
		selectWindow("nuc");
	run("Duplicate...", "title=cells duplicate");
		run("8-bit");
		run("Gaussian Blur...", "sigma=4 stack");
		setAutoThreshold("Triangle dark stack");
		setOption("BlackBackground", true);
		run("Convert to Mask", "method=Triangle background=Default calculate black");
		run("Fill Holes", "stack");
		run("Watershed", "stack");

//Subtract background for clean PI signal 
selectWindow("C" + deadChannel + "-" + imageName);
rename("pi");
	setSlice(target);
	resetMinAndMax();
	run("Duplicate...", "title=blur");
	run("Gaussian Blur...", "sigma=100 stack");
	imageCalculator("Subtract stack", "pi","blur");
	close("blur");
	run("Clear Results");
	run("Measure");
	pi_Background = getResult("Mean", nResults-1);


// Make ROIs of nuclei, count totals
run("Clear Results");
roiManager("reset")
run("Set Measurements...", "mean min integrated stack display redirect=None decimal=3");
selectWindow("cells");
	run("Analyze Particles...", "size=100-Infinity show=Outlines summarize overlay add stack");
	selectWindow("Drawing of cells");
	close("Drawing of cells");
selectWindow("Summary of cells");
	total = Table.getColumn("Count");
	total = Array.concat("Total Nuclei",total);
	close("Summary of cells");


//Delete PI+ Dead cells from ROI index
selectWindow("pi");
for (roi=0;roi<roiManager("Count");roi++){
	roiManager("Select",roi); 
	run("Clear Results");
	run("Measure");
	piMean = getResult("Mean", nResults-1);
	if (piMean >= 10*pi_Background){ 
	roiManager("delete");
	roi = roi-1; //de-iterate the roi variable, otherwise you'll skip some rois
	}
}

//Make dilated cells from live nuclei signal & count
newImage("live_cells", "8-bit black", width, height, frames);
	selectWindow("live_cells");
	roiManager("Select",Array.sequence(roiManager("count")));
	roiManager("fill");
selectWindow("live_cells");
	roiManager("reset");
	run("Analyze Particles...", "size=100-Infinity summarize stack");
selectWindow("Summary of live_cells");
	live_total = Table.getColumn("Count");
	live_total = Array.concat("Live Nuclei",live_total);
	close("Summary of live_cells");

selectWindow("live_cells");
	run("Options...", "iterations=40 count=2 black do=Dilate stack");

//To better measure bacterial area, not only vacuoles
//	run("Options...", "iterations=100 count=2 black do=Dilate stack"); 

	run("Watershed", "stack");
		run("Clear Results");
		roiManager("reset")
selectWindow("live_cells");
	run("Analyze Particles...", "size=100-Infinity show=Outlines summarize overlay add stack");
	selectWindow("Drawing of live_cells");
	close("Drawing of live_cells");
selectWindow("Summary of live_cells");
	live_cells = Table.getColumn("Count");
	live_cells = Array.concat("Live Cells",live_cells);
	close("Summary of live_cells");
	run("Clear Results");


//Remove un-infected cells from Total Cells, Save GFP data from infected cells
	
	selectWindow("gfp");
		run("Set Measurements...", "mean integrated area_fraction stack display redirect=None decimal=3");
	for (roi=0;roi<roiManager("Count");roi++){
		roiManager("Select",roi); 
		run("Measure");
		gfpSlice = getResult("Slice", nResults-1);
		gfpSignal = getResult("%Area",nResults-1);
		gfpInt = getResult("IntDen", nResults-1);
		gfpMean = getResult("Mean", nResults-1);
		if (gfpSignal < 0.01){
			roiManager("delete"); //deletes all cells without thresholded GFP within boundaries.
			roi = roi-1; //de-iterate the roi variable, otherwise you'll skip some rois
			}
			else { File.append(imageName + "," + gfpSlice + "," + gfpSignal + "," + gfpInt + "," + gfpMean, CellsData);	//save data for infected cells
			}
	}

//Create new image with infected cells, count number
newImage("infected_cells", "8-bit black", width, height, frames);
	selectWindow("infected_cells");
	roiManager("Select",Array.sequence(roiManager("count")));
	roiManager("fill");
selectWindow("infected_cells");
	roiManager("reset")
	run("Clear Results");
	run("Set Measurements...", "integrated stack display redirect=None decimal=3");
	run("Analyze Particles...", "size=100-Infinity show=Outlines summarize overlay add stack");
selectWindow("Summary of infected_cells");
	infected = Table.getColumn("Count");
	infected = Array.concat("Infected Cells",infected);
	close("Summary of infected_cells");
	close("Drawing of infected_cells");
	run("Clear Results");


//Measuring Thresholded GFP Circularity, Area, A.R., Skew, Kurt over all time points ONLY WITHIN live cells!
run("Clear Results");
run("Set Measurements...", "area center fit shape skewness kurtosis stack display redirect=None decimal=3");
selectWindow("gfp");
	for (roi=0;roi<roiManager("Count");roi++){
		roiManager("Select",roi); 
		run("Analyze Particles...", "size=0.50-Infinity display slice"); //Measures all thresholded particles within each ROI
	}

	selectWindow("Results");
	for (i = 0; i < nResults(); i++) {
	    l = getResult("Slice", i); 
	    v = getResult("Circ.", i);
	    r = getResult("Area", i);
		ar = getResult("AR", i);	
	    s = getResult("Skew", i);	
	    k = getResult("Kurt", i);	
	    x = getResult("XM", i);
	    y = getResult("YM", i);
	    if ( ((l<6) && (v>0.8) && (r<20)) || ( (l>=6) && (v > 0.8) ) || ( ((0.5 < v) && (v <=0.8)) && ((1.8 <= ar) && (ar < 3.6)) ) ) 
	//Early, circular, small particles OR late circular particles OR somewhat circular particles with Aspect Ratio between 1.8 and 3.6
	    	{File.append(imageName + "," + "Vac" + "," + l + "," + v + "," + r + "," + ar + "," + s + "," + k + "," + x + "," + y, ResultsPath);
	    	
	 	}	else {File.append(imageName + "," + "Spread" + "," + l + "," + v + "," + r + "," + ar + "," + s + "," + k + "," + x + "," + y, ResultsPath);}
			
	}
	
run("Clear Results");


//Place counts of total, live, and infected cells per slice into arrays and save to a new file
selectWindow("infected_cells");
	n = nSlices();


total2 = String.join(total);
live_total2 = String.join(live_total);
live_cells2 = String.join(live_cells);
infected2 = String.join(infected);

total_cell = Array.concat(total2, live_total2, live_cells2, infected2);
total_cell = String.join(total_cell);
File.append(imageName + "," + n + "," + total_cell, CellsPath);
	//This saves a file with a string of measurements extending horizontally. 
	//Each row is one field of view or image.
	//For each row, there will be 4 sets of 'n' measurements, with n = #slices in timelapse. This # is saved as the second column after image name.
		//First 'n' of measurements is total nuclei counted in each frame.  "Total Cells"
		//Second is all nuclei PI signal below threshold. "Live Total"
		//Third is a dilated mask of all live nuclei, representing the intracellular region. "Live Cells"
		//Fourth is all live cells that have any GFP particles over threshold (see line 118). "Infected Cells"

	//Save the name of the sample in the sample name array
	sampleNameArray[currentND2] = replace(fileList[a], ".nd2", "");

    //Add 1 to the ND2 count
	currentND2 += 1;

	close("*");

	}
}

setBatchMode(false);


close("*");
