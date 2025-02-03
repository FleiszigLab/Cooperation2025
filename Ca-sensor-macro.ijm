//Throughout this macro the cell death channel will be referred to as piChannel
//Throughout this macro Ca-Sensor or Bacteria channel will be referred to as gfpChannel

//Close all images that are currently open
close("*");

//Prevent display of images between calls
setBatchMode(true);

//Start the bio-formats extensions function
run("Bio-Formats Macro Extensions");

//Get list of all files
inputDirectory = getDirectory("Choose the input directory");
fileList = getFileList(inputDirectory);

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



//Ask the user to specify which channels are the DAPI andpi channels
//Reject any invalid entries
dapiChannel = channelCount + 1;
piChannel = channelCount + 1;

while(dapiChannel > channelCount || dapiChannel < 1 || dapiChannel != round(dapiChannel)) dapiChannel = getNumber("Which channel is the DAPI channel: 1-" + channelCount, 1);
while(piChannel > channelCount ||piChannel == dapiChannel ||piChannel < 1 ||piChannel != round(piChannel) )piChannel = getNumber("Which channel is the cell death channel: 1-" + channelCount, 3);

//Ask the user if they want to use GFP channel for proximity analysis
useGFP = getBoolean("Do you want to use GFP channel for proximity analysis?");


if (useGFP) {
    // Ask the user to specify which channel is the Ca-sensor channel
    // Reject any invalid entries
    gfpChannel = channelCount + 1;

    while(gfpChannel > channelCount || gfpChannel == dapiChannel || gfpChannel ==piChannel || gfpChannel < 1 || gfpChannel != round(gfpChannel)) gfpChannel = getNumber("Which channel is the Ca-sensor/bacteria channel: 1-" + channelCount, 2);
}

//Create a hyperstack to store the DAPI count and pi count
newImage("DAPI and pi counts per sample.tif", "16-bit black", maxFrames, nd2Count, 2);

//Label each slice accordingly
setSlice(1);
setMetadata("Label", "DAPI counts - Samples by row");
setSlice(2);
setMetadata("Label", "pi counts - Samples by row");

//Create an array for saving each sample name
sampleNameArray = newArray(nd2Count);

//Open each nd2 file, count the DAPI andpi nuclei, and output the counts to the final image
run("Set Measurements...", "area mean stack redirect=None decimal=3");
currentND2 = 0;

for(a=0; a<fileList.length; a++){
	//If the file is an nd2 file, then open it
	if(endsWith(fileList[a], ".nd2")){
		Ext.setId(inputDirectory + fileList[a]);
		Ext.getSizeC(channelCount);

        run("Bio-Formats Importer", "open=[" + inputDirectory + fileList[a] + "] autoscale color_mode=Default concatenate_series open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

        //If opened as series, this fixes the name
        rename(fileList[a]);

        run("Split Channels");

        //Show current progress
        showProgress(currentND2/nd2Count);
        showStatus("Processing sample " + currentND2+1 + " of " + nd2Count + ": " + fileList[a]);

        //Close the channels that are neither DAPI or pi
        if(useGFP){
        for(b=1; b<=channelCount; b++){
            if(b != dapiChannel && b !=piChannel && b != gfpChannel)
             {close("C" + b + "-" + fileList[a]);
            }}}
            else{
           for(b=1; b<=channelCount; b++){
            if(b != dapiChannel && b !=piChannel)
             {close("C" + b + "-" + fileList[a]);}}}
        

        //Set the filenames for the dapi,pi and GFP stacks
        dapiStack = "C" + dapiChannel + "-" + fileList[a];
       piStack = "C" +piChannel + "-" + fileList[a];
        if (useGFP){
        gfpStack = "C" + gfpChannel + "-" + fileList[a];
        }

		//Normalize background intensity using high pass filter
		selectWindow(dapiStack);
		run("Duplicate...", "title=blur duplicate");
		run("Gaussian Blur...", "sigma=100 stack");
		imageCalculator("Subtract stack", dapiStack,"blur");
		close("blur");
		selectWindow(piStack);
		run("Duplicate...", "title=blur duplicate");
		run("Gaussian Blur...", "sigma=100 stack");
		imageCalculator("Subtract stack",piStack,"blur");
		close("blur");
		if(useGFP){
		selectWindow(gfpStack);
		run("Measure");
		gfpBackground = getResult("Mean")/2;
		run("Subtract...","value=gfpBackground");
		}
		
		//GFP mask
		if(useGFP){
selectWindow(gfpStack);
run("Duplicate...", "title=gfpMask duplicate");
run("8-bit");
setOption("BlackBackground", true);
setAutoThreshold("Triangle dark stack");
run("Convert to Mask", "method=Traingle background=Default black");
run("Fill Holes", "stack");
		}
		
		//Segment the DAPI nuclei using a watershed
selectWindow(dapiStack);
run("Duplicate...", "title=dapiMask duplicate");
run("8-bit");
run("Gaussian Blur...", "sigma=4 stack");
setAutoThreshold("Triangle dark stack");
setOption("BlackBackground", true);
run("Convert to Mask", "method=Triangle background=Default black");
run("Fill Holes", "stack");
if(useGFP){
    run("Options...", "iterations=15 count=1 black do=Nothing");
    run("Dilate", "stack");
}
run("Watershed", "stack");
close(dapiStack);

// Set input parameters
if(useGFP){
roiManager("reset")
run("Set Measurements...", "area mean integrated stack redirect=None decimal=3");
run("Analyze Particles...", "  show=Outlines overlay add stack");
selectWindow("Drawing of dapiMask");
close("Drawing of dapiMask");
selectWindow("gfpMask");
roiManager("multi-measure measure_all");

//NewImage
dapiMaskWidth = getWidth();
dapiMaskHeight = getHeight();
dapiMaskSlices = nSlices;
newImage("FilteredROIs", "8-bit black", dapiMaskWidth, dapiMaskHeight, dapiMaskSlices);

//Check ROIs
selectWindow("gfpMask");
for (roi=0;roi<roiManager("Count");roi++){//for each roi
roiManager("Select",roi); //select it and measure it
run("Measure");
gfpSignal = getResult("IntDen",nResults-1);
if (gfpSignal < 5){ //if it's below threshold (you need to replace threshold with your own value
roiManager("delete"); //delete it
roi = roi-1; //de-iterate the roi variable, otherwise you'll skip some rois
}
}
//Fill new image
selectWindow("FilteredROIs");
roiManager("Select",Array.sequence(roiManager("count")));
roiManager("fill");
close("dapiMask");
selectImage("FilteredROIs");
run("Watershed", "stack");
rename("dapiMask");
close("gfpMask");
    }
   
   //Create a mask to count DAPI and pi objects
		selectWindow("dapiMask");
		run("Divide...", "value=255 stack");

		//Add 1 to the pi image and then multiply by mask
		//This way, even DAPI nuclei with no pi will still have an
		//intensity of 1
		selectWindow(piStack);
		run("8-bit");
		run("Add...", "value=1 stack");
		imageCalculator("Multiply stack",piStack, "dapiMask");
		close("dapiMask");
				
		//Measure the pi intensity of each nucleus
		setThreshold(1, 255);
		run("Analyze Particles...", "size=50-5000 circularity=0.20-1.00 show=Nothing display clear stack");

		/// Find the first frame with nuclei
		firstFrameWithNuclei = -1;
		for (b = 0; b < nResults; b++) {
    	frame = getResult("Slice", b);
    	if (frame != firstFrameWithNuclei && frame > firstFrameWithNuclei) {
        firstFrameWithNuclei = frame;
        break;
    }
}

		if (firstFrameWithNuclei != -1) {
    // Count the number of nuclei in the first frame with nuclei
    nucleiCount = 0;
    for (b = 0; b < nResults; b++) {
        frame = getResult("Slice", b);
        if (frame == firstFrameWithNuclei) {
            nucleiCount += 1;
        } else if (frame > firstFrameWithNuclei) {
            break;
        }
    }
  

    // Build an array for saving the nuclei intensity in the first frame with nuclei
    frame1Array = newArray(nucleiCount);
    idx = 0;
    for (b = 0; b < nResults; b++) {
        frame = getResult("Slice", b);
        if (frame == firstFrameWithNuclei) {
            frame1Array[idx++] = getResult("Mean", b);
        } else if (frame > firstFrameWithNuclei) {
            break;
        }
    }


    // Find the median nucleus intensity
    Array.sort(frame1Array);
    
    if (frame1Array.length == 1) {
    backgroundInt = frame1Array[0];
	} 
	else {
    backgroundInt = frame1Array[round(frame1Array.length / 2)];
	}


		//Define positive nuclei as nuclei with an intensity that is 10x initial median
		piThreshold = 10*backgroundInt;

		//Count the number of nuclei and pi positive nuclei in each frame
		totalCountArray = newArray(maxFrames);
		PICountArray = newArray(maxFrames);

		for(b=0; b<nResults; b++){
			currentSlice = getResult("Slice", b);
			PIIntensity = getResult("Mean", b);
			totalCountArray[currentSlice-1] += 1;
			if(piIntensity >piThreshold)piCountArray[currentSlice-1] += 1;
		}

		//Export the counts to the count image
		selectWindow("DAPI and pi counts per sample.tif");
		setSlice(1);
		for(b=0; b<PICountArray.length; b++) setPixel(b, currentND2, totalCountArray[b]);
		setSlice(2);
		for(b=0; b<PICountArray.length; b++) setPixel(b, currentND2,piCountArray[b]);			

		//Save the name of the sample in the sample name array
		sampleNameArray[currentND2] = replace(fileList[a], ".nd2", "");
		
		//Add 1 to the ND2 count
		currentND2 += 1;

		//close the data stack
		close(piStack); 
	}
		
}
}

//Convert the count image into a set of spreadsheets
run("Clear Results");
if(useGFP){
	File.makeDirectory(inputDirectory + "Macro Results Infected cells");
	resultsDirectory = inputDirectory + File.separator + "Macro Results Infected cells" + File.separator;

//Export the DAPI count results
for(a=0; a<nd2Count; a++){
	setResult("Sample", a, sampleNameArray[a]);
	selectWindow("DAPI and pi counts per sample.tif");
	setSlice(1);
	for(b=0; b<maxFrames; b++) setResult("Frame " + b, a, getPixel(b,a));
}

updateResults();
saveAs("Results", resultsDirectory + "DAPI nuclei counts per frame.csv");

run("Clear Results");
//Export the DAPI count results
for(a=0; a<nd2Count; a++){
	setResult("Sample", a, sampleNameArray[a]);
	selectWindow("DAPI and pi counts per sample.tif");
	setSlice(2);
	for(b=0; b<maxFrames; b++) setResult("Frame " + b, a, getPixel(b,a));
}
updateResults();
saveAs("Results", resultsDirectory + "pi nuclei counts per frame.csv");

//Calculate the ratio of pi nuclei to total nuclei
selectWindow("DAPI and pi counts per sample.tif");
run("Stack to Images");
imageCalculator("Divide create 32-bit", "pi counts - Samples by row","DAPI counts - Samples by row");
selectWindow("pi counts - Samples by row");
saveAs("Tiff", resultsDirectory + "pi counts - Samples by row.tif");
close("pi counts - Samples by row.tif");
selectWindow("DAPI counts - Samples by row");
saveAs("Tiff", resultsDirectory + "DAPI counts - Samples by row.tif");
close("DAPI counts - Samples by row.tif");
selectWindow("Result of pi counts - Samples by row");
saveAs("Tiff", resultsDirectory + "pi ratio image.tif");

//Export the ratio results to a spreadsheet
run("Clear Results");
//Export the pi ratio results
for(a=0; a<nd2Count; a++){
	setResult("Sample", a, sampleNameArray[a]);
	selectWindow("pi ratio image.tif");
	for(b=0; b<maxFrames; b++) setResult("Frame " + b, a, getPixel(b,a));
}

updateResults();
saveAs("Results", resultsDirectory + "pi ratio per frame.csv");
}
	

else { 
	File.makeDirectory(inputDirectory + "Macro Results Total");
	resultsDirectory = inputDirectory + File.separator + "Macro Results Total" + File.separator;

//Export the DAPI count results
for(a=0; a<nd2Count; a++){
	setResult("Sample", a, sampleNameArray[a]);
	selectWindow("DAPI and pi counts per sample.tif");
	setSlice(1);
	for(b=0; b<maxFrames; b++) setResult("Frame " + b, a, getPixel(b,a));
}

updateResults();
saveAs("Results", resultsDirectory +"DAPI nuclei counts per frame.csv");

run("Clear Results");
//Export the DAPI count results
for(a=0; a<nd2Count; a++){
	setResult("Sample", a, sampleNameArray[a]);
	selectWindow("DAPI andpi counts per sample.tif");
	setSlice(2);
	for(b=0; b<maxFrames; b++) setResult("Frame " + b, a, getPixel(b,a));
}
updateResults();
saveAs("Results", resultsDirectory+ "pi nuclei counts per frame.csv");

//Calculate the ratio of pi nuclei to total nuclei
selectWindow("DAPI and pi counts per sample.tif");
run("Stack to Images");
imageCalculator("Divide create 32-bit", "pi counts - Samples by row","DAPI counts - Samples by row");
selectWindow("pi counts - Samples by row");
saveAs("Tiff", resultsDirectory+ "pi counts - Samples by row.tif");
close("pi counts - Samples by row.tif");
selectWindow("DAPI counts - Samples by row");
saveAs("Tiff", resultsDirectory + "DAPI counts - Samples by row.tif");
close("DAPI counts - Samples by row.tif");
selectWindow("Result of pi counts - Samples by row");
saveAs("Tiff", resultsDirectory+ "pi ratio image.tif");

//Export the ratio results to a spreadsheet
run("Clear Results");
//Export the pi ratio results
for(a=0; a<nd2Count; a++){
	setResult("Sample", a, sampleNameArray[a]);
	selectWindow("pi ratio image.tif");
	for(b=0; b<maxFrames; b++) setResult("Frame " + b, a, getPixel(b,a));
}

updateResults();
saveAs("Results", resultsDirectory + "pi ratio per frame.csv");
}
selectWindow("Results");
run("Close");
setBatchMode(false);
close("*");
//Array.print(totalCountArray);
//Array.print(piCountArray);
//Array.print(sampleNameArray);