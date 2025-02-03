//Open - split channels - background subtract - merge channels - select GFP - find maxima - go to RFP - Measure - save Mean value - next image

//Close all images that are currently open
close("*");

//Prevent display of images between calls
setBatchMode(true);

//Start the bio-formats extensions function
run("Bio-Formats Macro Extensions");

//Get list of all files and set output directory
inputDirectory = getDirectory("Choose the input directory");
fileList = getFileList(inputDirectory);
	File.makeDirectory(inputDirectory + "Macro Results");
	resultsDirectory = inputDirectory + File.separator + "Macro Results" + File.separator;

//Find the maximum number of frames in all the files, and the total number of ND2 files
nd2Count = 0;
for(a=0; a<fileList.length; a++){
	if(endsWith(fileList[a], ".nd2")){
		
		//Set the current file
		//Since there are no series in an nd2 file, do NOT set series - causes an error
		Ext.setId(inputDirectory + fileList[a]);
		//Ext.setSeries(1);

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
			
			exit("Inconsistent channel count in sample: " + fileList[a]);
		}

		//Add one to the nd2 count
		nd2Count += 1;
	}
}

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

//Define Results Outputs
//resultsDirectory = getDirectory("downloads");
Results = resultsDirectory + "GFP on PI Results.csv";

        //Macro Here:
imageName = getTitle();
	run("Split Channels");
	run("Set Measurements...", "mean median redirect=None decimal=3");
	run("Options...", "iterations=1 count=1 black do=Nothing");
	run("Clear Results");
	roiManager("reset");

selectWindow("C1-" + imageName);
	rename("gfp");
	resetMinAndMax();
	run("Measure");
	gfpBackground = getResult("Median");
	run("Subtract...", "value=gfpBackground");
	run("Clear Results");

selectWindow("C2-" + imageName);
	rename("rfp");
	resetMinAndMax();
	run("Measure");
	rfpBackground = getResult("Median");
	run("Subtract...", "value=rfpBackground");
	run("Clear Results");

run("Merge Channels...", "c1=rfp c2=gfp create");
	Stack.setDisplayMode("color");
	Stack.setChannel(1);
	run("Find Maxima...", "prominence=10 output=[Point Selection]");
	Stack.setChannel(2);
	run("Measure");
	selectWindow("Results");
	gfp_data = Table.getColumn("Mean");
	gfp_data = Array.concat(imageName, gfp_data);
	run("Clear Results");

gfp_data2 = String.join(gfp_data);
File.append(gfp_data2, Results);




    //Add 1 to the ND2 count
	currentND2 += 1;

	close("*");

	}
}

setBatchMode(false);


close("*");