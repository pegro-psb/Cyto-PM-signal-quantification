//@Boolean (label="Delete original ROIs", value=false, description="Delete the original Cell selection ROIs from the ROI Manager") deleteOriginalROis
//@Float (label="Median Radius Size", value=0.5) medianRadiusSize
//@Integer (label="top X percentile", value=10) xPercentile
//@Integer (label="AutoLocal Threshold radius", value=15) autoLocaThresholdRadius
//Boolean (Label="Record intensity statistic in the log", value=false, description="Record statistic in the log") doRecordStatistic
//String(label="Threshold Method for Tissue Detection",choices={'Phansalkar','Mean','Otsu','Niblack', 'Sauvola'}) thresholdMethod


//From a plant fluorescent microscopy images with cell wall in channel 1 and GFP in channel 2
//from a collection of manually annotated cells stored in the ROI Manager, detect the membrane and the cytoplasm
//of each annotated cell and measure the top 10 percent mean intensty of the channel 2 in the detected membrane and cytoplasm
//output a result table with the measurement, 1 row per cell and store the detected membrane and cytoplasm in the ROI Manager


//2020/12/10
/**
Requirements:
Fiji
Auto Local Threshold Plugin

@Author: Benjamin Pavie - benjamin.pavie@vib.be VIB BioImaging Core
Please acknowledge us if you use this script in a publication.
**/


thresholdMethod="Phansalkar";
//Clear the log
print("\\Clear");
setOption("ExpandableArrays", true);//needed by getFilesListByExtension
setBatchMode(false);
originalID=getImageID();
run("Colors...", "foreground=white background=black selection=yellow");
setOption("BlackBackground", true);

//Store the ID of the open image
originalID=getImageID();
//Get the cell number from the ROIManager
cellNumber=RoiManager.size;
cellSegmentationSucceedArray=newArray(cellNumber);
numberOfSuccessfullCell=0;

for(i=0;i<cellNumber;i++)
{  
  selectImage(originalID);
  RoiManager.select(i);
  getSelectionBounds(x, y, width, height);  
  run("Duplicate...", "duplicate");
  singleCellID=getImageID();
  //Get the first channel
  run("Duplicate...", "duplicate channels=1");
  singleCellChannel1ID=getImageID();
  // 1- Detect the cytoplasm
  //Filter the membrane
  run("Unsharp Mask...", "radius=5 mask=0.60");
  //threshold it using a local threshold, more sensitive and get a mask
  run("Auto Local Threshold", "method="+thresholdMethod+" radius="+autoLocaThresholdRadius+" parameter_1=0 parameter_2=0 white");
  //Apply a median filter to smooth it
  run("Median...", "radius="+medianRadiusSize); 
  //CLear everything outisde the original selection
  run("Clear Outside");
  //Deselect everything
  run("Select None");
  run("Invert");
  run("Grays");
  run("Invert LUT");
  //Detect the cytoplasm and add it to the ROI Manager
  setOption("BlackBackground", false);
  run("Analyze Particles...", "size=20-Infinity pixel exclude add"); //Add it as index 0 to the ROI Manager
  
  if(RoiManager.size!=(cellNumber+numberOfSuccessfullCell*2+1))
  {
    cellSegmentationSucceedArray[i]=false;
  }
  else
  {
    //Clear the leftover particle inside the cytoplasm  
    RoiManager.select(cellNumber+numberOfSuccessfullCell*2);
    run("Clear", "slice");
    // 2- Detect the membrane
    //Deselect everything
    run("Select None");
    run("Invert");
    run("Create Selection");
    roiManager("Add"); //Add it as index 1 to the ROI Manager
    //reset the location of the detection to match the image
    RoiManager.select(cellNumber+numberOfSuccessfullCell*2);
    getSelectionBounds(x1, y1, width1, height1);
    setSelectionLocation(x1+x, y1+y);
    RoiManager.select(cellNumber+numberOfSuccessfullCell*2+1);
    getSelectionBounds(x2, y2, width2, height2);
    setSelectionLocation(x2+x, y2+y);
    
    cellSegmentationSucceedArray[i]=true;
    numberOfSuccessfullCell=numberOfSuccessfullCell+1;
  }
  //Close the chanel 1 cell image
  selectImage(singleCellChannel1ID);
  close();
  //Close the chanel 2 cell image
  selectImage(singleCellID);
  close();
}
//Select the 2nd channel
setSlice(2);
run("Clear Results");
numberOfSuccessfullCell=0;
for(i=0;i<cellNumber;i++)
{
  row = nResults;
  setResult("Cell Number ", row, (i+1));
  if(cellSegmentationSucceedArray[i]==true)
  {
    //Select the cytoplasm
    RoiManager.select(cellNumber+numberOfSuccessfullCell*2);
    getRawStatistics(nPixels, mean, min, max, std, histogram);
    meantopXPercent= getMeanTopXPercent(histogram,nPixels,xPercentile);
    setResult("Cyto Top "+xPercentile+"% Mean", row, meantopXPercent);
    //measure the intensity on channel 2
    //Select the membrane
    RoiManager.select(cellNumber+numberOfSuccessfullCell*2+1);
    getRawStatistics(nPixels, mean, min, max, std, histogram);
    meantopXPercent= getMeanTopXPercent(histogram,nPixels,xPercentile);
    setResult("Membrane Top "+xPercentile+"% Mean", row, meantopXPercent);
    numberOfSuccessfullCell=numberOfSuccessfullCell+1;
  }
  else
  {
    setResult("Cyto Top "+xPercentile+"% Mean", row, -1);
    setResult("Membrane Top "+xPercentile+"% Mean", row, -1);
  }
}
if(deleteOriginalROis)
{
  //optionally, clear the original ROIs to keep only the detection
  for(i=cellNumber-1;i>=0;i--)
  {
    RoiManager.select(i);
    roiManager("delete");
  }
}
setBatchMode(false);
updateResults();

function getMeanTopXPercent(histogram,nPixels,xPercentile)
{
  //Sort the histogram from Max intensity to bottom, each row of the list is an intensity 
  //and value is the number of pixel having this intensity
  Array.reverse(histogram);
  //Get the number of pixel we need to sum the intensity to reach the xPercentile,
  //starting from max intensity pixels
  tenPercentIndex=round(nPixels*xPercentile/100);
  nrPixelSummed=0;
  totalIntensity=0;
  done=false;
  for(i=0;i<histogram.length;i++)
  {
    if(done==false)
    {
      pixelNumber=histogram[i];
      if(nrPixelSummed+pixelNumber<tenPercentIndex)
      {
        totalIntensity=totalIntensity+pixelNumber*(histogram.length-i-1);//pixelNumber*(histogram.length-i-1);
        nrPixelSummed=nrPixelSummed+pixelNumber;
      }
      else
      {        
        //print("pixelNumber:"+tenPercentIndex-nrPixelSummed);
        totalIntensity=totalIntensity+(tenPercentIndex-nrPixelSummed)*(histogram.length-i-1);
        nrPixelSummed=nrPixelSummed+(tenPercentIndex-nrPixelSummed);
        done = true; // break 
      }
    }
  }
  mean = totalIntensity/tenPercentIndex;
  //print("mean:"+mean);

  return mean;
}
