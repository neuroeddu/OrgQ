#
#		OrgQ		- ImageJ Macro written in  Python
#
#		Input 		- Folders containing split channel images for Organoids
#					- The filenames of the images must end in ch00, ch01, ch02, ch03 and have merged in the names ex. "Merged_ch00"
#                   - Folder must also contained the merged image file ending with Merged.tif
#					- Optional thresholding can be loaded as well

#		Output		- CSV file containing nuclei counts and marker colocalization data for each image
#					- Intensity per area
#					- Count number of pixels for each marker

#		Written by: 						Eddie Cai, Rhalena A. Thomas & Harris Nami
#		Incorporating analysis ideas from Vincent Soubanier and Valerio Piscopo


import os, sys, math, csv, datetime
from ij import IJ, Prefs, ImagePlus
from ij.io import DirectoryChooser
from ij.io import OpenDialog
from ij.measure import ResultsTable
from ij.measure import Measurements
from ij.process import ImageProcessor
from ij.process import ImageConverter
from ij.plugin.frame import RoiManager
from ij.plugin.filter import ParticleAnalyzer
from ij.gui import GenericDialog
from ij.gui import WaitForUserDialog
from ij.plugin.filter import ThresholdToSelection
from ij.plugin import ImageCalculator
from ij.WindowManager import getCurrentImage
import xml.etree.ElementTree as ET


# To enable displayImages mode (such as for testing thresholds), make displayImages = True
displayImages = True

#Enable using the wand tool to manually select the organoid ROI in cases where auto-threshold does not work
enableWand = True

# Function to get the markers needed with a generic dialog for each subfolder, as well as the name of the output for that subfolder
def getChannels(subFolder):
    gd = GenericDialog("Channel Options")

    gd.addMessage("Name the markers associated with this directory:")
    gd.addMessage(inputDirectory + subFolder)
    gd.addMessage("(Leave empty to ignore)")
    gd.addMessage("")
    gd.addStringField("Channel ch00:", "Dapi")
    gd.addStringField("Channel ch01:", "pSYN")
    gd.addStringField("Channel ch02:", "MAP2")
    gd.addStringField("Channel ch03:", "SYN")
    gd.addMessage("")
    gd.addStringField("What would you like the output file to be named:", subFolder)

    gd.showDialog()

    channelNames = []

    channelNames.append([gd.getNextString(), 0])
    channelNames.append([gd.getNextString(), 1])
    channelNames.append([gd.getNextString(), 2])
    channelNames.append([gd.getNextString(), 3])
    outputName = gd.getNextString()

    channels = []
    for i, v in enumerate(channelNames):
        if v[0] != "":
            channels.append(v)

    if gd.wasCanceled():
        print
        "User canceled dialog!"
        return

    return channels, outputName


# Function to get the thresholds.

def getThresholds():
    thresholds = {}

    gd = GenericDialog("Threshold options")
    gd.addChoice("How would you like to set your thresholds?", ["default", "use threshold csv file"], "default")
    gd.showDialog()

    choice = gd.getNextChoice()
    log.write("Option: " + choice + "\n")

    if choice == "use threshold csv file":
        path = OpenDialog("Open the thresholds csv file")
        log.write("File used: " + path.getPath() + "\n")
        with open(path.getPath()) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                thresholds = row
    return thresholds


def rreplace(s, old, new):
    return (s[::-1].replace(old[::-1], new[::-1], 1))[::-1]


############# Main loop, will run for every image. ##############

def process(subFolder, outputDirectory, filename):
    #IJ.close()
    imp = IJ.openImage(inputDirectory + subFolder + '/' + rreplace(filename, "_ch00.tif", ".tif"))
    imp.show()


    # Finds the pixel length in microns from the xml metadata file
    file_list = [file for file in os.listdir(inputDirectory + subFolder) if file.endswith('.xml')]
    if len(file_list) > 0:
        xml = os.path.join(inputDirectory + subFolder, file_list[0])
        element_tree = ET.parse(xml)
        root = element_tree.getroot()
        for dimensions in root.iter('DimensionDescription'):
            num_pixels = int(dimensions.attrib['NumberOfElements'])
            if dimensions.attrib['Unit'] == "m":
                length = float(dimensions.attrib['Length']) * 1000000
            else:
                length = float(dimensions.attrib['Length'])
        pixel_length = length / num_pixels
    else:
        pixel_length = 0.877017

    log.write("Pixel Length:" + str(pixel_length) + "\n")

    IJ.run(imp, "Properties...",
           "channels=1 slices=1 frames=1 unit=um pixel_width=" + str(pixel_length) + " pixel_height=" + str(pixel_length) + " voxel_depth=25400.0508001")
    ic = ImageConverter(imp);
    ic.convertToGray8();
    #IJ.setThreshold(imp, 2, 255)


    # If wand tool is enabled, then this will prompt that to be used
    if enableWand:
        # Call threshold function to adjust threshold and select Organoid ROI
        IJ.run("Threshold...")
        WaitForUserDialog("Adjust Threshold to create mask").show()
        IJ.setTool("Wand")
        WaitForUserDialog("Click on Organoid Area for it to be selected. Best selection will be at the edge of the organoid to get entire organoid shape.").show()
        IJ.run("Clear Outside")

    if not enableWand:
        IJ.setAutoThreshold(imp, "Mean dark no-reset")
        IJ.run(imp, "Convert to Mask", "")
        IJ.run(imp, "Analyze Particles...", "size=100000-Infinity add select")
        rm = RoiManager.getInstance()
        imp = getCurrentImage()
        rm.select(imp, 0)
        IJ.setBackgroundColor(0, 0, 0)
        IJ.run(imp, "Clear Outside", "")

    IJ.run(imp, "Convert to Mask", "")
    IJ.run(imp, "Remove Outliers...", "radius=5" + " threshold=50" + " which=Dark")
    IJ.run(imp, "Remove Outliers...", "radius=5" + " threshold=50" + " which=Bright")

    # #Save the mask and open it
    IJ.saveAs("tiff", inputDirectory + '/mask')
    mask = IJ.openImage(inputDirectory + '/mask.tif')

    if enableWand:
        #Select ROI again to add it to the the ROI manager so that intensities and area is saved
        #IJ.run("Threshold...")
        IJ.setTool("Wand")
        WaitForUserDialog("Select Organoid area again for it to register within the ROI manager").show()
        rm = RoiManager()
        boundroi = ThresholdToSelection.run(mask)
        rm.addRoi(boundroi)


    if not displayImages:
        imp.changes = False
        imp.close()

    images = [None] * 5
    intensities = [None] * 5
    blobsarea = [None] * 5
    blobsnuclei = [None] * 5
    bigAreas = [None] * 5

    imp.close()

    #Loop to open all the channel images
    for chan in channels:
        v, x = chan
        images[x] = IJ.openImage(
            inputDirectory + subFolder + '/' + rreplace(filename, "_ch00.tif", "_ch0" + str(x) + ".tif"))


        # Apply Mask on all the images and save them into an array
        apply_mask = ImageCalculator()
        images[x] = apply_mask.run("Multiply create 32 bit", mask, images[x])
        ic = ImageConverter(images[x])
        ic.convertToGray8()
        imp = images[x]

        # Calculate the intensities for each channel as well as the organoid area
        for roi in rm.getRoisAsArray():
            imp.setRoi(roi)
            stats = imp.getStatistics(Measurements.MEAN | Measurements.AREA)
            intensities[x] = stats.mean
            bigAreas[x] = stats.area

    rm.close()

    # Opens the ch00 image and sets default properties
    apply_mask = ImageCalculator()
    imp = IJ.openImage(inputDirectory + subFolder + '/' + filename)
    imp = apply_mask.run("Multiply create 32 bit", mask, imp)
    IJ.run(imp, "Properties...",
           "channels=1 slices=1 frames=1 unit=um pixel_width=" + str(pixel_length) + " pixel_height=" + str(pixel_length) + " voxel_depth=25400.0508001")


    # Sets the threshold and watersheds. for more details on image processing, see https://imagej.nih.gov/ij/developer/api/ij/process/ImageProcessor.html

    ic = ImageConverter(imp);
    ic.convertToGray8();

    IJ.run(imp, "Remove Outliers...", "radius=2" + " threshold=50" + " which=Dark")

    IJ.run(imp, "Gaussian Blur...", "sigma=" + str(blur))

    IJ.setThreshold(imp, lowerBounds[0], 255)

    if displayImages:
        imp.show()
    IJ.run(imp, "Convert to Mask", "")
    IJ.run(imp, "Watershed", "")

    if not displayImages:
        imp.changes = False
        imp.close()

    # Counts and measures the area of particles and adds them to a table called areas. Also adds them to the ROI manager

    table = ResultsTable()
    roim = RoiManager(True)
    ParticleAnalyzer.setRoiManager(roim);
    pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, Measurements.AREA, table, 15, 9999999999999999, 0.2, 1.0)
    pa.setHideOutputImage(True)
    # imp = impM

    # imp.getProcessor().invert()
    pa.analyze(imp)

    areas = table.getColumn(0)

    # This loop goes through the remaining channels for the other markers, by replacing the ch00 at the end with its corresponding channel
    # It will save all the area fractions into a 2d array called areaFractionsArray

    areaFractionsArray = [None] * 5
    for chan in channels:
        v, x = chan
        # Opens each image and thresholds

        imp = images[x]

        IJ.run(imp, "Properties...",
               "channels=1 slices=1 frames=1 unit=um pixel_width=" + str(pixel_length) + " pixel_height=" + str(pixel_length) + " voxel_depth=25400.0508001")

        ic = ImageConverter(imp);
        ic.convertToGray8();
        IJ.setThreshold(imp, lowerBounds[x], 255)


        if displayImages:
            imp.show()
            WaitForUserDialog("Title", "Adjust Threshold for Marker " + v).show()

        IJ.run(imp, "Convert to Mask", "")

        # Measures the area fraction of the new image for each ROI from the ROI manager.
        areaFractions = []
        for roi in roim.getRoisAsArray():
            imp.setRoi(roi)
            stats = imp.getStatistics(Measurements.AREA_FRACTION)
            areaFractions.append(stats.areaFraction)

        # Saves the results in areaFractionArray

        areaFractionsArray[x] = areaFractions

    roim.close()

    for chan in channels:
        v, x = chan

        imp = images[x]
        imp.deleteRoi()
        roim = RoiManager(True)
        ParticleAnalyzer.setRoiManager(roim);
        pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, Measurements.AREA, table, 15, 9999999999999999, 0.2, 1.0)
        pa.analyze(imp)

        blobs = []
        for roi in roim.getRoisAsArray():
            imp.setRoi(roi)
            stats = imp.getStatistics(Measurements.AREA)
            blobs.append(stats.area)

        blobsarea[x] = sum(blobs) #take this out and use intial mask tissue area from the beginning
        blobsnuclei[x] = len(blobs)


        if not displayImages:
            imp.changes = False
            imp.close()
        roim.reset()
        roim.close()

        imp.close()

    # Creates the summary dictionary which will correspond to a single row in the output csv, with each key being a column

    summary = {}

    summary['Image'] = filename
    summary['Directory'] = subFolder

    # Adds usual columns

    summary['size-average'] = 0
    summary['#nuclei'] = 0
    summary['all-negative'] = 0

    summary['too-big-(>' + str(tooBigThreshold) + ')'] = 0
    summary['too-small-(<' + str(tooSmallThreshold) + ')'] = 0

    # Creates the fieldnames variable needed to create the csv file at the end.

    fieldnames = ['Name', 'Directory', 'Image', 'size-average', 'too-big-(>' + str(tooBigThreshold) + ')',
                  'too-small-(<' + str(tooSmallThreshold) + ')', '#nuclei', 'all-negative']

    # Adds the columns for each individual marker (ignoring Dapi since it was used to count nuclei)

    summary["organoid-area"] = bigAreas[x]
    fieldnames.append("organoid-area")

    for chan in channels:
        v, x = chan
        summary[v + "-positive"] = 0
        fieldnames.append(v + "-positive")

        summary[v + "-intensity"] = intensities[x]
        fieldnames.append(v + "-intensity")

        summary[v + "-blobsarea"] = blobsarea[x]
        fieldnames.append(v + "-blobsarea")

        summary[v + "-blobsnuclei"] = blobsnuclei[x]
        fieldnames.append(v + "-blobsnuclei")

    # Adds the column for colocalization between first and second marker

    if len(channels) > 2:
        summary[channels[1][0] + '-' + channels[2][0] + '-positive'] = 0
        fieldnames.append(channels[1][0] + '-' + channels[2][0] + '-positive')

    # Adds the columns for colocalization between all three markers

    if len(channels) > 3:
        summary[channels[1][0] + '-' + channels[3][0] + '-positive'] = 0
        summary[channels[2][0] + '-' + channels[3][0] + '-positive'] = 0
        summary[channels[1][0] + '-' + channels[2][0] + '-' + channels[3][0] + '-positive'] = 0

        fieldnames.append(channels[1][0] + '-' + channels[3][0] + '-positive')
        fieldnames.append(channels[2][0] + '-' + channels[3][0] + '-positive')
        fieldnames.append(channels[1][0] + '-' + channels[2][0] + '-' + channels[3][0] + '-positive')

    # Loops through each particle and adds it to each field that it is True for.

    areaCounter = 0
    for z, area in enumerate(areas):

        log.write(str(area))
        log.write("\n")

        if area > tooBigThreshold:
            summary['too-big-(>' + str(tooBigThreshold) + ')'] += 1
        elif area < tooSmallThreshold:
            summary['too-small-(<' + str(tooSmallThreshold) + ')'] += 1
        else:

            summary['#nuclei'] += 1
            areaCounter += area

            temp = 0
            for chan in channels:
                v, x = chan
                if areaFractionsArray[x][z] > areaFractionThreshold[
                    0]:  # theres an error here im not sure why. i remember fixing it before
                    summary[chan[0] + '-positive'] += 1
                    if x != 0:
                        temp += 1

            if temp == 0:
                summary['all-negative'] += 1

            if len(channels) > 2:
                if areaFractionsArray[1][z] > areaFractionThreshold[1]:
                    if areaFractionsArray[2][z] > areaFractionThreshold[2]:
                        summary[channels[1][0] + '-' + channels[2][0] + '-positive'] += 1

            if len(channels) > 3:
                if areaFractionsArray[1][z] > areaFractionThreshold[1]:
                    if areaFractionsArray[3][z] > areaFractionThreshold[3]:
                        summary[channels[1][0] + '-' + channels[3][0] + '-positive'] += 1
                if areaFractionsArray[2][z] > areaFractionThreshold[2]:
                    if areaFractionsArray[3][z] > areaFractionThreshold[3]:
                        summary[channels[2][0] + '-' + channels[3][0] + '-positive'] += 1
                        if areaFractionsArray[1][z] > areaFractionThreshold[1]:
                            summary[channels[1][0] + '-' + channels[2][0] + '-' + channels[3][0] + '-positive'] += 1

    # Calculate the average of the particles sizes

    if float(summary['#nuclei']) > 0:
        summary['size-average'] = round(areaCounter / summary['#nuclei'], 2)

    # Opens and appends one line on the final csv file for the subfolder (remember that this is still inside the loop that goes through each image)

    with open(outputDirectory + "/" + outputName + ".csv", 'a') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction='ignore', lineterminator='\n')
        if os.path.getsize(outputDirectory + "/" + outputName + ".csv") < 1:
            writer.writeheader()
        writer.writerow(summary)

    IJ.run(imp, "Close All", "")

########################## code begins running here ##############################


# Get input and output directories

dc = DirectoryChooser("Choose an input directory")
inputDirectory = dc.getDirectory()

dc = DirectoryChooser("Choose an output directory")
outputDirectory = dc.getDirectory()

# Opens log file

with open(outputDirectory + "log.txt", "w") as log:
    log.write("log: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    log.write("\n")

    log.write("________________________\n")
    log.write("Input directory selected: " + inputDirectory + "\n")

    log.write("________________________\n")
    log.write("Output directory selected: " + outputDirectory + "\n")

    # Finds all the subfolders in the main directory

    directories = []

    for subFolder in os.listdir(inputDirectory):
        if os.path.isdir(inputDirectory + subFolder):
            directories.append(subFolder)
            print("subfolder found")

    # A few default options

    areaFractionThreshold = [0.1, 0.1, 0.1, 0.1, 0.1]  # you can change these
    tooSmallThreshold = 50
    tooBigThreshold = 500
    blur = 1

    log.write("________________________\n")
    log.write("Default calculation thresholds: \n")
    log.write("	areaFractionThreshold:" + str(areaFractionThreshold) + "\n")
    log.write("	tooSmallThreshold:" + str(tooSmallThreshold) + "\n")
    log.write("	tooBigThreshold:" + str(tooBigThreshold) + "\n")

    # Get options from user. (see functions written on top)

    log.write("________________________\n")
    log.write("Getting thresholds...\n")
    thresholds = getThresholds()

    # Set arrays to store data for each subfolder

    allChannels = []
    allOutputNames = []
    for subFolder in directories:
        chan, outputName = getChannels(subFolder)
        allChannels.append(chan)
        allOutputNames.append(outputName)

    # Loop that goes through each sub folder.

    log.write("_______________________________________________________________________\n")
    log.write("Beginning main directory loop: \n")
    log.write("\n")
    for inde, subFolder in enumerate(directories):

        log.write("______________________________________\n")
        log.write("Subfolder: " + subFolder + "\n")
        log.write("\n")

        channels = allChannels[inde]
        outputName = allOutputNames[inde]

        log.write("Channels: " + str(channels) + "\n")
        log.write("Output Name: " + outputName + "\n")

        open(outputDirectory + "/" + outputName + ".csv", 'w').close
        lowerBounds = [40, 20, 35, 50, 40]
        for chan in channels:
            v, x = chan
            if v in thresholds:
                lowerBounds[x] = int(thresholds[v])

        log.write("Lower Bound Thresholds: " + str(lowerBounds) + "\n")

        # Finds all correct EVOS split channel ch0 files and runs through them one at a time (see main loop process() on top)

        log.write("_________________________\n")
        log.write("Begining loop for each image \n")

        for filename in os.listdir(inputDirectory + subFolder):
            if filename.endswith("ch00.tif"):
                log.write("Processing: " + filename + " \n")
                process(subFolder, outputDirectory, filename);
        log.write("_________________________\n")
        log.write("Completed subfolder " + subFolder + ".  \n")
        log.write("\n")

    cat = """

      \    /\           Macro completed!
       )  ( ')   meow!
      (  /  )
       \(__)|"""

    log.write(cat)

print(cat)
