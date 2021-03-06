"""
Module that contains the ARCSIHyperion class.
"""
############################################################################
#  arcsisensorhyperion.py
############################################################################

# Import updated print function into python 2.7
from __future__ import print_function
# Import updated division operator into python 2.7
from __future__ import division
# import abstract base class stuff
from .arcsisensor import ARCSIAbstractSensor
# Import the ARCSI exception class
from .arcsiexception import ARCSIException
# Import the ARCSI utilities class
from .arcsiutils import ARCSIUtils
# Import the datetime module
import datetime
# Import the GDAL/OGR spatial reference library
from osgeo import osr
from osgeo import ogr
# Import OS path module for manipulating the file system
import os.path
# Import the RSGISLib Module.
import rsgislib
# Import the RSGISLib Image Calibration Module.
import rsgislib.imagecalibration
# Import the RSGISLib Image Utilities Module.
import rsgislib.imageutils
# Import the RSGISLib Image Calculations Module.
import rsgislib.imagecalc
# Import the RSGISLib segmentation Module
import rsgislib.segmentation
# Import the RSGISLib segmentation Module
import rsgislib.segmentation.segutils
# Import the RSGISLib Raster GIS Module
import rsgislib.rastergis
# Import the collections module
import collections
# Import the py6s module for running 6S from python.
import Py6S
# Import the python maths library
import math
# Import the RIOS RAT library
from rios import rat
# Import the GDAL python library
import osgeo.gdal as gdal
# Import the scipy optimisation library - used for finding AOD values form the imagery.
from scipy.optimize import minimize
# Import the numpy library
import numpy
# Import JSON module
import json
# Import the shutil module
import shutil
# Import the solar angle tools from RSGISLib
import rsgislib.imagecalibration.solarangles
# Using python-fmask (http://pythonfmask.org)
import fmask.landsatangles
import fmask.config
import fmask.fmask
import rios.fileinfo


class ARCSIHyperionSensor (ARCSIAbstractSensor):
    """
    A class which represents the landsat 8 sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        #Change this to array

        self.bandFile = [""]*242
        self.wave = [0.0]*242
        self.irr = [0.0]*242
        self.sensor = "hyp"

        self.vnir_rad_scale = 0.0
        self.swir_rad_scale = 0.0
        self.FWHM = [0.0]*242
        self.sensorID = ""
        self.spacecraftID = ""
        self.cloudCover = 0.0
        self.cloudCoverLand = 0.0
        self.earthSunDistance = 0.0
        self.gridCellSizePan = 0.0
        self.gridCellSizeRefl = 0.0
        self.gridCellSizeTherm = 0.0

    def extractHeaderParameters(self, inputHeader, wktStr):
        try:
            if not self.userSpInputImage is None:
                raise ARCSIException("Hyperion sensor cannot accept a user specified image file - only the images in the header file will be used.")
            self.headerFileName = os.path.split(inputHeader)[1]

            arcsiUtils = ARCSIUtils()
            print(inputHeader)
            print("Reading header file")
            pPath = os.path.dirname(inputHeader)+"/properties.txt"
            pFile = open(pPath, 'r')
            hFile = open(inputHeader, 'r')
            prop = dict()
            headerParams = dict()
            #opens the properties file for wavelength and irradiance data.
            for line in pFile:
                line = line.strip()
                if line:
                    lineVals = line.split('=')
                    if len(lineVals) == 2:
                        if (lineVals[0].strip() != "GROUP") or (lineVals[0].strip() != "END_GROUP"):
                            prop[lineVals[0].strip()] = lineVals[1].strip().replace('"','')
            for line in hFile:
                line = line.strip()
                if line:
                    lineVals = line.split('=')
                    if len(lineVals) == 2:
                        if (lineVals[0].strip() != "GROUP") or (lineVals[0].strip() != "END_GROUP"):
                            headerParams[lineVals[0].strip()] = lineVals[1].strip().replace('"','')
            hFile.close()
            print("Extracting Header Values")
            # Get the sensor info.
            if ((headerParams["SPACECRAFT_ID"].upper() == "EO1") or (headerParams["SPACECRAFT_ID"].upper() == "EO-1")) and (headerParams["SENSOR_ID"].upper() == "HYPERION"):
                self.sensor = "hyp"
            else:
                raise ARCSIException("Do no recognise the spacecraft and sensor or combination.")

            self.sensorID = headerParams["SENSOR_ID"]
            self.spacecraftID = headerParams["SPACECRAFT_ID"]

            # Get date and time of the acquisition
            acData = headerParams["ACQUISITION_DATE"].split('-')
            acTime = headerParams["START_TIME"].split(' ')
            secsTime = acTime[2].split(':')
            self.acquisitionTime = datetime.datetime(year=int(acData[0]),month=int(acData[1]),day= int(acData[2]),hour=int(secsTime[0]),minute=int(secsTime[1]),second=int(secsTime[2]))
            #
            #TODO: Go Ahead and make sure the metadata here is correct from the file
            self.solarZenith = arcsiUtils.str2Float(headerParams["SENSOR_LOOK_ANGLE"])
            self.solarAzimuth = 98 

            # Get the geographic lat/long corners of the image.
            self.latTL = arcsiUtils.str2Float(headerParams["PRODUCT_UL_CORNER_LAT"])
            self.lonTL = arcsiUtils.str2Float(headerParams["PRODUCT_UL_CORNER_LON"])
            self.latTR = arcsiUtils.str2Float(headerParams["PRODUCT_UR_CORNER_LAT"])
            self.lonTR = arcsiUtils.str2Float(headerParams["PRODUCT_UR_CORNER_LON"])
            self.latBL = arcsiUtils.str2Float(headerParams["PRODUCT_LL_CORNER_LAT"])
            self.lonBL = arcsiUtils.str2Float(headerParams["PRODUCT_LL_CORNER_LON"])
            self.latBR = arcsiUtils.str2Float(headerParams["PRODUCT_LR_CORNER_LAT"])
            self.lonBR = arcsiUtils.str2Float(headerParams["PRODUCT_LR_CORNER_LON"])

            # Get the projected X/Y corners of the image
            self.xTL = arcsiUtils.str2Float(headerParams["PRODUCT_UL_CORNER_MAPX"])
            self.yTL = arcsiUtils.str2Float(headerParams["PRODUCT_UL_CORNER_MAPY"])
            self.xTR = arcsiUtils.str2Float(headerParams["PRODUCT_UR_CORNER_MAPX"])
            self.yTR = arcsiUtils.str2Float(headerParams["PRODUCT_UR_CORNER_MAPY"])
            self.xBL = arcsiUtils.str2Float(headerParams["PRODUCT_LL_CORNER_MAPX"])
            self.yBL = arcsiUtils.str2Float(headerParams["PRODUCT_LL_CORNER_MAPY"])
            self.xBR = arcsiUtils.str2Float(headerParams["PRODUCT_LR_CORNER_MAPX"])
            self.yBR = arcsiUtils.str2Float(headerParams["PRODUCT_LR_CORNER_MAPY"])

            # Get projection
            inProj = osr.SpatialReference()
            if (headerParams["MAP_PROJECTION"] == "UTM") and (headerParams["REFERENCE_DATUM"] == "WGS84") and (headerParams["REFERENCE_ELLIPSOID"] == "WGS84"):
                utmZone = int(headerParams["ZONE_NUMBER"])
                utmCode = "WGS84UTM" + str(utmZone) + str("N")
                #print("UTM: ", utmCode)
                inProj.ImportFromEPSG(self.epsgCodes[utmCode])
            elif (headerParams["MAP_PROJECTION"] == "PS") and (headerParams["REFERENCE_DATUM"] == "WGS84") and (headerParams["REFERENCE_ELLIPSOID"] == "WGS84"):
                inProj.ImportFromWkt("PROJCS[\"PS WGS84\", GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563, AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Polar_Stereographic\"],PARAMETER[\"latitude_of_origin\",-71],PARAMETER[\"central_meridian\",0],PARAMETER[\"scale_factor\",1],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]]]")
            else:
                raise ARCSIException("Expecting Hyperion to be projected in UTM or PolarStereographic (PS) with datum=WGS84 and ellipsoid=WGS84.")

            if self.inWKT is "":
                self.inWKT = inProj.ExportToWkt()

            # Check image is square!
            if not ((self.xTL == self.xBL) and (self.yTL == self.yTR) and (self.xTR == self.xBR) and (self.yBL == self.yBR)):
                raise ARCSIException("Image is not square in projected coordinates.")

            self.xCentre = self.xTL + ((self.xTR - self.xTL)/2)
            self.yCentre = self.yBR + ((self.yTL - self.yBR)/2)

            self.lonCentre, self.latCentre = arcsiUtils.getLongLat(inProj, self.xCentre, self.yCentre)

            #print("Lat: " + str(self.latCentre) + " Long: " + str(self.lonCentre))
            filesDIR = os.path.dirname(inputHeader)
            #opens the Properties.txt and extracts wavelength and other parameters from it.
            for i in range(242):
                self.bandFile[i] = os.path.join(filesDIR, headerParams["BAND%s_FILE_NAME" %str(i+1)])
                self.wave[i] = float(prop["Wavelengths %s"%str(i+1)])
                self.wave[i] = self.wave[i]/1000
                self.irr[i] = float(prop["Irradiance %s"%str(i+1)])
                self.FWHM[i] = float(prop["FWHM %s"%str(i+1)])/1000
            #Gets the scaling factors from the header file
            self.vnir_rad_scale = 1/arcsiUtils.str2Float(headerParams["SCALING_FACTOR_VNIR"])
            self.swir_rad_scale = 1/arcsiUtils.str2Float(headerParams["SCALING_FACTOR_SWIR"])

            if "CLOUD_COVER" in headerParams:
                self.cloudCover = arcsiUtils.str2Float(headerParams["CLOUD_COVER"], 0.0)
            if "CLOUD_COVER_LAND" in headerParams:
                self.cloudCoverLand = arcsiUtils.str2Float(headerParams["CLOUD_COVER_LAND"], 0.0)
            if "EARTH_SUN_DISTANCE" in headerParams:
                self.earthSunDistance = arcsiUtils.str2Float(headerParams["EARTH_SUN_DISTANCE"], 0.0)
            if "GRID_CELL_SIZE" in headerParams:
                self.gridCellSizeRefl = arcsiUtils.str2Float(headerParams["GRID_CELL_SIZE"], 30.0)


            # Read MTL header into python dict for python-fmask
            self.fmaskMTLInfo = fmask.config.readMTLFile(inputHeader)

            fileDateStr = headerParams["PRODUCT_CREATION_TIME"].strip()
            fileDateStr = fileDateStr.replace('Z', '')
            self.fileDateObj = datetime.datetime.strptime(fileDateStr, "%Y-%m-%dT%H:%M:%S")

        except Exception as e:
            raise e

    def getSolarIrrStdSolarGeom(self):
        """
        Get Solar Azimuth and Zenith as standard geometry.
        Azimuth: N=0, E=90, S=180, W=270.
        """
        solarAz = rsgislib.imagecalibration.solarangles.getSolarIrrConventionSolarAzimuthFromUSGS(self.solarAzimuth)
        return (solarAz, self.solarZenith)

    def getSensorViewGeom(self):
        """
        Get sensor viewing angles
        returns (viewAzimuth, viewZenith)
        """
        return (0.0, 0.0)

    def generateOutputBaseName(self):
        """
        Provides an implementation for the Hyperion sensor
        """
        #rowpath = "r" + str(self.row) + "p" + str(self.path)
        name = self.bandFile[0].split("\\")
        name = name[len(name)-1].split("_")[0]
        outname = name
        return outname

    def generateMetaDataFile(self, outputPath, outputFileName, productsStr, validMaskImage="", footprintCalc=False, calcdValuesDict=dict(), outFilesDict=dict()):
        """
        Generate file metadata.
        """
        outJSONFilePath = os.path.join(outputPath, outputFileName)
        jsonData = self.getJSONDictDefaultMetaData(productsStr, validMaskImage, footprintCalc, calcdValuesDict, outFilesDict)
        sensorInfo = jsonData['SensorInfo']
        # sensorInfo['Row'] = self.row
        # sensorInfo['Path'] = self.path
        sensorInfo['SensorID'] = self.sensorID
        sensorInfo['SpacecraftID'] = self.spacecraftID
        acqDict = jsonData['AcquasitionInfo']
        imgInfo = dict()
        imgInfo['CellSizePan'] = self.gridCellSizePan
        jsonData['ImageInfo'] = imgInfo

        with open(outJSONFilePath, 'w') as outfile:
            json.dump(jsonData, outfile, sort_keys=True,indent=4, separators=(',', ': '), ensure_ascii=False)

    def expectedImageDataPresent(self):
        imageDataPresent = True
        for i in range(242):
            if not os.path.exists(self.bandFile[0]):
                imageDataPresent = False
        #if not os.path.exists(self.bandQAFile):
        #    imageDataPresent = False

        return imageDataPresent

    def hasThermal(self):
        return False

    def applyImageDataMask(self, inputHeader, outputPath, outputMaskName, outputImgName, outFormat, outWKTFile):
        raise ARCSIException("Hyperion does not provide any image masks, do not use the MASK option.")

    def mosaicImageTiles(self, outputPath):
        raise ARCSIException("Image data does not need mosaicking")

    def resampleImgRes(self, outputPath, resampleToLowResImg, resampleMethod='cubic', multicore=False):
        raise ARCSIException("Image data does not need resampling")

    def sharpenLowResRadImgBands(self, inputImg, outputImage, outFormat):
        raise ARCSIException("Image sharpening is not available for this sensor.")

    # def generateValidImageDataMask(self, outputPath, outputMaskName, viewAngleImg, outFormat):
    #     def sunAnglesForExtent(imgInfo, mtlInfo):
    #         cornerLatLong = imgInfo.getCorners(outEPSG=4326)
    #         (ul_long, ul_lat, ur_long, ur_lat, lr_long, lr_lat, ll_long, ll_lat) = cornerLatLong
    #         pts = numpy.array([
    #             [ul_long, ul_lat],
    #             [ur_long, ur_lat],
    #             [ll_long, ll_lat],
    #             [lr_long, lr_lat]
    #         ])
    #         longDeg = pts[:, 0]
    #         latDeg = pts[:, 1]

    #         # Date/time in UTC
    #         dateStr = mtlInfo['DATE_ACQUIRED']
    #         timeStr = (mtlInfo['START_TIME'].replace('Z', '')).split(" ")[2]
    #         print(timeStr)
    #         ymd = [int(i) for i in dateStr.split('-')]
    #         dateObj = datetime.date(ymd[0], ymd[1], ymd[2])
    #         julianDay = (dateObj - datetime.date(ymd[0], 1, 1)).days + 1
    #         juldayYearEnd = (datetime.date(ymd[0], 12, 31) - datetime.date(ymd[0], 1, 1)).days + 1
    #         # Julian day as a proportion of the year
    #         jdp = julianDay / juldayYearEnd
    #         # Hour in UTC
    #         hms = [float(x) for x in timeStr.split(':')]
    #         hourGMT = hms[0] + hms[1] / 60.0 + hms[2] / 3600.0

    #         (sunAz, sunZen) = fmask.landsatangles.sunAnglesForPoints(latDeg, longDeg, hourGMT, jdp)


    #         sunAngles = numpy.vstack((sunAz, sunZen)).T
    #         return sunAngles
    #     print("Create the valid data mask")
    #     tmpBaseName = os.path.splitext(outputMaskName)[0]
    #     tmpValidPxlMsk = os.path.join(outputPath, tmpBaseName+'vldpxlmsk.kea')
    #     outputImage = os.path.join(outputPath, outputMaskName)
    #     inImages = self.bandFile
    #     rsgislib.imageutils.genValidMask(inimages=inImages, outimage=tmpValidPxlMsk, gdalformat='KEA', nodata=0.0)
    #     rsgislib.rastergis.populateStats(tmpValidPxlMsk, True, False, True)
    #     # Check there is valid data
    #     ratDS = gdal.Open(tmpValidPxlMsk, gdal.GA_ReadOnly)
    #     Histogram = rat.readColumn(ratDS, "Histogram")
    #     ratDS = None
    #     if Histogram.shape[0] < 2:
    #         raise ARCSIException("There is no valid data in this image.")
    #     if not os.path.exists(viewAngleImg):
    #         print("Calculate Image Angles.")
    #         imgInfo = rios.fileinfo.ImageInfo(tmpValidPxlMsk)
    #         corners = fmask.landsatangles.findImgCorners(tmpValidPxlMsk, imgInfo)
    #         nadirLine = fmask.landsatangles.findNadirLine(corners)
    #         extentSunAngles = sunAnglesForExtent(imgInfo, self.fmaskMTLInfo)
    #         satAzimuth = fmask.landsatangles.satAzLeftRight(nadirLine)
    #         fmask.landsatangles.makeAnglesImage(tmpValidPxlMsk, viewAngleImg, nadirLine, extentSunAngles, satAzimuth, imgInfo)
    #         dataset = gdal.Open(viewAngleImg, gdal.GA_Update)
    #         if not dataset is None:
    #             dataset.GetRasterBand(1).SetDescription("SatelliteAzimuth")
    #             dataset.GetRasterBand(2).SetDescription("SatelliteZenith")
    #             dataset.GetRasterBand(3).SetDescription("SolorAzimuth")
    #             dataset.GetRasterBand(4).SetDescription("SolorZenith")
    #         dataset = None
    #     rsgislib.imagecalc.bandMath(outputImage, '(VA<14)&&(VM==1)?1:0', outFormat, rsgislib.TYPE_8UINT, [rsgislib.imagecalc.BandDefn('VA', viewAngleImg, 2), rsgislib.imagecalc.BandDefn('VM', tmpValidPxlMsk, 1)])
    #     rsgisUtils = rsgislib.RSGISPyUtils()
    #     rsgisUtils.deleteFileWithBasename(tmpValidPxlMsk)
    #     return outputImage

    def convertImageToRadiance(self, outputPath, outputReflName, outputThermalName, outFormat):
        print("Converting to Radiance")
        outputReflImage = os.path.join(outputPath, outputReflName)
        outputThermalImage = None
        bandDefnSeq = list()
        lsBand = collections.namedtuple('LSBand', ['bandName', 'fileName', 'bandIndex', 'addVal', 'multiVal'])
        #Change this to for loop? TODO:ForLoop
        for i in range(242):
            if(i < 70):
                bandDefnSeq.append(lsBand(bandName="Band%s"%str(i+1), fileName=self.bandFile[i], bandIndex=1, addVal=0.0, multiVal=self.vnir_rad_scale))
            else:
                bandDefnSeq.append(lsBand(bandName="Band%s"%str(i+1), fileName=self.bandFile[i], bandIndex=1, addVal=0.0, multiVal=self.swir_rad_scale))
        print("ok")

        rsgislib.imagecalibration.landsat2RadianceMultiAdd(outputReflImage, outFormat, bandDefnSeq)

        return outputReflImage, outputThermalImage
        
    #Try this out and see how it changes it
    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        lsBand = collections.namedtuple('LSBand', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
        for i in range(242):
            bandDefnSeq.append(lsBand(bandName="Band%s"%str(i+1), fileName=self.bandFile[i], bandIndex=1, satVal=self.b1CalMax))

        rsgislib.imagecalibration.saturatedPixelsMask(outputImage, outFormat, bandDefnSeq)

        return outputImage

    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        for i in range(242):
            solarIrradianceVals.append(IrrVal(irradiance=self.irr[i]))

        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage

    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, inputViewAngleImg, inputValidImg, outputPath, outputName, outFormat, tmpPath, scaleFactor, cloud_msk_methods=None):
        try:
            arcsiUtils = ARCSIUtils()
            rsgisUtils = rsgislib.RSGISPyUtils()
            outputImage = os.path.join(outputPath, outputName)
            tmpBaseName = os.path.splitext(outputName)[0]
            imgExtension = arcsiUtils.getFileExtension(outFormat)
            tmpBaseDIR = os.path.join(tmpPath, tmpBaseName)

            tmpDIRExisted = True
            if not os.path.exists(tmpBaseDIR):
                os.makedirs(tmpBaseDIR)
                tmpDIRExisted = False

            if (cloud_msk_methods is None) or (cloud_msk_methods == 'FMASK'):
                tmpFMaskOut = os.path.join(tmpBaseDIR, tmpBaseName+'_pyfmaskout.kea')

                # Create tmp TOA stack with Band 9 (Cirrus)
                tmpTOAB9Out = os.path.join(tmpBaseDIR, tmpBaseName+'_TOAB9.kea')
                toaExp = '((b1*'+str(self.b9ReflMulti)+')+'+str(self.b9ReflAdd)+')*'+str(scaleFactor)
                rsgislib.imagecalc.imageMath(self.band9File, tmpTOAB9Out, toaExp, 'KEA', rsgislib.TYPE_16UINT, False)
                if not rsgisUtils.doGDALLayersHaveSameProj(inputReflImage, tmpTOAB9Out):
                    tmpTOAB9OutNotProj = tmpTOAB9Out
                    tmpTOAB9Out = os.path.join(tmpBaseDIR, tmpBaseName+'_TOAB9_reproj.kea')
                    rsgislib.imageutils.resampleImage2Match(inputReflImage, tmpTOAB9OutNotProj, tmpTOAB9Out, 'KEA', 'cubic', rsgislib.TYPE_16UINT)

                tmpReflStackOut = os.path.join(tmpBaseDIR, tmpBaseName+'_TOAreflStack.kea')
                rsgislib.imageutils.stackImageBands([inputReflImage, tmpTOAB9Out], None, tmpReflStackOut, None, 0, 'KEA', rsgislib.TYPE_16UINT)

                tmpThermStack = os.path.join(tmpBaseDIR, tmpBaseName + '_ThermB10B11Stack.kea')
                rsgislib.imageutils.stackImageBands([self.band10File, self.band11File], None, tmpThermStack, None, 0, 'KEA', rsgislib.TYPE_16UINT)

                tmpThermalLayer = tmpThermStack
                if not rsgisUtils.doGDALLayersHaveSameProj(inputThermalImage, tmpThermStack):
                    tmpThermalLayer = os.path.join(tmpBaseDIR, tmpBaseName+'_thermalresample.kea')
                    rsgislib.imageutils.resampleImage2Match(inputThermalImage, tmpThermStack, tmpThermalLayer, 'KEA', 'cubic', rsgislib.TYPE_32FLOAT)

                minCloudSize = 0
                cloudBufferDistance = 150
                shadowBufferDistance = 300

                fmaskFilenames = fmask.config.FmaskFilenames()
                fmaskFilenames.setTOAReflectanceFile(tmpReflStackOut)
                fmaskFilenames.setThermalFile(tmpThermalLayer)
                fmaskFilenames.setSaturationMask(inputSatImage)
                fmaskFilenames.setOutputCloudMaskFile(tmpFMaskOut)

                thermalGain1040um = self.b10RadMulti
                thermalOffset1040um = self.b10RadAdd
                thermalBand1040um = 0
                thermalInfo = fmask.config.ThermalFileInfo(thermalBand1040um, thermalGain1040um, thermalOffset1040um, self.k1ConstB10, self.k2ConstB10)

                anglesInfo = fmask.config.AnglesFileInfo(inputViewAngleImg, 3, inputViewAngleImg, 2, inputViewAngleImg, 1, inputViewAngleImg, 0)

                fmaskConfig = fmask.config.FmaskConfig(fmask.config.FMASK_LANDSAT8)
                fmaskConfig.setTOARefScaling(float(scaleFactor))
                fmaskConfig.setThermalInfo(thermalInfo)
                fmaskConfig.setAnglesInfo(anglesInfo)
                fmaskConfig.setKeepIntermediates(False)
                fmaskConfig.setVerbose(True)
                fmaskConfig.setTempDir(tmpBaseDIR)
                fmaskConfig.setMinCloudSize(minCloudSize)
                fmaskConfig.setEqn17CloudProbThresh(fmask.config.FmaskConfig.Eqn17CloudProbThresh)
                fmaskConfig.setEqn20NirSnowThresh(fmask.config.FmaskConfig.Eqn20NirSnowThresh)
                fmaskConfig.setEqn20GreenSnowThresh(fmask.config.FmaskConfig.Eqn20GreenSnowThresh)

                # Work out a suitable buffer size, in pixels, dependent on the resolution of the input TOA image
                toaImgInfo = rios.fileinfo.ImageInfo(inputReflImage)
                fmaskConfig.setCloudBufferSize(int(cloudBufferDistance / toaImgInfo.xRes))
                fmaskConfig.setShadowBufferSize(int(shadowBufferDistance / toaImgInfo.xRes))

                fmask.fmask.doFmask(fmaskFilenames, fmaskConfig)

                rsgislib.imagecalc.imageMath(tmpFMaskOut, outputImage, '(b1==2)?1:(b1==3)?2:0', outFormat, rsgislib.TYPE_8UINT)

            elif (cloud_msk_methods == 'LSMSK'):
                if (self.bandQAFile == "") or (not os.path.exists(self.bandQAFile)):
                    raise ARCSIException("The QA band is not present - cannot use this for cloud masking.")

                bqa_img_file = self.bandQAFile
                if not rsgisUtils.doGDALLayersHaveSameProj(bqa_img_file, inputReflImage):
                    bqa_img_file = os.path.join(tmpBaseDIR, tmpBaseName+'_BQA.kea')
                    rsgislib.imageutils.resampleImage2Match(inputReflImage, self.bandQAFile, bqa_img_file, 'KEA',
                                                            'nearestneighbour', rsgislib.TYPE_16UINT, noDataVal=0,
                                                            multicore=False)

                exp = '(b1==2800)||(b1==2804)||(b1==2808)||(b1==2812)||(b1==6896)||(b1==6900)||(b1==6904)||(b1==6908)?1:' \
                      '(b1==2976)||(b1==2980)||(b1==2984)||(b1==2988)||(b1==3008)||(b1==3012)||(b1==3016)||(b1==3020)||' \
                      '(b1==7072)||(b1==7076)||(b1==7080)||(b1==7084)||(b1==7104)||(b1==7108)||(b1==7112)||(b1==7116)?2:0'
                rsgislib.imagecalc.imageMath(bqa_img_file, outputImage, exp, outFormat, rsgislib.TYPE_8UINT)

            else:
                raise ARCSIException("Landsat only has FMASK and LSMSK cloud masking options; option provided is unknown.")

            if outFormat == 'KEA':
                rsgislib.rastergis.populateStats(outputImage, True, True)
                ratDataset = gdal.Open(outputImage, gdal.GA_Update)
                red = rat.readColumn(ratDataset, 'Red')
                green = rat.readColumn(ratDataset, 'Green')
                blue = rat.readColumn(ratDataset, 'Blue')
                ClassName = numpy.empty_like(red, dtype=numpy.dtype('a255'))

                red[0] = 0
                green[0] = 0
                blue[0] = 0

                if (red.shape[0] == 2) or (red.shape[0] == 3):
                    red[1] = 0
                    green[1] = 0
                    blue[1] = 255
                    ClassName[1] = 'Clouds'

                    if (red.shape[0] == 3):
                        red[2] = 0
                        green[2] = 255
                        blue[2] = 255
                        ClassName[2] = 'Shadows'

                rat.writeColumn(ratDataset, "Red", red)
                rat.writeColumn(ratDataset, "Green", green)
                rat.writeColumn(ratDataset, "Blue", blue)
                rat.writeColumn(ratDataset, "ClassName", ClassName)
                ratDataset = None
            rsgislib.imageutils.copyProjFromImage(outputImage, inputReflImage)

            if not self.debugMode:
                if not tmpDIRExisted:
                    shutil.rmtree(tmpBaseDIR, ignore_errors=True)

            return outputImage
        except Exception as e:
            raise e

    def createCloudMaskDataArray(self, inImgDataArr):
        return inImgDataArr

    def defineDarkShadowImageBand(self):
        return 5

    def calc6SCoefficients(self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        sixsCoeffs = numpy.zeros((242, 6), dtype=numpy.float32)
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        #s.ground_reflectance = Py6S.GroundReflectance.HomogeneousHapke(0.101, -0.263, 0.589, 0.046)
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.Landsat_TM()
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute)/60.0
        s.geometry.latitude = self.latCentre
        s.geometry.longitude = self.lonCentre
        s.altitudes = Py6S.Altitudes()
        s.altitudes.set_target_custom_altitude(surfaceAltitude)
        s.altitudes.set_sensor_satellite_level()
        if useBRDF:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrBRDFFromRadiance(200)
        else:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal

        # Band 1
        for i in range(242):
            s.wavelength = Py6S.Wavelength(((self.wave[i]+self.wave[i]+.01)/2-self.FWHM[i]/2),((self.wave[i]+self.wave[i]+.01)/2+self.FWHM[i]/2))

            #s.wavelength = Py6S.Wavelength(self.wave[i], self.wave[i]+.01,  [0.000010, 0.000117, 0.000455, 0.001197, 0.006869])
            s.run()
            sixsCoeffs[i,0] = float(s.outputs.values['coef_xa'])
            sixsCoeffs[i,1] = float(s.outputs.values['coef_xb'])
            sixsCoeffs[i,2] = float(s.outputs.values['coef_xc'])
            sixsCoeffs[i,3] = float(s.outputs.values['direct_solar_irradiance'])
            sixsCoeffs[i,4] = float(s.outputs.values['diffuse_solar_irradiance'])
            sixsCoeffs[i,5] = float(s.outputs.values['environmental_irradiance'])

        return sixsCoeffs

    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF, scaleFactor):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        imgBandCoeffs = list()

        sixsCoeffs = self.calc6SCoefficients(aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF)
        for i in range(242):
            imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=(i+1), aX=float(sixsCoeffs[i,0]), bX=float(sixsCoeffs[i,1]), cX=float(sixsCoeffs[i,2]), DirIrr=float(sixsCoeffs[i,3]), DifIrr=float(sixsCoeffs[i,4]), EnvIrr=float(sixsCoeffs[i,5])))
        print("here")
        rsgislib.imagecalibration.apply6SCoeffSingleParam(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, imgBandCoeffs)
        
        return outputImage

    def convertImageToSurfaceReflDEMElevLUT(self, inputRadImage, inputDEMFile, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, scaleFactor, elevCoeffs=None):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevCoeffs is None:
            print("Build an LUT for elevation values.")
            elev6SCoeffsLUT = self.buildElevation6SCoeffLUT(aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax)
            print("LUT has been built.")

            elevCoeffs = list()
            for elevLUT in elev6SCoeffsLUT:
                imgBandCoeffs = list()
                sixsCoeffs = elevLUT.Coeffs
                elevVal = elevLUT.Elev
                for i in range(242):
                    imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=(i+1), aX=float(sixsCoeffs[i,0]), bX=float(sixsCoeffs[i,1]), cX=float(sixsCoeffs[i,2]), DirIrr=float(sixsCoeffs[i,3]), DifIrr=float(sixsCoeffs[i,4]), EnvIrr=float(sixsCoeffs[i,5])))

                elevCoeffs.append(rsgislib.imagecalibration.ElevLUTFeat(Elev=float(elevVal), Coeffs=imgBandCoeffs))

        rsgislib.imagecalibration.apply6SCoeffElevLUTParam(inputRadImage, inputDEMFile, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, elevCoeffs)
        return outputImage, elevCoeffs


    def convertImageToSurfaceReflAOTDEMElevLUT(self, inputRadImage, inputDEMFile, inputAOTImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax, scaleFactor, elevAOTCoeffs=None):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevAOTCoeffs is None:
            print("Build an LUT for elevation and AOT values.")
            elevAOT6SCoeffsLUT = self.buildElevationAOT6SCoeffLUT(aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax)

            elevAOTCoeffs = list()
            for elevLUT in elevAOT6SCoeffsLUT:
                elevVal = elevLUT.Elev
                aotLUT = elevLUT.Coeffs
                aot6SCoeffsOut = list()
                for aotFeat in aotLUT:
                    sixsCoeffs = aotFeat.Coeffs
                    aotVal = aotFeat.AOT
                    imgBandCoeffs = list()
                    for i in range(242):
                        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=(i+1), aX=float(sixsCoeffs[i,0]), bX=float(sixsCoeffs[i,1]), cX=float(sixsCoeffs[i,2]), DirIrr=float(sixsCoeffs[i,3]), DifIrr=float(sixsCoeffs[i,4]), EnvIrr=float(sixsCoeffs[i,5])))

                    aot6SCoeffsOut.append(rsgislib.imagecalibration.AOTLUTFeat(AOT=float(aotVal), Coeffs=imgBandCoeffs))
                elevAOTCoeffs.append(rsgislib.imagecalibration.ElevLUTFeat(Elev=float(elevVal), Coeffs=aot6SCoeffsOut))

        rsgislib.imagecalibration.apply6SCoeffElevAOTLUTParam(inputRadImage, inputDEMFile, inputAOTImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, elevAOTCoeffs)

        return outputImage, elevAOTCoeffs

    def run6SToOptimiseAODValue(self, aotVal, radBlueVal, predBlueVal, aeroProfile, atmosProfile, grdRefl, surfaceAltitude):
        """Used as part of the optimastion for identifying values of AOD"""
        print("Testing AOD Val: ", aotVal,)
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.Landsat_TM()
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute)/60.0
        s.geometry.latitude = self.latCentre
        s.geometry.longitude = self.lonCentre
        s.altitudes = Py6S.Altitudes()
        s.altitudes.set_target_custom_altitude(surfaceAltitude)
        s.altitudes.set_sensor_satellite_level()
        s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal

        # Band 2 (Blue!) -- For Hyperion Band 16
        s.wavelength = Py6S.Wavelength(self.wave[15])
        s.run()
        aX = float(s.outputs.values['coef_xa'])
        bX = float(s.outputs.values['coef_xb'])
        cX = float(s.outputs.values['coef_xc'])

        tmpVal = (aX*radBlueVal)-bX
        reflBlueVal = tmpVal/(1.0+cX*tmpVal)
        outDist = math.sqrt(math.pow((reflBlueVal - predBlueVal),2))
        print("\taX: ", aX, " bX: ", bX, " cX: ", cX, "     Dist = ", outDist)
        return outDist

    def findDDVTargets(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath):
        try:
            print("Finding dark targets.")
            arcsiUtils = ARCSIUtils()
            tmpBaseName = os.path.splitext(outputName)[0]
            thresImage = os.path.join(tmpPath, tmpBaseName+"_thresd"+arcsiUtils.getFileExtension(outFormat))
            thresImageClumps = os.path.join(tmpPath, tmpBaseName+"_thresdclumps"+arcsiUtils.getFileExtension(outFormat))
            thresImageClumpsRMSmall = os.path.join(tmpPath, tmpBaseName+"_thresdclumpsgt10"+arcsiUtils.getFileExtension(outFormat))
            thresImageClumpsFinal = os.path.join(tmpPath, tmpBaseName+"_thresdclumpsFinal"+arcsiUtils.getFileExtension(outFormat))

            percentiles = rsgislib.imagecalc.bandPercentile(inputTOAImage, 0.05, 0)
            if percentiles[6] > 30:
                b7Thres = str(percentiles[6])
            else:
                b7Thres = "30.0"
            print("SWIR DDV Threshold = ", b7Thres)

            thresMathBands = list()
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b4', fileName=inputTOAImage, bandIndex=4))
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b5', fileName=inputTOAImage, bandIndex=5))
            thresMathBands.append(rsgislib.imagecalc.BandDefn(bandName='b7', fileName=inputTOAImage, bandIndex=7))
            rsgislib.imagecalc.bandMath(thresImage, "(b7<" + b7Thres + ")&&(b7!=0)&&(((b5-b4)/(b5+b4))>0.1)?1:0", outFormat, rsgislib.TYPE_8UINT, thresMathBands)
            rsgislib.segmentation.clump(thresImage, thresImageClumps, outFormat, False, 0.0)
            rsgislib.rastergis.populateStats(thresImageClumps, True, True)
            rsgislib.segmentation.rmSmallClumps(thresImageClumps, thresImageClumpsRMSmall, 100, outFormat)
            rsgislib.segmentation.relabelClumps(thresImageClumpsRMSmall, thresImageClumpsFinal, outFormat, False)
            rsgislib.rastergis.populateStats(thresImageClumpsFinal, True, True)

            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName(outFormat)
                gdalDriver.Delete(thresImage)
                gdalDriver.Delete(thresImageClumps)
                gdalDriver.Delete(thresImageClumpsRMSmall)

            return thresImageClumpsFinal
        except Exception as e:
            raise e

    def estimateImageToAODUsingDDV(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax):
        print("Estimating AOD through Blue - SWIR relationship.")
        try:
            outputAOTImage = os.path.join(outputPath, outputName)

            thresImageClumpsFinal = self.findDDVTargets(inputTOAImage, outputPath, outputName, "KEA", tmpPath)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanElev"))
            rsgislib.rastergis.populateRATWithStats(inputDEMFile, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=2, minField="MinB2TOA", meanField="MeanB2TOA"))
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=7, minField="MinB7TOA", meanField="MeanB7TOA"))
            rsgislib.rastergis.populateRATWithStats(inputTOAImage, thresImageClumpsFinal, stats2CalcTOA)
            stats2CalcRad = list()
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=2, minField="MinB2RAD", meanField="MeanB2RAD"))
            rsgislib.rastergis.populateRATWithStats(inputRADImage, thresImageClumpsFinal, stats2CalcRad)

            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            Histogram = rat.readColumn(ratDS, "Histogram")
            MeanElev = rat.readColumn(ratDS, "MeanElev")

            selected = Histogram * 2
            selected[...] = 1
            selected[0] = 0
            rat.writeColumn(ratDS, "Selected", selected)
            ratDS = None

            rsgislib.rastergis.spatialLocation(thresImageClumpsFinal, "Eastings", "Northings")
            rsgislib.rastergis.selectClumpsOnGrid(thresImageClumpsFinal, "Selected", "PredictAOTFor", "Eastings", "Northings", "MinB7TOA", "min", 10, 10)

            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            MeanB2TOA = rat.readColumn(ratDS, "MeanB2TOA")
            MeanB7TOA = rat.readColumn(ratDS, "MeanB7TOA")
            MeanB2RAD = rat.readColumn(ratDS, "MeanB2RAD")
            PredictAOTFor = rat.readColumn(ratDS, "PredictAOTFor")

            PredB2Refl = (MeanB7TOA/1000) * 0.33

            rat.writeColumn(ratDS, "PredB2Refl", PredB2Refl)

            numAOTValTests = int(math.ceil((aotValMax - aotValMin)/0.05))+1

            if not numAOTValTests >= 1:
                raise ARCSIException("min and max AOT range are too close together, they need to be at least 0.05 apart.")

            cAOT = aotValMin
            cDist = 0.0
            minAOT = 0.0
            minDist = 0.0

            aotVals = numpy.zeros_like(MeanB7TOA, dtype=numpy.float)

            for i in range(len(PredB2Refl)):
                if PredictAOTFor[i] == 1:
                    print("Predicting AOD for Segment ", i)
                    for j in range(numAOTValTests):
                        cAOT = aotValMin + (0.05 * j)
                        cDist = self.run6SToOptimiseAODValue(cAOT, MeanB2RAD[i], PredB2Refl[i], aeroProfile, atmosProfile, grdRefl, MeanElev[i]/1000)
                        if j == 0:
                            minAOT = cAOT
                            minDist = cDist
                        elif cDist < minDist:
                            minAOT = cAOT
                            minDist = cDist
                    #predAOTArgs = (MinB2RAD[i], PredB2Refl[i], aeroProfile, atmosProfile, grdRefl, MeanElev[i]/1000)
                    #res = minimize(self.run6SToOptimiseAODValue, minAOT, method='nelder-mead', options={'maxiter': 20, 'xtol': 0.001, 'disp': True}, args=predAOTArgs)
                    #aotVals[i] = res.x[0]
                    aotVals[i] = minAOT
                    print("IDENTIFIED AOT: ", aotVals[i])
                else:
                    aotVals[i] = 0
            rat.writeColumn(ratDS, "AOT", aotVals)

            Eastings = rat.readColumn(ratDS, "Eastings")
            Northings = rat.readColumn(ratDS, "Northings")
            ratDS = None

            Eastings = Eastings[PredictAOTFor!=0]
            Northings = Northings[PredictAOTFor!=0]
            aotVals = aotVals[PredictAOTFor!=0]

            interpSmoothing = 10.0
            self.interpolateImageFromPointData(inputTOAImage, Eastings, Northings, aotVals, outputAOTImage, outFormat, interpSmoothing, True, 0.05)

            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName("KEA")
                gdalDriver.Delete(thresImageClumpsFinal)

            return outputAOTImage
        except Exception as e:
            raise e

    def estimateImageToAODUsingDOS(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax, globalDOS, simpleDOS, dosOutRefl):
        try:
            print("Estimating AOD Using DOS")
            arcsiUtils = ARCSIUtils()
            outputAOTImage = os.path.join(outputPath, outputName)
            tmpBaseName = os.path.splitext(outputName)[0]
            imgExtension = arcsiUtils.getFileExtension(outFormat)

            dosBlueImage = ""
            minObjSize = 5
            darkPxlPercentile = 0.01
            blockSize = 1000
            if simpleDOS:
                outputDOSBlueName = tmpBaseName + "DOSBlue" + imgExtension
                dosBlueImage, bandOff = self.convertImageBandToReflectanceSimpleDarkSubtract(inputTOAImage, outputPath, outputDOSBlueName, outFormat, dosOutRefl, 2)
            elif globalDOS:
                dosBlueImage = self.performDOSOnSingleBand(inputTOAImage, 2, outputPath, tmpBaseName, "Blue", "KEA", tmpPath, minObjSize, darkPxlPercentile, dosOutRefl)
            else:
                dosBlueImage = self.performLocalDOSOnSingleBand(inputTOAImage, 2, outputPath, tmpBaseName, "Blue", "KEA", tmpPath, minObjSize, darkPxlPercentile, blockSize, dosOutRefl)

            thresImageClumpsFinal = os.path.join(tmpPath, tmpBaseName + "_clumps" + imgExtension)
            rsgislib.segmentation.segutils.runShepherdSegmentation(inputTOAImage, thresImageClumpsFinal, tmpath=tmpPath, gdalformat="KEA", numClusters=20, minPxls=10, bands=[5,6,4], processInMem=True)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanElev"))
            rsgislib.rastergis.populateRATWithStats(inputDEMFile, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanB2DOS"))
            rsgislib.rastergis.populateRATWithStats(dosBlueImage, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcRad = list()
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=2, meanField="MeanB2RAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=5, meanField="MeanB5RAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=4, meanField="MeanB4RAD"))
            rsgislib.rastergis.populateRATWithStats(inputRADImage, thresImageClumpsFinal, stats2CalcRad)

            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            Histogram = rat.readColumn(ratDS, "Histogram")
            MeanElev = rat.readColumn(ratDS, "MeanElev")

            MeanB5RAD = rat.readColumn(ratDS, "MeanB5RAD")
            MeanB4RAD = rat.readColumn(ratDS, "MeanB4RAD")

            radNDVI = (MeanB5RAD - MeanB4RAD)/(MeanB5RAD + MeanB4RAD)

            selected = Histogram * 2
            selected[...] = 0
            selected[radNDVI>0.2] = 1
            rat.writeColumn(ratDS, "Selected", selected)
            ratDS = None

            rsgislib.rastergis.spatialLocation(thresImageClumpsFinal, "Eastings", "Northings")
            rsgislib.rastergis.selectClumpsOnGrid(thresImageClumpsFinal, "Selected", "PredictAOTFor", "Eastings", "Northings", "MeanB2DOS", "min", 10, 10)

            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            MeanB2DOS = rat.readColumn(ratDS, "MeanB2DOS")
            MeanB2DOS = MeanB2DOS / 1000
            MeanB2RAD = rat.readColumn(ratDS, "MeanB2RAD")
            PredictAOTFor = rat.readColumn(ratDS, "PredictAOTFor")

            numAOTValTests = int(math.ceil((aotValMax - aotValMin)/0.05))+1

            if not numAOTValTests >= 1:
                raise ARCSIException("min and max AOT range are too close together, they need to be at least 0.05 apart.")

            cAOT = aotValMin
            cDist = 0.0
            minAOT = 0.0
            minDist = 0.0

            aotVals = numpy.zeros_like(MeanB2RAD, dtype=numpy.float)

            for i in range(len(MeanB2RAD)):
                if PredictAOTFor[i] == 1:
                    print("Predicting AOD for Segment ", i)
                    for j in range(numAOTValTests):
                        cAOT = aotValMin + (0.05 * j)
                        cDist = self.run6SToOptimiseAODValue(cAOT, MeanB2RAD[i], MeanB2DOS[i], aeroProfile, atmosProfile, grdRefl, MeanElev[i]/1000)
                        if j == 0:
                            minAOT = cAOT
                            minDist = cDist
                        elif cDist < minDist:
                            minAOT = cAOT
                            minDist = cDist
                    #predAOTArgs = (MinB2RAD[i], MeanB2DOS[i], aeroProfile, atmosProfile, grdRefl, MeanElev[i]/1000)
                    #res = minimize(self.run6SToOptimiseAODValue, minAOT, method='nelder-mead', options={'maxiter': 20, 'xtol': 0.001, 'disp': True}, args=predAOTArgs)
                    #aotVals[i] = res.x[0]
                    aotVals[i] = minAOT
                    print("IDENTIFIED AOT: ", aotVals[i])
                else:
                    aotVals[i] = 0
            rat.writeColumn(ratDS, "AOT", aotVals)

            Eastings = rat.readColumn(ratDS, "Eastings")
            Northings = rat.readColumn(ratDS, "Northings")
            ratDS = None

            Eastings = Eastings[PredictAOTFor!=0]
            Northings = Northings[PredictAOTFor!=0]
            aotVals = aotVals[PredictAOTFor!=0]

            interpSmoothing = 10.0
            self.interpolateImageFromPointData(inputTOAImage, Eastings, Northings, aotVals, outputAOTImage, outFormat, interpSmoothing, True, 0.05)

            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName(outFormat)
                gdalDriver.Delete(thresImageClumpsFinal)
                gdalDriver.Delete(dosBlueImage)

            return outputAOTImage
        except Exception as e:
            raise e

    def estimateSingleAOTFromDOS(self, radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl):
        try:
            return self.estimateSingleAOTFromDOSBandImpl(radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl, 2)
        except Exception as e:
            raise e

    def setBandNames(self, imageFile):
        dataset = gdal.Open(imageFile, gdal.GA_Update)
        for i in range(242):
            dataset.GetRasterBand(i+1).SetDescription("Band%s"%str(i+1))

    def cleanLocalFollowProcessing(self):
        print("")



