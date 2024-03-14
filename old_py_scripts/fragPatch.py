#PatchMorph Code, modified to FragPatch algorithm (see comments below)
#Yellow-highlighted text is where the user needs to change the text to fit their data/needs

import arcpy
from arcpy import env
from arcpy.sa import *

arcpy.CheckOutExtension("spatial")

arcpy.OverWriteOutput = True

Input_raster = "MogllonRimCCEx.tif" #0 = habitat, 100 = non-habitat
Change_missing_values_to_NoData = "false"

suitLow = 0
suitHigh = 100
suitNum = 2
suitBy = (suitHigh - suitLow) / suitNum

gapLow = 3 # number in map units
gapHigh = 15 # number in map units
gapNum = 5 # How many gaps do you have between Min. and Max.?
gapBy = (gapHigh - gapLow) / gapNum

spurLow = 5 # number in map units
spurHigh = 20 # number in map units
spurNum = 4 # How many spurs do you have between Min. and Max.?
spurBy = (spurHigh - spurLow) / spurNum

outDir = "D:\\fragPatch\\output"
arcpy.Workspace = outDir

suitList = range(suitLow, suitHigh, suitBy)
gapList = range(gapLow, gapHigh, gapBy)
spurList = range(spurLow, spurHigh, spurBy)

for cntSuit in range(2, suitNum+1):
	suit = suitList[cntSuit-1]
	arcpy.CreateFolder_management(outDir, "\\suit" + str(suit))
	outDir1 = outDir + "\\suit" + str(suit)
# In PatchMorph, the next few lines reference “Gaps”. In FragPatch, we changed the next few
# lines to reference “Spurs”.
	for cntSpur in range(1, spurNum+1):
		spurDist = spurList[cntSpur-1]
		arcpy.CreateFolder_management(outDir1, "\\spur" + str(spurDist))
		outDir2 = outDir1 + "\\spur" + str(spurDist)
 
# In PatchMorph, the next few lines reference “Spurs”. In FragPatch, we changed the next few
# lines to reference “Gaps”.
		for cntGap in range(1, gapNum+1):
			gapDist = gapList[cntGap-1]
			arcpy.CreateFolder_management(outDir2, "\\gap" + str(gapDist))
			outDir3 = outDir2 + "\\gap" + str(gapDist)
			suitRemap = RemapRange([[suit, suitHigh, 1]])
			temp1 = Reclassify(Input_raster, "VALUE", suitRemap, 
				"NODATA")
			temp4 = EucDistance(temp1, "", "", "")
			temp5 = FocalStatistics(temp4, NbrCircle(spurDist/2, "MAP"),
			"MAXIMUM", "DATA")
			spurRemap = RemapRange([[spurDist/2, 99999999999, 1]])
			temp6 = Reclassify(temp5, "VALUE", spurRemap, "NODATA")
			temp7 = EucDistance(temp6, "", "", "")
			temp8 = FocalStatistics(temp7, NbrCircle(gapDist/2, "MAP"),
			"MAXIMUM", "DATA")
			gapRemap = RemapRange([[gapDist/2, 99999999999, 1]])
			outputRaster = Reclassify(temp8, "VALUE", gapRemap, 
				"NODATA")
			outputRaster.save(outDir3 + "\\su" + str(cntSuit) + "sp" + 
				str(cntSpur) + "gp" + str(cntSpur))
			remapClass = RemapValue([[1,1],["NoData",0]])
			temp9 = Reclassify(outputRaster, "VALUE", remapClass, 
				"DATA")
			if ((cntGap == 1) & (cntSpur == 1)):
				sumRastSuiti = Plus(temp9, 0)
			else:
				sumRastSuiti = Plus(sumRastSuiti, temp9)
	sumRastSuiti.save(outDir1 + "\\psum_su" + str(cntSuit))
	if (cntSuit == 2):
		sumRast = Plus(sumRastSuiti, 0)
	else:
		sumRast = Plus(sumRast, sumRastSuiti)
sumRast.save(outDir + "\\patch_sum")

