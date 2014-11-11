

setwd('/data/jag/dpustina/APHASIA/GroupAnalyses_PerMalToPen025Thickness/thickness/')
file.copy('/data/jag/dpustina/templates/pennTemplate/templateBrain.nii.gz', 'templateBrain.nii.gz')
cores=4
imglist <- Sys.glob('./*PerMal*.nii.gz')
library(foreach)
library(doParallel)
#setup parallel backend to use 4 processor
cl <- makeCluster(cores)
registerDoParallel(cl)

foreach (i=1:length(imglist)) %dopar% {

    newimg <- gsub('./', './smo_3_', imglist[i])
    SmoothImage('3', imglist[i], '3', newimg)

}
stopCluster(cl)

thickaverage <- AverageImages('/data/jag/dpustina/APHASIA/GroupAnalyses_PerMalToPen025Thickness/thickness/smo_3_*PerMal*.nii.gz')
antsImageWrite(thickaverage, 'smo_3_Average.nii.gz')

thickaverage <- antsImageRead('smo_3_Average.nii.gz', 3)
template <- antsImageRead('templateBrain.nii.gz', 3)
mask <- antsImageClone(thickaverage)
ThresholdImage("3", thickaverage, mask, '0.5', '50')
mask[template==0] <- 0
antsImageWrite(mask, 'smo_3_groupmask.nii.gz')
#AverageImages('3', 'smo_3_Average.nii.gz', '0', 'smo_3_*PerMal*.nii.gz')




# load data
library(XLConnect)

# load excel sheet and data for strokes
wk = loadWorkbook("../../newdorianIMAGING_SCAN_BRAKEDOWN found dicoms.xlsm") 
df = readWorksheet(wk, sheet="Sheet1", endRow=139, forceConversion=TRUE)
df = subset(df, select=c("NameOnDisk", "Gender", "Ed", "DicomDate", "AgeatMRIval", "MonthsStrokeToMrival"), subset=(CortThick.Round.1 == 1))
df$DicomDate = as.Date(df$DicomDate)
Strokdata <- df
Strokdata$Gender <- as.factor(Strokdata$Gender)
rownames(Strokdata) <- Strokdata$NameOnDisk

SimageList <- Sys.glob("smo_3_Strok*PerMal*.nii.gz")  # get the filenames
imageList <- rep("", nrow(Strokdata))
# order files same as in excel sheet
for (i in 1:nrow(Strokdata)) {
    subject <- Strokdata[i,]$NameOnDisk
    index <- grep(paste("_",  subject, "_", sep=''), SimageList)
    imageList[i] <- SimageList[index]
}


# load excel sheet and data for controls
wk = loadWorkbook("../../70ContrMatch55Stroke.xlsx") 
df = readWorksheet(wk, sheet="Sheet1", endRow=75, forceConversion=TRUE)
df = subset(df, select=c("INDDID", "Sex", "Education", "MRIDate", "AgeatMRI"), subset=(subjects68 == 1))
df$DicomDate = as.Date(df$DicomDate)
Contrdata <- df
Contrdata$Sex <- as.factor(Contrdata$Sex)
rownames(Contrdata) <- Contrdata$INDDID

SimageList <- Sys.glob("smo_3_Contr*PerMal*.nii.gz")  # get the filenames
imageList <- c(imageList, rep("", nrow(Contrdata)))
# order files same as in excel sheet
for (i in 1:nrow(Contrdata)) {
    subject <- Contrdata[i,]$INDDID
    index <- grep(paste("_",  subject, "_", sep=''), SimageList)
    imageList[nrow(Strokdata) + i] <- SimageList[index]
}

# WE GOT THE DATA AND FILENAMES, LOAD AND RUN GLM
mask <- antsImageRead("smo_3_groupmask.nii.gz",3)  # get the mask

rightmask <- as.array(antsImageClone(mask))
rightmask[82:192, , ] <- 0

mask[mask > -1] <- rightmask*as.array(mask)  # merge mask with rightmask to remove left hemisphere and save time

mat <- imagesToMatrix(imageList,mask)

group <- as.factor(c(rep('Stroke', nrow(Strokdata)), rep('Control', nrow(Contrdata))))

glm <- lm(mat ~ group)
glmstat <- bigLMStats( glm )


thick_pval <- p.adjust(glmstat$beta.pval, method = "fdr")
pval <- antsImageClone(mask)
pval[mask==1] <- thick_pval
pval[pval > -1] <- as.array(pval) * rightmask
antsImageWrite(pval,"smo_3_thick_pval.nii.gz")


thick_beta <- glmstat$beta.t
betas <- antsImageClone(mask)
betas[mask==1] <- thick_beta
betas[pval > 0.05] <- 0  # threshold betas
betas[pval > -1] <- as.array(betas) * rightmask
antsImageWrite(betas,"smo_3_thick_betas.nii.gz")

# ADD MORE VARIABLES THAN JUST GROUP PERTINENCE
commonsex <- as.factor(c(as.character(Strokdata$Gender), as.character(Contrdata$Sex)))  # gender factor of both groups
commonage <- as.numeric(c(Strokdata$AgeatMRIval, Contrdata$AgeatMRI))  # common age factor of both groups
commonstrokeinterval <- as.numeric(c(Strokdata$MonthsStrokeToMrival, rep(0, nrow(Contrdata))))  # common stroke interval of both groups, 0 for controls

glm <- lm(mat ~ group*commonage)
glmstat <- bigLMStats( glm )

rownames(glmstat$beta.t)

groupbetas <- glmstat$beta.t[ "groupStroke" , ]
groupbetas[p.adjust(glmstat$beta.pval[ "groupStroke" , ], method = "fdr") > 0.05] <- 0
agebetas <- glmstat$beta.t["commonage", ]
agebetas[p.adjust(glmstat$beta.pval[ "commonage" , ], method = "fdr") > 0.05] <- 0
ageXgroupbetas <- glmstat$beta.t[ "groupStroke:commonage" , ]
ageXgroupbetas[p.adjust(glmstat$beta.pval[ "groupStroke:commonage" , ], method = "fdr") > 0.05] <- 0

# group
thisimg <- antsImageClone(mask)
thisimg[mask==1] <- groupbetas
thisimg[mask > -1] <- as.array(thisimg) * rightmask
antsImageWrite(thisimg, "smo_3_thick_groupbetas.nii.gz")

# age
thisimg <- antsImageClone(mask)
thisimg[mask==1] <- agebetas
thisimg[mask > -1] <- as.array(thisimg) * rightmask
antsImageWrite(thisimg, "smo_3_thick_agebetas.nii.gz")

# age x group
thisimg <- antsImageClone(mask)
thisimg[mask==1] <- ageXgroupbetas
thisimg[mask > -1] <- as.array(thisimg) * rightmask
antsImageWrite(thisimg, "smo_3_thick_agexgroupbetas.nii.gz")

# test gender pvalgender <- chisq.test(table(commonsex, group))$p.value


# voxel with highest correlation with age
minbeta <- which.min(agebetas)
plot(mat[ , minbeta] ~ commonage,
    main="Voxel with highest correlation with age",
    xlab="Age",
    col=c("red","blue")[group],
    pch=c(1,2)[group])
    
    legend( "topright", legend = levels(group), text.col = c("red", "blue"), col = c("red", "blue"), pt.bg = c("red","blue"), pch = c(1,2) )
    
    
savePlot(filename="LowestBetaAge.png", type="png")

# Average of all voxels correlating with age
temp <- rowMeans(mat[ ,agebetas < -1 ])
plot(temp ~ commonage,
    main="Average of all voxels correlating with age",
    xlab="Age",
    ylab="Average Smoothed Thickness",
    col=c("red","blue")[group],
    pch=c(1,2)[group])
    
    legend( "topright", legend = levels(group), text.col = c("red", "blue"), col = c("red", "blue"), pt.bg = c("red","blue"), pch = c(1,2) )
    
    savePlot(filename="AverageBetaAge.png", type="png")
    
    
    
    
# NOW LETS DO ANOTHER GLM TO INCLUDE MONTHS POST STROKE
glm <- lm(mat ~ commonage*commonstrokeinterval)
glmstat <- bigLMStats( glm )


intervalbetas <- glmstat$beta.t[ "commonstrokeinterval" , ]
intervalbetas[p.adjust(glmstat$beta.pval[ "commonstrokeinterval" , ], method = "fdr") > 0.05] <- 0
intervalXagebeteas <- glmstat$beta.t[ "commonage:commonstrokeinterval" , ]
intervalXagebeteas[p.adjust(glmstat$beta.pval[ "commonage:commonstrokeinterval" , ], method = "fdr") > 0.05] <- 0
agebetas <- glmstat$beta.t[ "commonage" , ]
agebetas[p.adjust(glmstat$beta.pval[ "commonage" , ], method = "fdr") > 0.05] <- 0

# interval
thisimg <- antsImageClone(mask)
thisimg[mask==1] <- intervalbetas
thisimg[mask > -1] <- as.array(thisimg) * rightmask
antsImageWrite(thisimg, "smo_3_glm3_intervalbetas.nii.gz")


# interval X age
thisimg <- antsImageClone(mask)
thisimg[mask==1] <- intervalXagebeteas
thisimg[mask > -1] <- as.array(thisimg) * rightmask
antsImageWrite(thisimg, "smo_3_glm3_intervalXagebetas.nii.gz")

# interval X age
thisimg <- antsImageClone(mask)
thisimg[mask==1] <- agebetas
thisimg[mask > -1] <- as.array(thisimg) * rightmask
antsImageWrite(thisimg, "smo_3_glm3_agebetas.nii.gz")
