# calling libraries
library(raster)
library(rgdal)
library(prettymapr)

#---------------------------------------------#
#             Calling functions
#---------------------------------------------#

# activation functions
activation<-function(v, afunction){
    a<-1
    if(afunction==1){ #logistic
        y <- 1./(1.+ exp(-a * v))
    }else if(afunction==2){ #tanget
        y <- (1. - exp(-v))/(1. + exp(-v))
        #y<-tanh(v)
    }else if(afunction==3){ #gauss
        y <- exp(-v)
    }
    y
}

# normalize data
normalize<-function(x){
    maxi<-max(x,na.rm=TRUE)
    mini<-min(x,na.rm=TRUE)
    x<-(x-mini)/(maxi-mini)
    x
}


#----------------------------------------------------------#
#                       ARRANGE DATA                       #
#----------------------------------------------------------#

#path for folder with activation files
pathfiles<-"PATH_TO_ACTIVATION_FILES"

#path for boundary of study area
limite<-readOGR('PATH_STUDY_AREA_BOUNDARY')

#list od files in folder
lista<-list.files(pathfiles)

#read fire data 
f2017<-readOGR('PATH_TO_FILE')

#call trained weights
pathpeso<-"PATH_TO_WEIGHTS"
wh1<-as.matrix(read.table(file.path(pathpeso,'wh1.txt'),sep=' '))
ws<-as.numeric(read.table(file.path(pathpeso,'ws.txt'),sep=' ')[,1])
bh1<-as.numeric(read.table(file.path(pathpeso,'bh1.txt'),sep=' ')[,1])
bs<-as.numeric(read.table(file.path(pathpeso,'bs.txt'),sep=' ')[,1])

#create percentage matrix
percentage<-matrix(0,nrow=12,ncol=6)

#----------------------------------------------------------#
#                START RNA                                 #
#----------------------------------------------------------#


#configurations
nInputs<-11
nOutputs<-1
neuronsLayer1<- 6
activationFunction<-2

#calculate the resut for each month
for ( me in 1:12){


    #grep the files of the month
    mest<-paste0('_',me,'.tif')
    lista1<-grep(mest,lista,value=TRUE)


    #call rasters for each variable
    x1<-raster(file.path(pathfiles,grep('tempmax',lista1,value=TRUE)))
    x2<-raster(file.path(pathfiles,grep(paste0('UR_',me,'.tif'),lista1,value=TRUE)))
    x3<-raster(file.path(pathfiles,grep('vp',lista,value=TRUE)))
    x4<-raster(file.path(pathfiles,grep('slope',lista,value=TRUE)))
    x5<-raster(file.path(pathfiles,grep('disturb',lista,value=TRUE)))
    x6<-raster(file.path(pathfiles,grep('aspect',lista,value=TRUE)))
    x7<-raster(file.path(pathfiles,grep('uso',lista,value=TRUE)))
    x8<-raster(file.path(pathfiles,grep('vento',lista1,value=TRUE)))
    x9<-raster(file.path(pathfiles,grep('total',lista1,value=TRUE)))
    x10<-raster(file.path(pathfiles,grep('rad',lista1,value=TRUE)))
    x11<-raster(file.path(pathfiles,grep('pressao',lista1,value=TRUE)))
    x12<-raster(file.path(pathfiles,grep('ndvi',lista1,value=TRUE)))


    #crop all files with boundary to make sure that they have the same size

    x2<-crop(x2,limite)
    x3<-crop(x3,limite)
    #writeRaster(x3,file.path(pathfiles,grep('vp',lista,value=TRUE)),overwrite=TRUE)
    x4<-crop(x4,limite)
    #writeRaster(x4,file.path(pathfiles,grep('slope',lista,value=TRUE)),overwrite=TRUE)
    x5<-crop(x5,limite)
    #writeRaster(x5,file.path(pathfiles,grep('disturb',lista,value=TRUE)),overwrite=TRUE)
    x6<-crop(x6,limite)
    #writeRaster(x6,file.path(pathfiles,grep('aspect',lista,value=TRUE)),overwrite=TRUE)
    #x7<-crop(x7,limite)
    #x8<-crop(x8,limite)
    #x9<-crop(x9,limite)
    #x10<-crop(x10,limite)
    #x11<-crop(x11,limite)
    x12<-crop(x12,limite)
    #writeRaster(x12,file.path(pathfiles,grep('ndvi',lista1,value=TRUE)),overwrite=TRUE)


    #put all files to a list, and change each variable to matrix 
    x<-list(as.matrix(x1),as.matrix(x2),as.matrix(x3),as.matrix(x4),as.matrix(x5),as.matrix(x6),as.matrix(x7),as.matrix(x8),as.matrix(x9), as.matrix(x10), as.matrix(x11),as.matrix(x12))

    #normalize all the data
    for (i in 1:12){
        x[[i]]<-normalize(x[[i]])
    }

    #create empty variables to be use
    s<-matrix(0,nrow=dim(x1)[1],ncol=dim(x1)[2])
    vh1<-list()
    yh1<-list()


    #run the RNA
    for ( j in 1:neuronsLayer1){
        for (i in 1:nInputs){
            m<-x[[i]]*wh1[i,j]
            s<-s+m
        }
        vh1[[j]]<-s-bh1[j]
        yh1[[j]]<-activation(vh1[[j]],activationFunction)
        s<-matrix(0,nrow=dim(x1)[1],ncol=dim(x1)[2])
    }

    s<-matrix(0,nrow=dim(x1)[1],ncol=dim(x1)[2])
    for ( j in 1:neuronsLayer1){
            m<-yh1[[j]]*ws[j]
            s<-s+m
    }
    vs<-s-bs
    ys<- activation(vs, activationFunction)

    #create a raster with the output
    fogo<-raster(ys)
    extent(fogo)<-extent(x1)
    projection(fogo)<-projection(x1)

    #get the fire points from the selected month
    date<-as.Date(f2017$ACQ_DATE)
    mes<- format(date,'%m')
    fire<-f2017[grep(sprintf("%02d", me),mes),]
    fmes<-extract(fogo,coordinates(fire))
    fmes<-fmes[!is.na(fmes)]

    #Get the quantities of each categories
    VH<-length(fmes[fmes>0.8])
    H<-length(fmes[fmes>0.6])-length(fmes[fmes>0.8])
    M<-length(fmes[fmes>0.4])-length(fmes[fmes>0.6])
    L<-length(fmes[fmes>0.2])-length(fmes[fmes>0.4])
    VL<-length(fmes[fmes<0.2])

    #Get the percentage of each categories
    pVH<-VH/length(fmes)
    pH<-H/length(fmes)
    pM<-M/length(fmes)
    pL<-L/length(fmes)
    pVL<-VL/length(fmes)

    #update percentage table
    percentage[me,1]<-me
    percentage[me,2]<-pVH
    percentage[me,3]<-pH
    percentage[me,4]<-pM
    percentage[me,5]<-pL
    percentage[me,6]<-pVL

    #write the tiff of raster
    writeRaster(fogo,paste0("modelo_",me,".tif"),overwrite=TRUE)

    #create a image of the result
    jpeg(paste0("modelo_",me,".jpeg"), width = 6 , height = 8, units = 'in', res = 300)

    mini<-min(values(fogo),na.rm=TRUE)

    #reclassify the values
    if(mini<0){
    reclass_df <- c(mini, 0.2, 1, 0.2, 0.4, 2, 0.4, 0.6, 3, 0.6, 0.8, 4, 0.8, 1, 5)
    }else{
    reclass_df <- c(0, 0.2, 1, 0.2, 0.4, 2, 0.4, 0.6, 3, 0.6, 0.8, 4, 0.8, 1, 5)
    }
    reclass_m <- matrix(reclass_df, ncol = 3, byrow = TRUE)
    chm_classified <- reclassify(fogo,reclass_m)
    plot(chm_classified, legend = FALSE, col = c("dark green", "green", "yellow","orange","red"), axes = TRUE, main = "Risco de Incêndios")

    #create legend
    legend("topleft", legend = c("muito baixo","baixo", "médio","alto", "muito alto","focos de incêndio"),
    fill = c("dark green", "green", "yellow","orange", "red", "black"), border = FALSE, bty = "n")
    #create scale
    scalebar(40, xy=c( -43.9,-20.55), type = "bar", divs = 4, below = "km", lonlat = TRUE,  adj=c(0.5, -1.3), lwd = 2,cex=0.8)
    #create northarrow
    addnortharrow(scale = 0.6, text.col = 'black', cols = c('black', 'black'),pos="topright")

    plot(fire, add=TRUE, col='black',pch=16)

    dev.off()

    colnames(percentage)<-c("mes","Muito Alto","Alto","Médio","Baixo", "Muito Baixo")

}

#write the final table
write.table(percentage,'porcentagem.txt',col.names=TRUE,row.names=FALSE)