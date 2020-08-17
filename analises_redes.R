library(NeuralNetTools)
library(nnet)
source('https://gist.githubusercontent.com/fawda123/6860630/raw/b8bf4a6c88d6b392b1bfa6ef24759ae98f31877c/lek_fun.r')
library(devtools)
source_gist('6206737')
source_url('https://gist.github.com/fawda123/7471137/raw/cd6e6a0b0bdb4e065c597e52165e5ac887f5fe95/nnet_plot_update.r')


#Ler arquivos de entrada e saída
x<-read.table("arquivo_entradas")
y<-read.table("arquivo_saidas")

#ler aqrquivo de pesos
wh1<-as.matrix(read.table(file.path(pathpeso,'wh1.txt'),sep=' '))
ws<-as.numeric(read.table(file.path(pathpeso,'ws.txt'),sep=' ')[,1])
bh1<-as.numeric(read.table(file.path(pathpeso,'bh1.txt'),sep=' ')[,1])
bs<-as.numeric(read.table(file.path(pathpeso,'bs.txt'),sep=' ')[,1])

#arrumar os dados
cbind(y,x)->data
colnames(data)[1]<-'Y1'
#normalizar os dados
normalize(data)[[1]]->data

#crriando modelo da rede Neural
mod <- nnet(Y1 ~ V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12, data = data, size = 6)
wts<-c(bh1[1],wh1[,1],bh1[2],wh1[,2],bh1[3],wh1[,3],bh1[4],wh1[,4],bh1[5],wh1[,5],bh1[6],wh1[,6],bs,ws)
mod$wts<-wts


#testando a rede pelo método Garson
jpeg("testeimportancia.jpg", width = 25, height = 10, units = 'in', res = 400)
teste<-t(garson(mod,x_lab=c("temp_max","UR_min","dist_rod","declividade","dist_urb", "orientação","uso do solo","vento_medio","prec_total", "radiação","pressão", "NDVI"),bar_plot=FALSE))
pal=colorRampPalette(colors=c("darkblue","lightblue"))(12)
barplot(teste[order(teste,decreasing=TRUE)],names.arg=colnames(teste)[order(teste,decreasing=TRUE)],col=pal,ylab='Importância Relativa',ylim=c(0,0.3),cex.axis=1.6,cex.names=1.6,axes=TRUE,border=TRUE,axis.lty = 1,cex.lab=1.5)
dev.off()

#Criando grafico da Rede Neural
jpeg("RNA.jpg", width = 10, height = 10, units = 'in', res = 300)
plot.nnet(mod,x.lab=c("temp_max (oC)","UR_min (%) ","dist_rod (km)","declividade (%)","dist_urb (km)", "orientação (graus)","uso do solo","vento_medio (m/s)","prec_total (mm)", "radiação (MJ/m2)","pressão (kPa)", "NDVI"))
dev.off()
