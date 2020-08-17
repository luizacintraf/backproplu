#abrir bibliotecas
library(raster) 
library(dismo)

#funcao para escolher aleatoriamente pixels em janelas definidas, entrada: tam janela(div), tamanho máximo da largura(size)
index<-function(size,div){
    ind_i<-c(rep(NA,(round(size/div))))
    c<-1
    for (i in seq(1,(size-div),div)){
        ind_i[c]<-sample(i:(i+div),1)
        c<-c+1
    }
    ind_i[!is.na(ind_i)]->ind_i
    ind_i
}

path='PATH_TO_FILES' #local dos dados

numvar<-14 #numero de entradas


#matriz com o caminho dos arquivos da rede
f<-matrix(0,nrow=18,ncol=(numvar+2)) 


#fazendo a matriz f
c<-1
for (i in 2014:2016){
    for(j in 1:4){
        lista<-list.files(file.path(path,i,j))
        for (n in 1:length(lista)){
            f[c,n]<-file.path(path,i,j,lista[n])
        }
        c<-c+1
    }
    for(j in 11:12){
        lista<-list.files(file.path(path,i,j))
        for (n in 1:length(lista)){
           f[c,n]<-file.path(path,i,j,lista[n])
        }
       c<-c+1
    }
    
}



#arrumando a ordem dos dados na tabela
t<-grep('fire',f[1,])
f[,t]->f[,dim(f)[2]]
f[,t:(dim(f)[2]-1)]<-f[,(t+1):dim(f)[2]]
f[,1:(dim(f)[2]-1)]->f

d<-grep('direcao',f[1,])  
f[,d:(dim(f)[2]-1)]<-f[,(d+1):dim(f)[2]]
f[,1:(dim(f)[2]-1)]->f

d<-grep('usodosolo_bh',f[1,])
f[,d:(dim(f)[2]-1)]<-f[,(d+1):dim(f)[2]]
f[,1:(dim(f)[2]-1)]->f

d<-grep('vu',f[1,])
f[,d:(dim(f)[2]-1)]<-f[,(d+1):dim(f)[2]]
f[,1:(dim(f)[2]-1)]->f

d<-grep('aspect',f[1,])
f[,d:(dim(f)[2]-1)]<-f[,(d+1):dim(f)[2]]
f[,1:(dim(f)[2]-1)]->f


d<-grep('slope',f[1,])
f[,d:(dim(f)[2]-1)]<-f[,(d+1):dim(f)[2]]
f[,1:(dim(f)[2]-1)]->f

d<-grep('vp',f[1,])
f[,d:(dim(f)[2]-1)]<-f[,(d+1):dim(f)[2]]
f[,1:(dim(f)[2]-1)]->f

d<-grep('disturb',f[1,])
f[,d:(dim(f)[2]-1)]<-f[,(d+1):dim(f)[2]]
f[,1:(dim(f)[2]-1)]->f

d<-grep('UR',f[1,])
f[,d:(dim(f)[2]-1)]<-f[,(d+1):dim(f)[2]]
f[,1:(dim(f)[2]-1)]->f

d<-grep('tempmedia',f[1,])
f[,d:(dim(f)[2]-1)]<-f[,(d+1):dim(f)[2]]
f[,1:(dim(f)[2]-1)]->f

#abrir tabela dados de fogo
table<-read.table('Path_fire_table', sep=';', header=TRUE) 
#espacializar os dados
coordinates(table)<-table[,2:1]

#tamanho da janela (pixels)
div<-130
#vetor nulo
data<-c() 

for (k in 1:dim(f)[1]){ #quantidade de meses
    ano=as.integer(strsplit(f[k,1],"/")[[1]][8]) #pegar os anos
    tab<-table[table$ano==ano,] #tabela fogo ano
    print(paste0('k:  ', k)) #controle de mes
    ta<-tab[tab$mes==as.integer(strsplit(f[k,1],"/")[[1]][9]),] #tabela fogo por mes
    mes=as.integer(strsplit(f[k,1],"/")[[1]][9])
    fogo<-raster(f[k,dim(f)[2]]) #arquivo de fogo
    #fogo<-crop(fogo,rmbh)
    larg<-nrow(fogo) #largura arquivo fogo
    altu<-ncol(fogo) #altura arquivo fogo
    ind_i<-index(larg,div) #escolher aleatoriamente indices i
    ind_j<-index(altu,div)#escolher aleatoriamente indices j
    
    #fire<-grep(1,values(fogo)) #verificar locais de pixels iguais a 1 (fogo efetivo)
    cellFromRowColCombine(fogo, ind_i, ind_j)->ind
    xyFromCell(fogo,ind)->cor
    fogo<-fogo[ind_i,ind_j] #selecionar pixels nao fogo
    firev<-rep(1.0,length(ta)) #dados  de fogo=1
    fogo<-append(fogo,firev)#juntando arquivos fogo
    for(j in 1:(dim(f)[2]-1)){ #variaveis
        xn<-raster(f[k,j]) #puxando a variavel
        val<-extract(xn,cor) #selecionando os pixels nao fogo
        v2<-extract(xn,coordinates(ta)) #selecionando os pixels fogo
        val<-append(val,v2) #juntando arquivos variavel
        fogo<-cbind(val,fogo) #juntando variavel com fogo
        colnames(fogo)[1]<-names(xn) #nomeando a coluna
        print(paste0('j:  ', j))
    }
        
    fogo[rowSums(is.na(fogo))!=14,]->fogo # retirando NA's
    fogo[rowSums(is.na(fogo))==0,]->fogo # retirando NA's
    data<-rbind(data,fogo) #juntando os dados
}

#ajeitando alguns dados
n<-c("vu","vp","disturb","radiacao")
ind<-grep(paste(n,collapse="|"),colnames(data))
data[,ind]<-data[,ind]/1000

ind<-grep("pressao",colnames(data))
data[,ind]<-data[,ind]/10
data[data[,ind]!=0,]->data

ind<-grep("total",colnames(data))
data[data[,ind]<0,]<-0

ind<-grep("radiacao",colnames(data))
data[data[,ind]>0,]->data



#dados de entrada
x<-data[,1:(dim(data)[2]-1)] 

#dados de saida
y<-data[,dim(data)[2]] 
y[y>0.5]<-1
y[y<0.5]<-0


#sorteio de pixels
nfolds<-3
kf<-kfold(data,nfolds) 

dtest<-data[kf==1,] #selecionar pixels teste dados
dtrain<-data[kf!=1,] #selecionar pixels treinamento dados

lf<-kfold(dtest,2) #sorteio validacao e generalizacao
dvalidate<-dtest[lf==1,] #selecionar pixels validacao dados
dgeneral<-dtest[lf==2,] #selecionar pixels geeralizacao dados

#dados de treinamento
x<-dtrain[,1:(dim(dtrain)[2]-1)] #dados de entrada
y<-dtrain[,dim(dtrain)[2]]
y[y>0.5]<-1
y[y<0.5]<-0


#dados de validação
x_v<-dvalidate[,1:(dim(dvalidate)[2]-1)] #dados de entrada
y_v<-dvalidate[,dim(dvalidate)[2]]
y_v[y_v>0.5]<-1
y_v[y_v<0.5]<-0


#dados de generalização
x_g<-dgeneral[,1:(dim(dgeneral)[2]-1)] #dados de entrada
y_g<-dgeneral[,dim(dgeneral)[2]]
y_g[y_g>0.5]<-1
y_g[y_g<0.5]<-0


#escrevendo tabelas
write.table(x,'x.txt', row.names=FALSE,col.names=FALSE) #escrever dados de treinamento
write.table(y,'y.txt', row.names=FALSE,col.names=FALSE) #escrever dados de treinamento de fogo
write.table(x_v,'x_valid.txt', row.names=FALSE,col.names=FALSE) #escrever dados de validacao
write.table(y_v,'y_valid.txt', row.names=FALSE,col.names=FALSE) #escrever dados de validacao de fogo
write.table(x_g,'x_gen.txt', row.names=FALSE,col.names=FALSE) #escrever dados de generalizacao
write.table(y_g,'y_gen.txt', row.names=FALSE,col.names=FALSE) #escrever dados de generalizacao de fogo
