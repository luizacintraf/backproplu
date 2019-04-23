
#--------------------------------------------------------------------#
#FUNCOES DA RNA
#--------------------------------------------------------------------#

#FUNCAO DE ATIVACAO
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

#GRADIENTE DA FUNCAO DE ATIVACAO
derivate<-function(v,y, afunction){
    a<-1
    if(afunction==1){ #logistic
        d <- ((a * exp(-a * v))/((1. + exp(-a * v))^2.))
    }else if(afunction==2){ #tanget
        d <- (2 * exp(-v)) / ((1. + (exp(-v))^2))
        #d<-1-(tanh(v))^2
    }else if(afunction==3){ #gauss
        d <- -y/a
    }
    d
}

#NORMALIZACAO DOS DADOS
normalize<-function(x){
    maxi<-c(rep(0,dim(x)[2]))
    mini<-c(rep(0,dim(x)[2]))
    for (i in 1:dim(x)[2]){
        maxi[i]<-max(x[,i])
        mini[i]<-min(x[,i])
        x[,i]<-(x[,i]-mini[i])/(maxi[i]-mini[i])
    }
    list(x,maxi,mini)
}


#--------------------------------------------------------------------#
#CARREGAR CONJUNTO DE DADOS
#--------------------------------------------------------------------#

#DADOS DE TREINAMENTO
x<-as.matrix(read.table('x.txt'))
y<-as.numeric(read.table('y.txt')[,1])

#DADOS DE VALIDAÇÃO
x_v<-as.matrix(read.table('x_valid.txt'))
y_v<-as.numeric(read.table('y_valid.txt')[,1])

#DADOS DE GENERALIZAÇÃO
x_g<-as.matrix(read.table('x_gen.txt'))
y_g<-as.numeric(read.table('y_gen.txt')[,1])

#--------------------------------------------------------------------#
#NORMALZIAÇÃO CONJUNTO DE DADOS
#--------------------------------------------------------------------#

x<-normalize(x)[[1]]
x<-t(x)

x_v<-normalize(x_v)[[1]]
x_v<-t(x_v)

x_g<-normalize(x_g)[[1]]
x_g<-t(x_g)

#--------------------------------------------------------------------#
#PARÂMETROS DA REDE
#--------------------------------------------------------------------#

nInputs<-12 #ENTRADAS
nOutputs<-1 #SAÍDAS
activationFunction<-2  #FUNÇÃO DE ATIVAÇÃO
#1 logistic
#2 tangent
#3 gauss
neuronsLayer1<- 6 #NEURONIOS NA CAMADA ESCONDIDA
eta<-  0.0674  #TAXA DE APRENDIZAGEM
alpha<-0.6304  #CONSTANTE DE MOMENTO

epocaN<-200 #NÚMERO DE ÉPOCAS
errod<-1E-5 #ERRO DESEJADO


nClasses<-dim(x)[2]
nClassesValidation<-dim(x_v)[2]
nClassesGeneral<-dim(x_g)[2]


#--------------------------------------------------------------------#
#INICIALIZAÇÃO PESOS E BIAS
#--------------------------------------------------------------------#


deltaWeightOutput <- 0
deltaWeightHiddenLayer1 <- matrix(0,nrow=nInputs,ncol=neuronsLayer1)
deltaBiasOutput <- 0
deltaBiasHiddenLayer1 <-rep(0,neuronsLayer1)


wh1<-matrix(rnorm(n=nInputs*neuronsLayer1,mean=0,sd=((nInputs)^(-1/2))),nrow=nInputs,ncol=neuronsLayer1) #inicializando pesos da camada de entrada
ws<-rnorm(neuronsLayer1,mean=0,sd=((neuronsLayer1)^(-1/2))) #inicializando pesos da camada oculta
bh1<-rnorm(n=neuronsLayer1,mean=0,sd=((nInputs)^(-1/2))) #inicializando bias da camada oculta
bs<-rnorm(n=nOutputs,mean=0,sd=((neuronsLayer1)^(-1/2)))  #inicializando bias da camada de saida


#----------------------------------------------------------------------#
# FEEDFORWARD: COMEÇO DA RNA
#----------------------------------------------------------------------#

#ZERANDO PARAMETROS
epoch <- 0
error<-rep(0,nClasses)
Erro<-0
errorClass<-rep(0,nClasses)
gradientHiddenLayer1<-rep(0,neuronsLayer1)
MeanSquaredError<-1

#EMBARALHANDO DADOS
sample(1:nClasses,nClasses)->ind
x[,ind]->x
y[ind]->y

#TEMPO INICIAL
t1<-Sys.time()

while (MeanSquaredError>=errod && epoch<epocaN){

epoch = epoch + 1

    for (i in 1:nClasses){
        

        deltaWeightOutputLast <- deltaWeightOutput
        deltaWeightHiddenLayer1Last <- deltaWeightHiddenLayer1
        deltaBiasOutputLast <- deltaBiasOutput
        deltaBiasHiddenLayer1Last <- deltaBiasHiddenLayer1
        
        # ATIVANDO CAMADA ESCONDIDA 1
        vh1 <- (x[,i]%*%wh1) - bh1
        yh1 <- activation(vh1, activationFunction)

        vs <- yh1%*%ws - bs
        
        #ATIVANDO SAIDA
        ys <- activation(vs, activationFunction)

        error[i] <- y[i] - ys

        #CALCULO PADRAO DO ERRO
        errorClass[i] <- 0.5 * (error[i]^2.)
        
        
    
        #-------------------------------------------------------------------------#
        #                        BACKPROPAGATION
        #-------------------------------------------------------------------------#
        #TREINAMENTO CAMADA DE SAÍDA

        dOutput <- derivate(vs, ys, activationFunction)
        
        gradientOutput <- as.numeric(error[i] * dOutput)
        
        deltaWeightOutput <- eta * as.numeric(gradientOutput) * yh1

        deltaBiasOutput<- (-1.) * eta * gradientOutput
        
        ws <-  as.numeric(ws + deltaWeightOutput + alpha * deltaWeightOutputLast)
        
        bs <- as.numeric(bs + deltaBiasOutput + alpha * deltaBiasOutputLast)

        #TREINANDO CAMADA OCULTA 1
        for (j in 1:neuronsLayer1){
            aux <- 0.
            aux <- aux + gradientOutput *  ws
            dOutput <- derivate(vh1, yh1, activationFunction)
            
            gradientHiddenLayer1<- dOutput * aux

            deltaWeightHiddenLayer1[,j] <-eta * gradientHiddenLayer1[j] * x[, i]
            deltaBiasHiddenLayer1[j] <-(-1) * eta * gradientHiddenLayer1[j]
            
            wh1[,j]<-wh1[, j] + deltaWeightHiddenLayer1[,j] + alpha * deltaWeightHiddenLayer1Last[,j]
            bh1[j] <- bh1[j] + deltaBiasHiddenLayer1[j] + alpha * deltaBiasHiddenLayer1Last[j]
        }
    
    }

MeanSquaredError = sum(errorClass) / nClasses
print(paste0('Epoca:  ',epoch))
print(paste0('Error:  ',MeanSquaredError))
Erro<-c(Erro,MeanSquaredError)
t2<-Sys.time()
}


#SALVANDO ERRO
Erro<-Erro[-1]

jpeg("Erro.jpg", width = 7, height = 5, units = 'in', res = 300)
plot(Erro,t="l",ylim=c(min(Erro),max(Erro)),xlab="Número de Épocas", ylab="Erro")
dev.off()

#TEMPO FINAL
t3<-t2-t1

write.table(Erro,'erro.txt',col.names=FALSE,row.names=FALSE) #ESCREVENDO ERRO

#ESCREVENDO PESOS
write.table(wh1,'wh1.txt',col.names=FALSE,row.names=FALSE)
write.table(ws,'ws.txt',col.names=FALSE,row.names=FALSE)
write.table(bh1,'bh1.txt',col.names=FALSE,row.names=FALSE)
write.table(bs,'bs.txt',col.names=FALSE,row.names=FALSE)

#--------------------------------------------------------------------#
#PRODUZIDO POR LUIZA CINTRA FERNANDES
#TREINAMENTO DA REDE PARA PRODUÇÃO DE MODELO DE INCÊNDIOS
#UNIVERSIDADE FEDERAL DE MINAS GERAIS (2019)
#--------------------------------------------------------------------#
