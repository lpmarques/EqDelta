#lê tabela com resultados compilados após todas as otimizações
frame <- read.table("frames+extras.txt.clean",head=T,sep="\t",quote="")

#inclui extensão
load("frame_extension.RData")

frame <- cbind(frame,sum_b/mean_b,sum_b,sumlog_b,mean_b,sd_b,skew_b,kurt_b,quant_b)
colnames(frame) <- c("tree","total_k","tips_number","likelihood","gamma_shape","seq_length","gap_prop","tree_length","rt_RFdist","rt_KFdist","rr_RFdist","rr_KFdist","rt_KFdist_diff","rt_delta","rr_delta","mean_rr_delta","intended_rr_KFdist","tree_k","sum_b","sumlog_b","mean_b","sd_b","skew_b","kurt_b","min_b","quant10_b","quant20_b","quant30_b","quant40_b","median_b","quant60_b","quant70_b","quant80_b","quant90_b","max_b")
error <- error[-1]
frame <- frame[error,]
remove(i,id,kurt_b,mean_b,quant_b,row,sd_b,skew_b,splitid,sum_b,sumlog_b,t,tree1,tree2,trioid,v)

h0 <- as.data.frame(frame)
attach(h0)

#corrige a distância em releção a arvore verdadeira para a média das distâncias do par
h0[,"mean_rt_KFdist"] <- rt_KFdist-(rt_KFdist_diff/2)
detach(h0)
attach(h0)

for(t in seq(30,180,50)){
	h0[which(tips_number==t),"tree_k"] <- t*2-3
}

#====== DAS ANÁLISES DE DISPERSÃO DE DELTA ======#

#há relação entre as variáveis compiladas e a dispersão de delta?
par(mfrow=c(2,3))
plot(gap_prop,mean_rr_delta,cex=0.1,xlab="Proporção de gaps (%)",ylab="∆lnL",main="a",cex.lab=1.4,cex.axis=1.35,cex.main=2)
plot(seq_length,mean_rr_delta,cex=0.1,xlab="Número de sítios",ylab="∆lnL",main="b",cex.lab=1.5,cex.axis=1.35,cex.main=2)
plot(tree_k,mean_rr_delta,cex=0.1,xlab="Número de ramos",ylab="∆lnL",main="c",cex.lab=1.5,cex.axis=1.35,cex.main=2)
plot(gamma_shape,mean_rr_delta,cex=0.1,xlim=c(0.57,0.80),xlab="Alfa de Gama",ylab="∆lnL",main="d",cex.lab=1.5,cex.axis=1.35,cex.main=2)  # seção que corresponde aproximadamente ao hpd 95% empírico de valores de alfa
plot(rr_KFdist,mean_rr_delta,cex=0.1,xlim=c(0,0.3),xlab="BSD alternativa-alternativa",ylab="∆lnL",main="e",cex.lab=1.5,cex.axis=1.35,cex.main=2)  # seção que corresponde aproximadamente ao hpd 95% empírico dos valores de distância entre alternativas
plot(mean_rt_KFdist,mean_rr_delta,cex=0.1,xlim=c(0,0.34),xlab="BSD alternativas-verdadeira",ylab="∆lnL",main="f",cex.lab=1.5,cex.axis=1.35,cex.main=2)  # seção que corresponde aproximadamente ao hpd 95% empírico de valores de distância entre alternativas e sua verdadeira


# Os plots mostram que os valores médios de delta entre árvores equidistantes da verdadeira tendem a se dispersar em torno de 0, não importando alterações em outras variáveis
# essa dispersão, no entanto, parece estar relacionada a gp, sl, rr_dist e rt_dist média. O alfa e o número de tips, por sua vez, não apresentaram nenhum efeito muito evidente sobre ela, a priori
# Como a rt_dist média se mostrou influente apesar de não ser um dos fatores parcial ou totalmente controlados e há uma relação lógica entre ela a rr_dist, fiz plots de uma dist a contra outra para avaliar essa relação.

plot(rr_KFdist,mean_rt_KFdist,cex=0.1,xlab="BSD alternativa-alternativa",ylab="BSD alternativas-verdadeira",main="a",cex.lab=1.5,cex.axis=1.5,cex.main=2)

# Como esperado, o caráter métrico da distância BSD impede que a rt_dist média (em situação de equidistância) seja menor que a metade da rr_dist (vide um triângulo retângulo).
# Ao mesmo tempo, não há limite superior para a rt_dist média, que se extende de forma relativamente independente da rr_dist.
# Para garantir que esta propriedade não estava mascarando o efeito dos outros fatores sobre a dispersão do delta médio, repeti o plot separando diferentes categorias desses outros fatores.
par(mfrow=c(3,2))
for(gp in seq(0,0.2,0.1)){
	w <- which(gap_prop>=gp-0.02&gap_prop<=gp+0.02)
	plot(rr_KFdist[w],mean_rt_KFdist[w],cex=0.1,ylim=c(0,1),main=sprintf("%d%% gp",gp*100))
}
for(gp in seq(0.3,0.4,0.1)){
	w <- which(gap_prop>=gp-0.02&gap_prop<=gp+0.02)
	plot(rr_KFdist[w],mean_rt_KFdist[w],cex=0.1,ylim=c(0,1),main=sprintf("%d%% gp",gp*100))
}
# as rt_dist médias se distribuíram de forma homogênea ao longo das rr_dist para as diferentes gp

par(mfrow=c(3,3))
for(sl in c(100,200,300,400,500,750,1000,1500,2000)){
	w <- which(seq_length==sl)
	plot(rr_KFdist[w],mean_rt_KFdist[w],cex=0.1,xlim=c(0,2),ylim=c(0,2),main=sprintf("%d sl",sl))
}
# para as diferentes sl também, apesar de, aqui, parecerem reduzir de forma relativamente homogênea quanto maior o sl.
# Dado que maior número de observações (sítios) leva a estimações, inclusive de tamanhos de ramos, mais precisas (Swartz & Muller, 2010),
# não é surpreendente que estimações com alinhamentos mais extensos diminuam a distância entre árvores alternativas e a verdadeira.

par(mfrow=c(2,2))
for(tn in seq(30,180,50)){
	w <- which(tips_number==tn)
	plot(rr_KFdist[w],mean_rt_KFdist[w],cex=0.1,xlim=c(0,2),ylim=c(0,2),main=sprintf("%d sl",tn))
}
# Surpreendentemente, as árvores com 130 e 180 taxons tendem a apresentar distâncias muito mais reduzidas de sua árvore verdadeira quanto maior ou menor, respectivamente, é a distância entre elas.
# Essa distorção pode se dever a alguma particularidade dessas árvores e ter mascarado o efeito de seu tamanho sobre a disperção do delta médio.

# Sendo desejável controlar a independência entre essas variáveis para melhor compreender o tamanho do efeito de cada uma - a ser contemplado pelo modelo de regressão -, estabeleci um limite para elas atrelado à rr_dist.
# Assim, foram mantidas apenas as observações relativas a árvores alternativas que distassem, de sua árvore verdadeira, no máximo o mesmo que distavam entre si.

h0lev <- h0[which(mean_rt_KFdist<=rr_KFdist),]
detach(h0)
attach(h0lev)
# Essa redução da amostra (de 2989720 para 1539734 deltas) veio a introduzir um pressuposto adicional para o teste com o modelo regredido, mas essencial para garantir, sob a atual proposta, seu poder preditivo da distribuição de delta.
plot(rr_KFdist,mean_rt_KFdist,cex=0.1,xlab="BSD alternativa-alternativa",ylab="BSD alternativas-verdadeira",main="b",cex.lab=1.5,cex.axis=1.5,cex.main=2)



#===== OBTENÇÃO DAS DISTRIBUIÇÕES DE DELTA =====#

for (g in seq(0,0.4,0.1)){
	for (s in unique(seq_length)){
		for (t in seq(30,180,50)){
			for (d in  seq(0.0325,0.30,0.05)){
				cat(sprintf("gp%.2f_%dbp_%dtips_dist%.2f: %d\n",g,s,t,d,length(h0lev[which(gap_prop>g-0.02&gap_prop<g+0.02&seq_length==s&tips_number==t&rr_KFdist>=d-0.015&rr_KFdist<d+0.015),1])))
			}
		}
	}
}

###nivela diretamente para o total desejado, sem perda de obs nas frações menores
row <- 0
v <- matrix(ncol=29,nrow=1080)
s <- array(dim=35)
for (g in seq(0,0.4,0.1)){
	for (l in unique(seq_length)){
		for (t in seq(30,180,50)){
			for (dset in seq(0.0325,0.30,0.05)){
				set <- array(dim=35)
				temp <- h0lev[which(gap_prop>g-0.02&gap_prop<g+0.02&seq_length==l&tips_number==t&rr_KFdist>=dset-0.015&rr_KFdist<dset+0.015),2:36]
				detach(h0lev)
				attach(temp)
				lngs<-array()
				fracs <- list()
				n<-0
				for (d in seq(dset-0.014,dset+0.014,0.002)){
					n<-n+1
					fracs[[n]] <- as.matrix(temp[which(rr_KFdist>=d-0.001&rr_KFdist<d+0.001),])
					lngs[n]  <- length(fracs[[n]][,1])
				}
				clngs<-0
				o<-as.numeric(sprintf("%.f",200/n))
				for (i in sort(lngs,index.return=T)$ix){
					n<-n-1
					if(lngs[i]>o){
						set <- rbind(set,fracs[[i]][sample(c(1:lngs[i]),o),])
						clngs<-clngs+o
					}
					else{
						set <- rbind(set,fracs[[i]])
						clngs<-clngs+lngs[i]
					}
					o<-as.numeric(sprintf("%.f",(200-clngs)/n))
				}
				set <- set[-1,]
				class(set) <- "numeric"
				row <- row+1
				s <- rbind(s,set)
				v[row,1] <- median(set[,6])
				v[row,2] <- l
				v[row,3] <- t
				v[row,4] <- median(set[,10])
				v[row,5] <- median(set[,11])
				v[row,6] <- mean(set[,15])
				v[row,7] <- sd(set[,15])
				v[row,8] <- var(set[,15])
				v[row,9] <- mean(set[,4])
				v[row,10] <- median(set[,1])
				v[row,11] <- median(set[,17])
				v[row,12] <- mean(set[,18])
				v[row,13] <- mean(set[,19])
				v[row,14] <- mean(set[,20])
				v[row,15] <- mean(set[,21])
				v[row,16] <- mean(set[,22])
				v[row,17] <- mean(set[,23])
				v[row,18] <- mean(set[,24])
				v[row,19] <- mean(set[,25])
				v[row,20] <- mean(set[,26])
				v[row,21] <- mean(set[,27])
				v[row,22] <- mean(set[,28])
				v[row,23] <- mean(set[,29])
				v[row,24] <- mean(set[,30])
				v[row,25] <- mean(set[,31])
				v[row,26] <- mean(set[,32])
				v[row,27] <- mean(set[,33])
				v[row,28] <- mean(set[,34])
				v[row,29] <- mean(set[,35])
				detach(temp)
				attach(h0lev)
			}
		}
		cat(sprintf("%.0f%%\n",row/length(v[,1])*100))
	}
}
s <- s[-1,]
class(s) <- "numeric"
s <- as.data.frame(s)
colnames(v) <- c("gap_prop","seq_length","tips_number","rr_RFdist","rr_KFdist","delta_mean","delta_sd","delta_var","gamma_shape","total_k","tree_k","sum_b","sumlog_b","mean_b","sd_b","skew_b","kurt_b","min_b","quant10_b","quant20_b","quant30_b","quant40_b","median_b","quant60_b","quant70_b","quant80_b","quant90_b","max_b","mean_rt_KFdist")
class(v) <- "numeric"
v <- as.data.frame(v)
detach(h0lev)
attach(v)


#===== DA RELAÇÃO ENTRE DESVIO PADRÃO DE DELTA E AS VARIÁVEIS =====#

par(mfrow=c(2,2))
drange <- seq(0.0325,0.30,0.05)
rrd_mat <- matrix(ncol=6,nrow=length(rr_KFdist)/6)
c <- 0
for (d in drange){
	c <- c+1
	rrd_mat[,c] <- delta_sd[which(rr_KFdist>=d-0.015&rr_KFdist<d+0.015)]
}
plot(rr_KFdist,delta_sd,cex=0,xlim=c(0,0.3),xlab="BSD alternativa-alternativa",ylab="SD de ∆lnL",main="a",cex.lab=1.4,cex.axis=1.35,cex.main=1.7)
boxplot(rrd_mat,add=T,at=drange,cex=0.5,boxwex=0.015,xaxt='n',yaxt='n')

lrange <- sort(unique(seq_length))
sl_mat <- matrix(ncol=9,nrow=length(seq_length)/9)
c <- 0
for (l in lrange){
	c <- c+1
	sl_mat[,c] <- delta_sd[which(seq_length==l)]
}
plot(seq_length,delta_sd,cex=0,xlim=c(0,2000),xlab="Número de sítios",ylab="SD de ∆lnL",main="b",cex.lab=1.4,cex.axis=1.35,cex.main=1.7)
boxplot(sl_mat,add=T,at=lrange,cex=0.5,boxwex=50,xaxt='n',yaxt='n')

grange <- seq(0,0.4,0.1)
gp_mat <- matrix(ncol=5,nrow=length(gap_prop)/5)
c <- 0
for (g in grange){
	c <- c+1
	gp_mat[,c] <- delta_sd[which(gap_prop>g-0.02&gap_prop<g+0.02)]
}
plot(gap_prop,delta_sd,cex=0,xlim=c(-0.05,0.45),xlab="Proporção de gaps (%)",ylab="SD de ∆lnL",main="c",cex.lab=1.4,cex.axis=1.35,cex.main=1.7,xaxt='n')
axis(1,at=grange,labels=grange*100,cex.axis=1.35)
boxplot(gp_mat,add=T,at=grange,cex=0.5,boxwex=0.025,xaxt='n',yaxt='n')

trange <- seq(57,357,100)
tk_mat <- matrix(ncol=4,nrow=length(tree_k)/4)
c <- 0
for (t in trange){
	c <- c+1
	tk_mat[,c] <- delta_sd[which(tree_k==t)]
}
plot(tree_k,delta_sd,cex=0,xlim=c(30,400),xlab="Número de ramos",ylab="SD de ∆lnL",main="d",cex.lab=1.4,cex.axis=1.35,cex.main=1.7)
boxplot(tk_mat,add=T,at=trange,cex=0.5,boxwex=20,xaxt='n',yaxt='n')

#após logaritimizações...

par(mfrow=c(2,2))
drange <- seq(0.0325,0.30,0.05)
rrd_mat <- matrix(ncol=6,nrow=length(rr_KFdist)/6)
c <- 0
for (d in drange){
	c <- c+1
	rrd_mat[,c] <- log(delta_sd)[which(rr_KFdist>=d-0.015&rr_KFdist<d+0.015)]
}
drange <- log(seq(0.0325,0.30,0.05))
plot(log(rr_KFdist),log(delta_sd),xlim=c(-3.5,-1.2),cex=0,xlab="Log do BSD alternativa-alternativa",ylab="Log do SD de ∆lnL",main="a",cex.lab=1.4,cex.axis=1.35,cex.main=1.7)
boxplot(rrd_mat,add=T,at=drange,cex=0.5,boxwex=0.15,xaxt='n',yaxt='n')

lrange <- sort(unique(log(seq_length)))
sl_mat <- matrix(ncol=9,nrow=length(seq_length)/9)
c <- 0
for (l in lrange){
	c <- c+1
	sl_mat[,c] <- log(delta_sd)[which(log(seq_length)==l)]
}
plot(log(seq_length),log(delta_sd),cex=0,xlab="Log do número de sítios",ylab="Log do SD de ∆lnL",main="b",cex.lab=1.4,cex.axis=1.35,cex.main=1.7)
boxplot(sl_mat,add=T,at=lrange,cex=0.5,boxwex=0.15,xaxt='n',yaxt='n')

grange <- seq(0,0.4,0.1)
gp_mat <- matrix(ncol=5,nrow=length(gap_prop)/5)
c <- 0
for (g in grange){
	c <- c+1
	gp_mat[,c] <- log(delta_sd)[which(gap_prop>g-0.02&gap_prop<g+0.02)]
}
plot(gap_prop,log(delta_sd),cex=0,xlim=c(-0.05,0.45),xlab="Proporção de gaps (%)",ylab="Log do SD de ∆lnL",main="c",cex.lab=1.4,cex.axis=1.35,cex.main=1.7,xaxt='n')
axis(1,at=grange,labels=grange*100,cex.axis=1.35)
boxplot(gp_mat,add=T,at=grange,cex=0.5,boxwex=0.025,xaxt='n',yaxt='n')

trange <- log(seq(57,357,100))
tk_mat <- matrix(ncol=4,nrow=length(tree_k)/4)
c <- 0
for (t in trange){
	c <- c+1
	tk_mat[,c] <- log(delta_sd)[which(log(tree_k)==t)]
}
plot(log(tree_k),log(delta_sd),cex=0,xlim=c(3.9,6),xlab="Log do número de ramos",ylab="Log do SD de ∆lnL",main="d",cex.lab=1.4,cex.axis=1.35,cex.main=1.7)
boxplot(tk_mat,add=T,at=trange,cex=0.5,boxwex=0.125,xaxt='n',yaxt='n')





#===== DO MODELO DE REGRESSÃO =====#

#modelo equidistant delta
ed <- lm(log(delta_sd)~log(rr_KFdist)+gap_prop+log(seq_length)+log(tree_k))
summary(ed)
#Call:
#lm(formula = log(delta_sd) ~ log(rr_KFdist) + log(tree_k) + log(seq_length) + 
#    gap_prop)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.49542 -0.10662 -0.01041  0.10230  0.47577 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     -1.292388   0.053193  -24.30   <2e-16 ***
#log(rr_KFdist)   0.944579   0.006942  136.07   <2e-16 ***
#log(tree_k)      0.294764   0.007235   40.74   <2e-16 ***
#log(seq_length) -0.201262   0.005492  -36.64   <2e-16 ***
#gap_prop        -1.008815   0.035497  -28.42   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.1648 on 1075 degrees of freedom
#Multiple R-squared:  0.9539,	Adjusted R-squared:  0.9537 
#F-statistic:  5563 on 4 and 1075 DF,  p-value: < 2.2e-16

#lnSD = 0.945 * lnBSDt1t2 + 0.295 * lnKt1t2 - 0.201 * lnSm - 1.009 * GPm -1.292

#plot homocedasticidade dos resíduos
plot(ed$fitted,ed$residuals,cex=0.2,xlab="lnSDs previstos",ylab=expression(italic("E")),main="a",cex.lab=1.4,cex.axis=1.35,cex.main=2.2)
segments(-5.2,0,-1.5,0,lty=2,col="red3",lwd=1.5)

#histograma e curva de probabilidade normal
hist(ed$residuals,nclass=20,prob=T,xlab=expression(italic("E")),ylab="Probabilidade",main="b",cex.lab=1.4,cex.axis=1.35,cex.main=2.2)
curve(dnorm(x,mean(ed$residuals),sd(ed$residuals)),add=T,col="red3",lwd=1.5,lty=2)

#plot quantil-quantil normal vs. residuos
qqnorm(ed$residuals,xlab="Quantis Normal Padrão",ylab=substitute(paste("Quantis ",italic("E"))),main="c",cex=0.2,cex.lab=1.4,cex.axis=1.35,cex.main=2.2)
qqline(ed$residuals,col="red3",lwd=1.5,lty=2)
