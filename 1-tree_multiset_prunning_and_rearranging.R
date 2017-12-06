library(phangorn)

BSrightpyramid <- function(toptree, max.tries=300, D.range=c(0.01,0.02), max.top.DDiff=0.01){
	moves.array <- seq(D.range[1]*100,D.range[2]*100,1)
 	for(moves in moves.array){
 		tries <- 0
 		basetrees <- list()
	  	repeat{
	  		b1 <- rNNI(toptree, moves=moves)
	  		if(is.binary.tree(b1)){
	  			tb1.RFdistance <- RF.dist(toptree, b1)
	  			tb1.KFdistance <- KF.dist(toptree, b1)
	  			basetrees[[1]] <- b1
	  			break
	  		}
	  	}
	  	repeat{
			b2 <- rNNI(toptree, moves=moves)
			if(is.binary.tree(b2)){
				tries <- tries+1
				b1b2.KFdistance <- KF.dist(b1,b2)
				tb2.RFdistance <- RF.dist(toptree,b2)
				tb2.KFdistance <- KF.dist(toptree,b2)
				if(abs(tb1.KFdistance-tb2.KFdistance)<max.top.DDiff&tb1.RFdistance==tb2.RFdistance&b1b2.KFdistance>D.range[1]&b1b2.KFdistance<=D.range[2]){
					basetrees[[2]] <- b2
					break
				}
				else if(tries>=max.tries){
					break
				}
			}
	 	}
	 	if (length(basetrees)<2){
	 		next
	 	}
	 	repeat{
	 		b3 <- rNNI(toptree, moves=moves)
	 		if(is.binary.tree(b3)){
				tries <- tries+1
				tb3.RFdistance <- RF.dist(toptree,b3)
				tb3.KFdistance <- KF.dist(toptree,b3)
				b1b3.KFdistance <- KF.dist(b1,b3)
				b2b3.KFdistance <- KF.dist(b2,b3)
				if(abs(tb1.KFdistance-tb3.KFdistance)<max.top.DDiff&abs(tb2.KFdistance-tb3.KFdistance)<max.top.DDiff&tb1.RFdistance==tb3.RFdistance&tb2.RFdistance==tb3.RFdistance&b1b3.KFdistance>D.range[1]&b1b3.KFdistance<=D.range[2]&b2b3.KFdistance>D.range[1]&b2b3.KFdistance<=D.range[2]){
					basetrees[[3]] <- b3
					basetrees[[4]] <- moves
					basetrees[[5]] <- matrix(nrow=3,ncol=3)
					for(i in c(1,2,3)){
						for(j in c(1,2,3)){
							if(i!=j){
								basetrees[[5]][i,j] <- KF.dist(basetrees[[i]],basetrees[[j]])
							}
							else{
								basetrees[[5]][i,j] <- 0
								break
							}
						}
					}
					break
				}
				else if(tries>=max.tries*2){
					break
				}
	 		}
	 	}
	 	if (length(basetrees)<3){
	 		next
	 	}
	 	break
	}
	return(basetrees)
}



#reag.trees <- list()
#openrun <- 1
models <- c("GTR")
for (mdl in 1:length(models)){
	fulltree <- read.tree(sprintf("FoleySpringerTeelingPRSB2016_%s+FO+G8_0.03cut.treefile",models[mdl]))
	D.seq <- seq(0.05,0.35,0.05)
	drop.seq <- seq(150,0,-50)
	reag.trees[[mdl]] <- list()
	for (tx in 3:length(drop.seq)){
		if(tx==length(drop.seq)){
			system(sprintf("cp FoleySpringerTeelingPRSB2016_%s+FO+G8_%dtaxa.treefile true.tree_%s_%dtips.tre",models[mdl],180-drop.seq[tx],models[mdl],180-drop.seq[tx]))
			openrun <- 1
		}
	  reag.trees[[mdl]][[tx]] <- list()
	  repeat{
	  	  if(openrun>0 & !is(try(read.tree(sprintf("FoleySpringerTeelingPRSB2016_%s+FO+G8_%dtaxa.ckp.gz",models[mdl],180-drop.seq[tx])),silent=TRUE),"try-error")){
	  	  	system(sprintf("iqtree-omp -s FoleySpringerTeelingPRSB2016_Alignment_%dtaxa.phy -pre FoleySpringerTeelingPRSB2016_%s+FO+G8_%dtaxa -m %s+FO+G8 -nt 2",180-drop.seq[tx],models[mdl],180-drop.seq[tx],models[mdl]))
	  	  	truet <- read.tree(sprintf("FoleySpringerTeelingPRSB2016_%s+FO+G8_%dtaxa.treefile",models[mdl],180-drop.seq[tx]))
	  	  	write.tree(truet,sprintf("true.tree_%s_%dtips.tre",models[mdl],180-drop.seq[tx]))
	  	  	cat(sprintf("defining true tree for reag.trees[[%s]][[%d]]\n",mdl,tx))
	  	  }
	  	  else{
	  	  	chopt <- drop.tip(fulltree,fulltree$tip.label[sample(c(1:length(fulltree$tip.label)),drop.seq[tx])],rooted=FALSE)
	  	  	write.tree(chopt,sprintf("choped.tree_%s_%dtips.tre",models[mdl],180-drop.seq[tx]))
	  	  	openrun <- 1
	  	  	system(sprintf("perl alignSeqSelector.pl choped.tree_GTR_%dtips.tre FoleySpringerTeelingPRSB2016_Alignment.phy > FoleySpringerTeelingPRSB2016_Alignment_%dtaxa.phy",180-drop.seq[tx],180-drop.seq[tx]))
	  	  	system(sprintf("iqtree-omp -s FoleySpringerTeelingPRSB2016_Alignment_%dtaxa.phy -pre FoleySpringerTeelingPRSB2016_%s+FO+G8_%dtaxa -m %s+FO+G8 -nt 2 -redo",180-drop.seq[tx],models[mdl],180-drop.seq[tx],models[mdl]))
	  	  	truet <- read.tree(sprintf("FoleySpringerTeelingPRSB2016_%s+FO+G8_%dtaxa.treefile",models[mdl],180-drop.seq[tx]))
	  	  	write.tree(truet,sprintf("true.tree_%s_%dtips.tre",models[mdl],180-drop.seq[tx]))
	  	  	cat(sprintf("defining true tree for reag.trees[[%s]][[%d]]\n",mdl,tx))
	  	  }
		  for (set in 1:3){
		    reag.trees[[mdl]][[tx]][[set]] <- list()
		    for (d in 1:(length(D.seq))){
		    	tries <- 1
		      	repeat{
			        reag.trees[[mdl]][[tx]][[set]][[d]] <- BSrightpyramid(truet,D.range=c(D.seq[[d]]-0.03,D.seq[[d]]+0.02))
			        if (length(reag.trees[[mdl]][[tx]][[set]][[d]])==5){
			        	if(set==2){
			        		for(i in 1:3){
			        			for(j in 1:3){
			        				if(RF.dist(reag.trees[[mdl]][[tx]][[set]][[d]][[i]],reag.trees[[mdl]][[tx]][[set-1]][[d]][[j]])==0){
			        					break
			        				}
			        			}
		        				if(RF.dist(reag.trees[[mdl]][[tx]][[set]][[d]][[i]],reag.trees[[mdl]][[tx]][[set-1]][[d]][[j]])==0){
			        				break
			        			}
			        		}
	        				if(RF.dist(reag.trees[[mdl]][[tx]][[set]][[d]][[i]],reag.trees[[mdl]][[tx]][[set-1]][[d]][[j]])>0){
								cat(sprintf("reag.trees[[%s]][[%d]][[%d]][[%d]] done\n",mdl,tx,set,d))
								for(i in 1:3){
									write.tree(reag.trees[[mdl]][[tx]][[set]][[d]][[i]],sprintf("reag.tree_%s_%dtips_dist%.2f_set%d_t%d.tre",models[mdl],180-drop.seq[tx],D.seq[d],set,i))
								}
								break
		        			}
		        			else{
		        				cat(sprintf("retrying reag.trees[[%s]][[%d]][[%d]][[%d]]\n",mdl,tx,set,d))
		        			}
			        	}
			        	else if(set==3){
			        		for(i in 1:3){
			        			for(j in 1:3){
			       					if(RF.dist(reag.trees[[mdl]][[tx]][[set]][[d]][[i]],reag.trees[[mdl]][[tx]][[set-1]][[d]][[j]])==0|RF.dist(reag.trees[[mdl]][[tx]][[set]][[d]][[i]],reag.trees[[mdl]][[tx]][[set-2]][[d]][[j]])==0){
			       						break
			       					}
			        			}
		        				if(RF.dist(reag.trees[[mdl]][[tx]][[set]][[d]][[i]],reag.trees[[mdl]][[tx]][[set-1]][[d]][[j]])==0|RF.dist(reag.trees[[mdl]][[tx]][[set]][[d]][[i]],reag.trees[[mdl]][[tx]][[set-2]][[d]][[j]])==0){
			        				break
			        			}
			        		}
	        				if(RF.dist(reag.trees[[mdl]][[tx]][[set]][[d]][[i]],reag.trees[[mdl]][[tx]][[set-1]][[d]][[j]])>0&RF.dist(reag.trees[[mdl]][[tx]][[set]][[d]][[i]],reag.trees[[mdl]][[tx]][[set-2]][[d]][[j]])>0){
								cat(sprintf("reag.trees[[%s]][[%d]][[%d]][[%d]] done\n",mdl,tx,set,d))
								for(i in 1:3){
									write.tree(reag.trees[[mdl]][[tx]][[set]][[d]][[i]],sprintf("reag.tree_%s_%dtips_dist%.2f_set%d_t%d.tre",models[mdl],180-drop.seq[tx],D.seq[d],set,i))
								}
								break
		        			}
		        			else{
		        				cat(sprintf("retrying reag.trees[[%s]][[%d]][[%d]][[%d]]\n",mdl,tx,set,d))
		        			}
			        	}
			        	else{
							cat(sprintf("reag.trees[[%s]][[%d]][[%d]][[%d]] done\n",mdl,tx,set,d))
							for(i in 1:3){
								write.tree(reag.trees[[mdl]][[tx]][[set]][[d]][[i]],sprintf("reag.tree_%s_%dtips_dist%.2f_set%d_t%d.tre",models[mdl],180-drop.seq[tx],D.seq[d],set,i))
							}
							break
			        	}
			        }
			        else{
			        	tries <- tries+1
			        	cat(sprintf("retrying reag.trees[[%s]][[%d]][[%d]][[%d]]\n",mdl,tx,set,d))
			        	if(tries>20*set && tx<length(drop.seq) && set<3){
			        		cat(sprintf("retrying reag.trees[[%s]][[%d]] with new subtree\n",mdl,tx))
			        		break
			        	}
			        }
		      	}
			    if(length(reag.trees[[mdl]][[tx]][[set]][[d]])<5){
			    	break
			    }
		    }
		    if(length(reag.trees[[mdl]][[tx]][[set]][[d]])<5){
		    	openrun <- 0
		    	break
		    }
		   	save.image("RearrangedTrees.RData")
			cat("saving RData\n")
		}
		if(length(reag.trees[[mdl]][[tx]][[set]][[d]])==5){
			openrun <- 0
			break
		}
	  }
	}
}


