netmatch<-function(net,name){
  pipei1<-match(net[,1],name)
  pipei2<-match(net[,2],name)
  newnet<-cbind(pipei1,pipei2,net[,3],deparse.level = 0 )
  index1 <- which(is.na(pipei1))
  index2<- which(is.na(pipei2))
  newnet<-newnet[-union(index1,index2),]
  return(newnet)
}
#initdata
initdata<-function(row,L){
  low<- matrix(runif(row*L,min=0,max=1),row,L)
  return(low)
}

#sigmoid


sigmoid = function(x,y) {
  1 / (1 + exp(-x%*%y))
}



#negtable
neg<-function(d){
  m<-c()
  pc<-rowSums(d)^0.75
  pc<-floor((pc/sum(pc)*1e6))
  index<-which(pc!=0)
  for(i in 1:length(index))
  {
    m<-c(m,(rep(1,pc[index[i]])*index[i]))
  }
  return(m)
}

#edge
edg<-function(d){
  dh<-as.vector(d)
  index<-which(dh!=0)
  dh2<-dh[c(index)]
  dh2m<-dh2/sum(as.numeric(dh2))
  dht<-floor(dh2m*1e7)
  dhcum<-cumsum(dht)
  ze=rep(0.1,dhcum[c(length(dhcum))])
  for(i in 2:length(dhcum)){
    ze[dhcum[i-1]:dhcum[i]]<-index[i]
  }
  ze[1:dhcum[1]]=index[1]
  return(ze)
}

#compedg
com<-function(index,v,T){
  tmp<-index[ceiling(runif(T,max=length(index)))]
  row<-ceiling(tmp/v)
  col<-tmp%%v
  posi<-list(row,col)
  return(posi)
}

matga<-function(path){
  ga<-path[,2]
  b1<-table(ga)
  b2<-b1^0.75
  b3<-floor(b2/sum(b2)*1e6)
  na<-names(b3)
  mat_ga<-c()
  for(i in 1:length(b3))
  {
    mat_ga<-c(mat_ga,(rep(1,b3[i])*as.numeric(na[i])))
  }
  return(mat_ga)
}


edge<-function(path){
  sp<-path[,3]
  spm<-floor(sp/sum(as.numeric(sp))*1e6)
  spmcum<-cumsum(spm)
  edgesum<-rep(0.1,spmcum[length(spmcum)])
  edgesum[1:spmcum[1]]<-1
  for(i in 2:length(spm)){
    edgesum[spmcum[i-1]:spmcum[i]]<-i
  }
  return(edgesum)
}

sam<-function(index,T,K){
  l<-index[ceiling(runif(T*K,max=length(index)))]
  l<-matrix(l,nrow = T,ncol = K)
}

compath<-function(index,path,T){
  suiji<-index[ceiling(runif(T,max=length(index)))]
  row<-as.numeric(path[suiji,1])
  col<-as.numeric(path[suiji,2])
  posi<-list(row,col)
  return(posi)
}

#main
#' Integrated networks embedding for scRNA-seq data
#' @usage scLINE(mat, network, L, K = 5, T, rho = 0.025, Random_weight)
#' @param mat A matrix of scRNA-seq data in which rows represent gene columns representing cells
#' @param network A list of gene networks with 3 columns: geneID, geneID, score
#' @param L Dim of the low-dimensional representations
#' @param K Number of negative samples, the default is 5
#' @param T Number of iterations
#' @param rho The initial learning rate, the default is 0.025
#' @param Random_weight True of False,Whether to add random weights to the network model
#' @return the low-dimensional matrix of input data
#' @export
#' @description A tool of dimension reduction for scRNA-seq data combined with multiple gene networks embedding model
#' @examples
#' load(system.file("data","Usoskin.Rdata",package = "scLINE"))
#'load(system.file("data","ppi.Rdata",package = "scLINE"))
#'load(system.file("data","humannet.Rdata",package = "scLINE"))
#'exp_mat<-Usoskin$rawdata
#'gene_network<-list(ppi = ppi,humannet = humannet)
#'lowdim_list<-scLINE(exp_mat, gene_network, L = 20, K = 5, T = 1e7, rho = 0.025, Random_weight = TRUE)


scLINE<-function(mat,network,L,K=5,T,rho=0.025,Random_weight){
  if(!is.matrix(mat)){
    stop( "'mat' must be a expression count matrix")
  }
  if(!is.list(network)){
    stop( "'network' must be a list")
  }
  N<-length(network)
  genename<-rownames(mat)
  message("Mathching geneID...")
  for (i in 1:N) {
    network[[i]]<-netmatch(network[[i]],genename)
  }
  for (i in 1:N) {
    if(length(network[[i]])==0){
      network<-network[-i]
      message("remove the network ",i,", no matched gene")
    }
  }
  N<-length(network)
  n<-ncol(mat)
  m<-nrow(mat)
  rowex<-ncol(mat)
  rowge<-nrow(mat)
  cell_low<-initdata(rowex,L)
  congene_low<-initdata(rowge,L)
  gene_low<-initdata(rowge,L)
  cell_low<-cell_low*2-1
  gene_low<-gene_low*2-1
  congene_low<-congene_low*2-1
  for(i in 1:n)
  {
    a<-cell_low[i,]
    cell_low[i,]<-cell_low[i,]/norm(a,type = "2")
  }
  for(i in 1:m)
  {
    b<-(congene_low[i,])
    congene_low[i,]<-congene_low[i,]/norm(b,type = "2")
  }
  for(i in 1:m)
  {
    c<-(gene_low[i,])
    gene_low[i,]<-gene_low[i,]/norm(b,type = "2")
  }
  message("Constructing mapping tables for sampling...")
  mat_ca<-neg(mat)
  mat_E<-edg(mat)
  edgesump <- c(paste0('edgesump',1:N))
  for (i in 1:N) {
    assign(edgesump[i],edge(network[[i]]))
  }

  r<-rho
  lm<-sam(mat_ca,T,K)

  comg<-com(mat_E,m,T)
  ww1<-comg[[1]]
  gg1<-comg[[2]]
  comp<-c(paste0('comp',1:N))
  ww<-c(paste0('ww',2:(N+1)))
  gg<-c(paste0('gg',2:(N+1)))
  for (i in 1:N) {
    edgeget<-get(paste("edgesump",i,sep=""))
    assign(comp[i],compath(edgeget,network[[i]],T))
    compget<-get(paste("comp",i,sep=""))
    assign(ww[i],compget[[1]])
    assign(gg[i],compget[[2]])
  }

  for(k in 1:T){
    ggs<-c(gg1[k])
    for (j in 2:(N+1)) {
      ggget<-get(paste("gg",j,sep=""))
      ggs<-c(ggs,ggget[k])
    }
    while (length(intersect(lm[k,],ggs))){
      lm[k,]<-sample(mat_ca,K)
    }
  }
  rho=0.025
  r<-rho
  message("Start iterating...")
  timestart<-Sys.time()
  if (Random_weight==TRUE){
    for(k in 1:T){
      weight<-runif(N,min=0,max=1)
      if((k>=(T/10))&(k%%(T/10)==0))
      {
        timeend<-Sys.time()
        tim<-timeend-timestart
        cat(k/T*100,"% completed,")
        print(tim)
      }

      #update c-g
      p1 = sigmoid((-1)*congene_low[gg1[k],],cell_low[ww1[k],])
      d_h<-c(p1)*cell_low[ww1[k],]
      d_u<-c(p1)*congene_low[gg1[k],]
      d_u<-d_u-colSums(sweep(congene_low[lm[k,],],1,sigmoid(congene_low[lm[k,],],cell_low[ww1[k],]),"*"))
      congene_low[lm[k,],]<-congene_low[lm[k,],]-matrix(outer(r*sigmoid(congene_low[lm[k,],],cell_low[ww1[k],]),cell_low[ww1[k],]),nrow = K)
      cell_low[ww1[k],]=cell_low[ww1[k],]+r*d_u
      congene_low[gg1[k],]=congene_low[gg1[k],]+r*d_h
      #update path1
      for (i in 2:N+1) {
        wwget<-get(paste("ww",i,sep=""))
        ggget<-get(paste("gg",i,sep=""))
        p<-weight[i-1]*sigmoid((-1)*congene_low[ggget[k],],gene_low[wwget[k],])
        d_h<-c(p)*gene_low[wwget[k],]
        d_v<-c(p)*congene_low[ggget[k],]
        d_v<-d_v-colSums(sweep(congene_low[lm[k,],],1,sigmoid(congene_low[lm[k,],],gene_low[wwget[k],]),"*"))
        congene_low[lm[k,],]<-congene_low[lm[k,],]-matrix(outer(r*sigmoid(congene_low[lm[k,],],gene_low[wwget[k],]),gene_low[wwget[k],]),nrow = K)

        gene_low[wwget[k],]=gene_low[wwget[k],]+r*d_v

        congene_low[ggget[k],]=congene_low[ggget[k],]+r*d_h
      }
      r<-rho*(1-k/T)
    }
  }
  if (Random_weight==FALSE){
    for(k in 1:T){
      if((k>=1e5)&(k%%1e5==0))
      {
        timeend<-Sys.time()
        tim<-timeend-timestart
        print(tim)
        timestart<-Sys.time()
        print(k)
      }


      #update c-g
      p1 <- sigmoid((-1)*congene_low[gg1[k],],cell_low[ww1[k],])
      d_h<-c(p1)*cell_low[ww1[k],]
      d_u<-c(p1)*congene_low[gg1[k],]
      d_u<-d_u-colSums(sweep(congene_low[lm[k,],],1,sigmoid(congene_low[lm[k,],],cell_low[ww1[k],]),"*"))
      congene_low[lm[k,],]<-congene_low[lm[k,],]-matrix(outer(r*sigmoid(congene_low[lm[k,],],cell_low[ww1[k],]),cell_low[ww1[k],]),nrow = K)
      cell_low[ww1[k],]=cell_low[ww1[k],]+r*d_u
      congene_low[gg1[k],]=congene_low[gg1[k],]+r*d_h
      #update path1
      for (i in 2:N+1) {
        wwget<-get(paste("ww",2,sep=""))
        ggget<-get(paste("gg",2,sep=""))
        p<-sigmoid((-1)*congene_low[ggget[k],],gene_low[wwget[k],])
        d_h<-c(p)*gene_low[wwget[k],]
        d_v<-c(p)*congene_low[ggget[k],]
        d_v<-d_v-colSums(sweep(congene_low[lm[k,],],1,sigmoid(congene_low[lm[k,],],gene_low[wwget[k],]),"*"))
        congene_low[lm[k,],]<-congene_low[lm[k,],]-matrix(outer(r*sigmoid(congene_low[lm[k,],],gene_low[wwget[k],]),gene_low[wwget[k],]),nrow = K)

        gene_low[wwget[k],]=gene_low[wwget[k],]+r*d_v

        congene_low[ggget[k],]=congene_low[ggget[k],]+r*d_h
      }

      r<-rho*(1-k/T)
    }
  }
  low_mat<-list(cell_low=cell_low,gene_low=gene_low)
  return(low_mat)
}
