# imports from multcomp::glht
# imports from PMCMR:posthoc.kruskal.nemenyi.test

easyDes=function(nc.g=NULL,nc.n=NULL,nc.f=NULL,dataIn=data,fisher=TRUE,aov=FALSE,t=FALSE){

  ###################################################
  ## function of numeric variables from two groups ##
  ###################################################

  t.wilcox.des <- function(nc.g=nc.g,nc.n=nc.n,dataInS=dataIn,t2=t){

    if(!t2){

    t.wilcox.single <- function(g=nc.g,i,data=dataInS){
      data[,i]=as.numeric(data[,i])
      m0=mean(data[,i],na.rm=TRUE)
      s0=sd(data[,i],na.rm=TRUE)
      ms0=paste(round(m0,3),"+/-",round(s0,3),sep="")
      m=aggregate(data[,i],by=list(data[,g]),mean,na.rm=TRUE)
      s=aggregate(data[,i],by=list(data[,g]),sd,na.rm=TRUE)
      rst.single=paste(round(m[,2],3),"+/-",round(s[,2],3),sep="")
      shapiro1=shapiro.test(data[,i][data[,g]==levels(data[,g])[1]])
      shapiro2=shapiro.test(data[,i][data[,g]==levels(data[,g])[2]])
      shapiro=min(shapiro1$p.value,shapiro2$p.value)
      bartlett=bartlett.test(data[,i]~data[,g])
      if(shapiro<0.1 | bartlett$p.value<0.1){
        temp=wilcox.test(data[,i]~data[,g],exact=FALSE)
        s=c("Wilcoxon test",round(temp$statistic,3))
        p=round(temp$p.value,3)}
      if(shapiro>=0.1 & bartlett$p.value<0.1){
        temp=t.test(data[,i]~data[,g],var.equal=FALSE)
        s=c("Adjust t test",round(temp$statistic,3))
        p=round(temp$p.value,3)}
      if(shapiro>=0.1 & bartlett$p.value>=0.1){
        temp=t.test(data[,i]~data[,g],var.equal=TRUE)
        s=c("t test",round(temp$statistic,3))
        p=round(temp$p.value,3)}
      rst.single=c(ms0,rst.single,s,p)
      rst.single=data.frame(matrix(rst.single,nrow=1))
      names(rst.single)=c("total",levels(data[,nc.g]),"method","statistic","p.value")
      row.names(rst.single)=names(data)[i]
      return(rst.single)
    }

    }#if(!t2)

    if(t2){

      t.wilcox.single <- function(g=nc.g,i,data=dataInS){
        data[,i]=as.numeric(data[,i])
        m0=mean(data[,i],na.rm=TRUE)
        s0=sd(data[,i],na.rm=TRUE)
        ms0=paste(round(m0,3),"+/-",round(s0,3),sep="")
        m=aggregate(data[,i],by=list(data[,g]),mean,na.rm=TRUE)
        s=aggregate(data[,i],by=list(data[,g]),sd,na.rm=TRUE)
        rst.single=paste(round(m[,2],3),"+/-",round(s[,2],3),sep="")
        temp=t.test(data[,i]~data[,g],var.equal=TRUE)
        s=c("t test",round(temp$statistic,3))
        p=round(temp$p.value,3)
        rst.single=c(ms0,rst.single,s,p)
        rst.single=data.frame(matrix(rst.single,nrow=1))
        names(rst.single)=c("total",levels(data[,nc.g]),"method","statistic","p.value")
        row.names(rst.single)=names(data)[i]
        return(rst.single)
      }

    }#if(t2)


    n.rst=t.wilcox.single(i=nc.n[1],data=dataInS)

    if(length(nc.n)==1){n.rst=n.rst}else{
      for(j in nc.n[2:length(nc.n)]){n.rst=rbind(n.rst,t.wilcox.single(i=j,data=dataInS))}
    }

    n.rst$p.value=as.numeric(as.character(n.rst$p.value))
    return(n.rst)
  }

  ##################################################
  ## function of factor variables from two groups ##
  ##################################################

  fisher.des=function(nc.g=nc.g,nc.f=nc.f,dataInS=dataIn,ftestS=fisher){

    fisher.des.single=function(g=nc.g,i,data=dataInS,test=ftestS){
      data[,i]=as.factor(data[,i])
      t=table(data[,i])
      p=prop.table(t)
      p=round(p*100,2)
      rst0=paste(t," (",p,")",sep="")
      t=table(data[,i],data[,g])
      p=prop.table(t,2)
      p=round(p*100,2)
      rst=data.frame(matrix(paste(t," (",p,")",sep=""),ncol=length(levels(data[,g]))))
      if(test){
        fisher=fisher.test(t)
        p=round(fisher$p.value,3)
        method=c("Fisher test",rep(NA,nrow(t)-1))
        statistic=rep(NA,nrow(t))
        p=c(p,rep(NA,nrow(t)-1))
        rst=cbind(rst0,rst,method,statistic,p)
        names(rst)=c("total",levels(data[,g]),"method","statistic","p.value")
        row.names(rst)=paste(names(data)[i],levels(data[,i]))}
      if(!test){
        chisq=chisq.test(t)
        p=round(chisq$p.value,3)
        method=c("Chi-square test",rep(NA,nrow(t)-1))
        statistic=rep(NA,nrow(t))
        p=c(p,rep(NA,nrow(t)-1))
        rst=cbind(rst0,rst,method,statistic,p)
        names.rst=paste("group",1:length(levels(data[,g])),sep="")
        names(rst)=c("total",names.rst,"method","statistic","p.value")
        row.names(rst)=paste(names(data)[i],levels(data[,i]))}

      return(rst)
    }

    f.rst=fisher.des.single(g=nc.g,i=nc.f[1],data=dataInS,test=ftestS)

    if(length(nc.f)==1){f.rst=f.rst}else{
      for(jj in 2:length(nc.f)){
        f.rst=rbind(f.rst,fisher.des.single(g=nc.g,i=nc.f[jj],data=dataInS,test=ftestS))
      }
    }

    return(f.rst)
  }

  ####################################################
  ## function of numeric variables from more groups ##
  ####################################################

  aov.kru=function(g=nc.g,j=nc.n,data=dataIn,aov=aov){

    names(data)[names(data)=="group"]="group_temp"

    if(aov){
      aov.single=function(data,i,g){
        data[,g]=as.factor(data[,g])
        data[,i]=as.numeric(data[,i])
        names(data)[g]="group"
        g.lev=levels(data[,g])
        p=NA
        model=formula(paste(names(data)[i],"~ group"))
        temp=aov(model,data)
        p[1]="ANOVA"
        p[2]=round(summary(temp)[[1]][1,4],3)
        p[3]=round(summary(temp)[[1]][1,5],3)
        temp=summary(glht(temp,linfct=mcp(group="Tukey")))#Dunnett
        temp=temp$test[[6]]
        p[4:(3+choose(length(g.lev),2))]=round(temp,3)
        p=data.frame(t(p))
        r.name=NA
        for(k in 1:(length(g.lev)-1)){r.name=c(r.name,paste(g.lev[k],"_vs_",g.lev[(k+1):length(g.lev)],sep=""))}
        rownames(p)=c(names(data)[i])
        colnames(p)=c("method","statistic","p.value",r.name[-1])

        m0=mean(data[,i],na.rm=T)
        s0=sd(data[,i],na.rm=T)
        ms0=paste(round(m0,3),"+/-",round(s0,3),sep="")
        m=aggregate(data[,i],by=list(data[,g]),mean,na.rm=TRUE)
        s=aggregate(data[,i],by=list(data[,g]),sd,na.rm=TRUE)
        ms=paste(as.matrix(round(m[,-1],3)),"+/-",as.matrix(round(s[,-1],3)),sep="")
        ms=matrix(ms,byrow=T,nrow=1)
        des=data.frame(cbind(ms0,ms))
        names(des)=c("total",levels(data[,g]))
        rownames(des)=names(data)[i]

        p=cbind(des,p)
        return(p)
      }

      if(length(j)>1){
        rst=aov.single(data,j[1],g)
        for(jj in 2:length(j)){
          rst=rbind(rst,aov.single(data,j[jj],g))
        }
      }
      if(length(j)==1){
        rst=aov.single(data,j[1],g)
      }

      return(rst)

    }#if(aov)

    if(!aov){

      aov.kru.single=function(data,i,g){
        data[,g]=as.factor(data[,g])
        data[,i]=as.numeric(data[,i])
        names(data)[g]="group"
        g.lev=levels(data[,g])
        p.shapiro=NA
        for(g.num in 1:length(g.lev)){
          shapiro=shapiro.test(data[,i][data[,g]==g.lev[g.num]])
          p.shapiro[g.num]=shapiro$p.value
        }
        bartlett=bartlett.test(data[,i]~data[,g])
        p=NA
        model=formula(paste(names(data)[i],"~ group"))
        if(min(p.shapiro)<0.1 | bartlett$p.value<0.1){
          temp=kruskal.test(model,data)
          p[1]="Kruskal-Wallis test"
          p[2]=round(temp$statistic,3)
          p[3]=round(temp$p.value,3)
          temp=posthoc.kruskal.nemenyi.test(x=data[,i],g=data[,g],method="Chisq")#Tukey
          p.post=as.numeric(temp$p.value)
          p.post=p.post[!is.na(p.post)]
          p[4:(3+choose(length(g.lev),2))]=round(p.post,3)
        }
        if(min(p.shapiro)>=0.1 &  bartlett$p.value>=0.1){
          temp=aov(model,data)
          p[1]="ANOVA"
          p[2]=round(summary(temp)[[1]][1,4],3)
          p[3]=round(summary(temp)[[1]][1,5],3)
          temp=summary(glht(temp,linfct=mcp(group="Tukey")))#Dunnett
          temp=temp$test[[6]]
          p[4:(3+choose(length(g.lev),2))]=round(temp,3)
        }
        p=data.frame(t(p))
        r.name=NA
        for(k in 1:(length(g.lev)-1)){r.name=c(r.name,paste(g.lev[k],"_vs_",g.lev[(k+1):length(g.lev)],sep=""))}
        rownames(p)=c(names(data)[i])
        colnames(p)=c("method","statistic","p.value",r.name[-1])

        if(p[1]=="Kruskal-Wallis test"){
          q0=quantile(data[,i],na.rm=T)
          qr0=paste(round(q0[3],3),"(",round(q0[2],3),",",round(q0[4],3),")",
                    sep="")
          q=aggregate(data[,i],by=list(data[,g]),quantile,na.rm=TRUE)
          qr=paste(round(q[,-1][,3],3),"(",round(q[,-1][,2],3),
                   ",",round(q[,-1][,4],3),")",sep="")
          qr=matrix(qr,byrow=T,nrow=1)
          des=data.frame(cbind(qr0,qr))
          names(des)=c("total",levels(data[,g]))
          rownames(des)=names(data)[i]
        }
        if(p[1]=="ANOVA"){
          m0=mean(data[,i],na.rm=T)
          s0=sd(data[,i],na.rm=T)
          ms0=paste(round(m0,3),"+/-",round(s0,3),sep="")
          m=aggregate(data[,i],by=list(data[,g]),mean,na.rm=TRUE)
          s=aggregate(data[,i],by=list(data[,g]),sd,na.rm=TRUE)
          ms=paste(as.matrix(round(m[,-1],3)),"+/-",as.matrix(round(s[,-1],3)),
                   sep="")
          ms=matrix(ms,byrow=T,nrow=1)
          des=data.frame(cbind(ms0,ms))
          names(des)=c("total",levels(data[,g]))
          rownames(des)=names(data)[i]
        }
        p=cbind(des,p)
        return(p)
      }

      if(length(j)>1){
        rst=aov.kru.single(data,j[1],g)
        for(jj in 2:length(j)){rst=rbind(rst,aov.kru.single(data,j[jj],g))}
      }

      if(length(j)==1){rst=aov.kru.single(data,j[1],g)}

      return(rst)
    }#if(!aov)
  }#aov.kru

  ####################################################
  ## function of factor variables from more groups ##
  ####################################################

  fisher.des.m=function(nc.g=nc.g,nc.f=nc.f,dataInS=dataIn,ftestS=fisher){

    fisher.des.msingle=function(g=nc.g,i,data=dataInS,test=ftestS){
      data[,g]=as.factor(data[,g])
      data[,i]=as.factor(data[,i])
      t=table(data[,i])
      p=prop.table(t)
      p=round(p*100,2)
      rst0=paste(t," (",p,")",sep="")
      t=table(data[,i],data[,g])
      p=prop.table(t,2)
      p=round(p*100,2)
      rst=data.frame(matrix(paste(t," (",p,")",sep=""),ncol=length(levels(data[,g]))))
      if(test){
        fisher=fisher.test(t)
        p=round(fisher$p.value,3)
        method=c("Fisher test",rep(NA,nrow(t)-1))
        statistic=rep(NA,nrow(t))
        p=c(p,rep(NA,nrow(t)-1))
        rst=cbind(rst0,rst,method,statistic,p)
        names(rst)=c("total",levels(data[,g]),"method","statistic","p.value")
        rownames(rst)=paste(names(data)[i],levels(data[,i]))

        pairwiseComb=combn(levels(data[,g]),2)
        results=matrix(NA,ncol=ncol(pairwiseComb),nrow=nrow(t))
        rownames(results)=rownames(rst)
        colnames(results)=apply(pairwiseComb,2,function(x){paste(x[1],"_vs_",x[2],sep="")})
        for(j in seq(ncol(pairwiseComb))){
          tempCol1=pairwiseComb[1,j]
          tempCol2=pairwiseComb[2,j]
          cols=c(tempCol1, tempCol2)
          tempMat=rbind(t[,cols])
          tempFisher=fisher.test(tempMat,alternative="two.sided")
          results[1,colnames(results)[j]]=tempFisher$p.value
        }#for(j in seq(ncol(pairwiseComb)))
        results[1,]=round(p.adjust(results[1,],method="fdr"),3)

        rst=cbind(rst,results)

      }#if(test)
      if(!test){
        chisq=chisq.test(t)
        p=round(chisq$p.value,3)
        method=c("Chi-square test",rep(NA,nrow(t)-1))
        statistic=c(chisq$statistic,rep(NA,nrow(t)-1))
        p=c(p,rep(NA,nrow(t)-1))
        rst=cbind(rst0,rst,method,statistic,p)
        names(rst)=c("total",levels(data[,g]),"method","statistic","p.value")
        rownames(rst)=paste(names(data)[i],levels(data[,i]))

        pairwiseComb=combn(levels(data[,g]),2)
        results=matrix(NA,ncol=ncol(pairwiseComb),nrow=nrow(t))
        rownames(results)=rownames(rst)
        colnames(results)=apply(pairwiseComb,2,function(x){paste(x[1],"_vs_",x[2],sep="")})
        for(j in seq(ncol(pairwiseComb))){
          tempCol1 <- pairwiseComb[1,j]
          tempCol2 <- pairwiseComb[2,j]
          cols <- c(tempCol1, tempCol2)
          tempMat=rbind(t[,cols])
          tempChisq=chisq.test(tempMat)
          results[1,colnames(results)[j]] <- tempChisq$p.value
        }#for(j in seq(ncol(pairwiseComb)))
        results[1,]=round(p.adjust(results[1,], method="fdr"),3)

        rst=cbind(rst,results)

      }#if(!test)
      return(rst)
    }#fisher.des.msingle

    f.rst=fisher.des.msingle(g=nc.g,i=nc.f[1],data=dataInS,test=ftestS)

    if(length(nc.f)==1){f.rst=f.rst}else{
      for(jj in 2:length(nc.f)){
        f.rst=rbind(f.rst,fisher.des.msingle(g=nc.g,i=nc.f[jj],data=dataInS,test=ftestS))
      }
    }
    for(i in 1:ncol(f.rst)){f.rst[,i]=as.character(f.rst[,i])}
    return(f.rst)
  }

  ###################
  ## main function ##
  ###################

  if(is.null(nc.g)){message("Please input the value of 'nc.g'!")}else{

    dataIn[,nc.g]=as.factor(dataIn[,nc.g])

    if(length(levels(dataIn[,nc.g]))==2){

      if(is.null(nc.n)){message("There is not any value of 'nc.n'!");n.rst=NULL}else{
        n.rst=t.wilcox.des(nc.g=nc.g,nc.n=nc.n,dataInS=dataIn)
      }

      if(is.null(nc.f)){message("There is not any value of 'nc.f'!");f.rst=NULL}else{
        f.rst=fisher.des(nc.g=nc.g,nc.f=nc.f,dataInS=dataIn,ftestS=fisher)
      }

      if(is.null(n.rst)){return(f.rst)}
      if(is.null(f.rst)){return(n.rst)}
      if(!is.null(n.rst) & !is.null(f.rst)){rst=rbind(n.rst,f.rst);return(rst)}

    }#if(length(levels(dataIn[,nc.g]))==2)

    if(length(levels(dataIn[,nc.g]))>2){

      if(is.null(nc.n)){message("There is not any value of 'nc.n'!");n.rst=NULL}else{
        n.rst=aov.kru(g=nc.g,j=nc.n,data=dataIn,aov=aov)
      }

      if(is.null(nc.f)){message("There is not any value of 'nc.f'!");f.rst=NULL}else{
        f.rst=fisher.des.m(nc.g=nc.g,nc.f=nc.f,dataInS=dataIn,ftestS=fisher)
      }

      if(is.null(n.rst)){return(f.rst)}
      if(is.null(f.rst)){return(n.rst)}
      if(!is.null(n.rst) & !is.null(f.rst)){rst=rbind(n.rst,f.rst);return(rst)}

    }#if(length(levels(dataIn[,nc.g]))>2){

  }#main function

}#easyDes
