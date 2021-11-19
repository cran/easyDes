# imports from multcomp::glht
# imports from PMCMR:posthoc.kruskal.nemenyi.test
# imports from PMCMRplus:kwAllPairsNemenyiTest

easyDes=function(
  nc.g=NULL,
  nc.n=NULL,
  nc.f=NULL,
  nc.of=NULL,
  dataIn=data,
  fisher=TRUE,
  aov=FALSE,
  t=FALSE,
  mean=FALSE,
  mcp.test.method="Tukey",
  mcp.stat=FALSE,
  mcp.t.test=FALSE,
  mcp.t.test.method="fdr",
  table.margin=2,
  decimal.p=3,
  decimal.prop=2){

  ###################################################
  ## function of numeric variables from two groups ##
  ###################################################

  t.wilcox.des <- function(nc.g=nc.g,nc.n=nc.n,dataInS=dataIn,t2=t,mean2=mean){

    if(t2==FALSE){

    t.wilcox.single <- function(g=nc.g,i,data=dataInS){
      data[,i]=as.numeric(data[,i])

      m0=mean(data[,i],na.rm=TRUE)
      s0=sd(data[,i],na.rm=TRUE)
      ms0=paste(sprintf("%.3f",round(m0,3)),"+/-",sprintf("%.3f",round(s0,3)),sep="")
      m=aggregate(data[,i],by=list(data[,g]),mean,na.rm=TRUE)
      s=aggregate(data[,i],by=list(data[,g]),sd,na.rm=TRUE)
      des.mean=paste(sprintf("%.3f",round(m[,2],3)),"+/-",sprintf("%.3f",round(s[,2],3)),sep="")

      med0=median(data[,i],na.rm=TRUE)
      qs0=quantile(data[,i],na.rm=TRUE)
      medq0=paste(sprintf("%.3f",round(med0,3)),"(",sprintf("%.3f",round(qs0[2],3)),",",sprintf("%.3f",round(qs0[4],3)),")")
      med=aggregate(data[,i],by=list(data[,g]),median,na.rm=TRUE)
      qs=aggregate(data[,i],by=list(data[,g]),quantile,na.rm=TRUE)
      des.median=paste(sprintf("%.3f",round(med[,2],3)),"(",sprintf("%.3f",round(qs[,2][,2],3)),",",sprintf("%.3f",round(qs[,2][,4],3)),")")

      if(length((data[,i][data[,g]==levels(data[,g])[1]]))>5000 |
         length((data[,i][data[,g]==levels(data[,g])[2]]))>5000){
        shapiro1=ks.test(data[,i][data[,g]==levels(data[,g])[1]],pnorm)
        shapiro2=ks.test(data[,i][data[,g]==levels(data[,g])[2]],pnorm)
      }else{
        shapiro1=shapiro.test(data[,i][data[,g]==levels(data[,g])[1]])
        shapiro2=shapiro.test(data[,i][data[,g]==levels(data[,g])[2]])
      }
      shapiro=min(shapiro1$p.value,shapiro2$p.value)
      bartlett=bartlett.test(data[,i]~data[,g])
      if(shapiro<0.1 | bartlett$p.value<0.1){
        temp=wilcox.test(data[,i]~data[,g],exact=FALSE)
        s=c("Wilcoxon test",sprintf("%.3f",round(temp$statistic,3)))
        p=sprintf(paste0("%.",decimal.p,"f"),temp$p.value)
        if(mean){rst.single=c(ms0,des.mean,s,p)}else{rst.single=c(medq0,des.median,s,p)}
        rst.single=data.frame(matrix(rst.single,nrow=1))}
      if(shapiro>=0.1 & bartlett$p.value<0.1){
        temp=t.test(data[,i]~data[,g],var.equal=FALSE)
        s=c("Adjust t test",sprintf("%.3f",round(temp$statistic,3)))
        p=sprintf(paste0("%.",decimal.p,"f"),temp$p.value)
        rst.single=c(ms0,des.mean,s,p)
        rst.single=data.frame(matrix(rst.single,nrow=1))}
      if(shapiro>=0.1 & bartlett$p.value>=0.1){
        temp=t.test(data[,i]~data[,g],var.equal=TRUE)
        s=c("t test",sprintf("%.3f",round(temp$statistic,3)))
        p=sprintf(paste0("%.",decimal.p,"f"),temp$p.value)
        rst.single=c(ms0,des.mean,s,p)
        rst.single=data.frame(matrix(rst.single,nrow=1))}
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
        ms0=paste(sprintf("%.3f",round(m0,3)),"+/-",sprintf("%.3f",round(s0,3)),sep="")
        m=aggregate(data[,i],by=list(data[,g]),mean,na.rm=TRUE)
        s=aggregate(data[,i],by=list(data[,g]),sd,na.rm=TRUE)
        rst.single=paste(sprintf("%.3f",round(m[,2],3)),"+/-",sprintf("%.3f",round(s[,2],3)),sep="")
        temp=t.test(data[,i]~data[,g],var.equal=TRUE)
        s=c("t test",sprintf("%.3f",round(temp$statistic,3)))
        p=sprintf(paste0("%.",decimal.p,"f"),temp$p.value)
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

  ####################################################
  ## function of numeric variables from more groups ##
  ####################################################

  aov.kru=function(g=nc.g,j=nc.n,data=dataIn,aov=aov){

    #fix the bug with the data having the 'group' variable already
    names(data)[names(data)=="group"]="group_temp"

    if(aov){

      aov.single=function(data,i,g){

        #fix the bug with the variable having the special character
        nc.name=names(data)[i]
        names(data)[i]="nc.key.name.only"

        data[,g]=as.factor(data[,g])
        data[,i]=as.numeric(data[,i])
        names(data)[g]="group"
        g.lev=levels(data[,g])
        p=NA
        model=formula(paste(names(data)[i],"~ group"))
        temp=aov(model,data)
        p[1]="ANOVA"
        p[2]=sprintf("%.3f",round(summary(temp)[[1]][1,4],3))
        p[3]=sprintf(paste0("%.",decimal.p,"f"),summary(temp)[[1]][1,5])

        if(mcp.t.test==TRUE){
          mcp.stat=FALSE
          temp=pairwise.t.test(data[,i],data[,g],p.adjust.method=mcp.t.test.method)
          temp=as.numeric(temp$p.value)[!is.na(as.numeric(temp$p.value))]
        }
        if(mcp.t.test==FALSE){
          temp=summary(glht(temp,linfct=mcp(group=mcp.test.method)))
          attr(temp$test$pvalues,"error")=NULL
          temp=temp$test$pvalues
        }

        if(mcp.stat==FALSE){
          p[4:(3+choose(length(g.lev),2))]=sprintf("%.3f",round(temp,3))
          p=data.frame(t(p))
          r.name=NULL
          for(k in 1:(length(g.lev)-1)){r.name=c(r.name,paste(g.lev[k],"_vs_",g.lev[(k+1):length(g.lev)],sep=""))}
          rownames(p)=nc.name
          colnames(p)=c("method","statistic","p.value",r.name)

          m0=mean(data[,i],na.rm=T)
          s0=sd(data[,i],na.rm=T)
          ms0=paste(sprintf("%.3f",round(m0,3)),"+/-",sprintf("%.3f",round(s0,3)),sep="")
          m=aggregate(data[,i],by=list(data[,g]),mean,na.rm=TRUE)
          s=aggregate(data[,i],by=list(data[,g]),sd,na.rm=TRUE)
          ms=paste(as.matrix(sprintf("%.3f",round(m[,-1],3))),"+/-",as.matrix(sprintf("%.3f",round(s[,-1],3))),sep="")
          ms=matrix(ms,byrow=T,nrow=1)
          des=data.frame(cbind(ms0,ms))
          names(des)=c("total",levels(data[,g]))
          rownames(des)=nc.name
          p=cbind(des,p)
        }

        if(mcp.stat){
          temp.s=sprintf("%.3f",temp$test$tstat)
          temp.p=sprintf(paste0("%.",decimal.p,"f"),temp$test$pvalues)
          p[4:(3+length(temp.s)+length(temp.p))]=c(temp.p,temp.s)
          p=data.frame(t(p))
          r.name=NULL
          for(k in 1:(length(g.lev)-1)){
            r.name=c(r.name,paste(g.lev[k],"_vs_",g.lev[(k+1):length(g.lev)],sep=""))
          }
          rownames(p)=nc.name
          colnames(p)=c("method","statistic","p.value",
                        paste("p.",r.name,sep=""),paste("stat.",r.name,sep=""))
          m0=mean(data[,i],na.rm=T)
          s0=sd(data[,i],na.rm=T)
          ms0=paste(sprintf("%.3f",round(m0,3)),"+/-",sprintf("%.3f",round(s0,3)),sep="")
          m=aggregate(data[,i],by=list(data[,g]),mean,na.rm=TRUE)
          s=aggregate(data[,i],by=list(data[,g]),sd,na.rm=TRUE)
          ms=paste(as.matrix(sprintf("%.3f",round(m[,-1],3))),"+/-",as.matrix(sprintf("%.3f",round(s[,-1],3))),sep="")
          ms=matrix(ms,byrow=T,nrow=1)
          des=data.frame(cbind(ms0,ms))
          names(des)=c("total",levels(data[,g]))
          rownames(des)=nc.name
          p=cbind(des,p)

        } #if(mcp.stat)

        return(p)
      }#aov.single

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

        #fix the bug with the variable having the special character
        nc.name=names(data)[i]
        names(data)[i]="nc.key.name.only"

        data[,g]=as.factor(data[,g])
        data[,i]=as.numeric(data[,i])
        names(data)[g]="group"
        g.lev=levels(data[,g])
        p.shapiro=NA
        for(g.num in 1:length(g.lev)){
          if(length(data[,i][data[,g]==g.lev[g.num]])>5000){
            shapiro=ks.test(data[,i][data[,g]==g.lev[g.num]])
          }else{
            shapiro=shapiro.test(data[,i][data[,g]==g.lev[g.num]])
          }
          p.shapiro[g.num]=shapiro$p.value
        }
        bartlett=bartlett.test(data[,i]~data[,g])
        p=NA
        model=formula(paste(names(data)[i],"~ group"))
        if(!mcp.stat){

          if(min(p.shapiro)<0.1 | bartlett$p.value<0.1){
            temp=kruskal.test(model,data)
            p[1]="Kruskal-Wallis test"
            p[2]=sprintf("%.3f",round(temp$statistic,3))
            p[3]=sprintf(paste0("%.",decimal.p,"f"),temp$p.value)
            temp=kwAllPairsNemenyiTest(x=data[,i],g=data[,g],method="Chisq")#Tukey
            temp.p=as.numeric(temp$p.value)
            temp.p=temp.p[!is.na(temp.p)]
            p[4:(3+choose(length(g.lev),2))]=sprintf(paste0("%.",decimal.p,"f"),temp.p)
          } #if(min(p.shapiro)<0.1 | bartlett$p.value<0.1)

          if(min(p.shapiro)>=0.1 &  bartlett$p.value>=0.1){
            temp=aov(model,data)
            p[1]="ANOVA"
            p[2]=sprintf("%.3f",round(summary(temp)[[1]][1,4],3))
            p[3]=sprintf(paste0("%.",decimal.p,"f"),summary(temp)[[1]][1,5])
            if(mcp.t.test==TRUE){
              mcp.stat=FALSE
              temp=pairwise.t.test(data[,i],data[,g],p.adjust.method=mcp.t.test.method)
              temp.p=as.numeric(temp$p.value)[!is.na(as.numeric(temp$p.value))]
            }
            if(mcp.t.test==FALSE){
              temp=summary(glht(temp,linfct=mcp(group=mcp.test.method)))
              temp.p=temp$test[[6]]
            }
            p[4:(3+choose(length(g.lev),2))]=sprintf(paste0("%.",decimal.p,"f"),temp.p)
          }#if(min(p.shapiro)>=0.1 &  bartlett$p.value>=0.1){

          p=data.frame(t(p))
          r.name=NULL
          for(k in 1:(length(g.lev)-1)){r.name=c(r.name,paste(g.lev[k],"_vs_",g.lev[(k+1):length(g.lev)],sep=""))}
          rownames(p)=nc.name
          colnames(p)=c("method","statistic","p.value",r.name)

        }

        if(mcp.stat){

          if(min(p.shapiro)<0.1 | bartlett$p.value<0.1){
            temp=kruskal.test(model,data)
            p[1]="Kruskal-Wallis test"
            p[2]=sprintf("%.3f",round(temp$statistic,3))
            p[3]=sprintf(paste0("%.",decimal.p,"f"),temp$p.value)
            temp=kwAllPairsNemenyiTest(x=data[,i],g=data[,g],method="Chisq")#Tukey
            temp.s=as.numeric(temp$statistic)
            temp.p=as.numeric(temp$p.value)
            temp.s=sprintf("%.3f",temp.s[!is.na(temp.s)])
            temp.p=sprintf(paste0("%.",decimal.p,"f"),temp.p[!is.na(temp.p)])
            p[4:(3+length(temp.s)+length(temp.p))]=c(temp.p,temp.s)
          } #if(min(p.shapiro)<0.1 | bartlett$p.value<0.1)

          if(min(p.shapiro)>=0.1 &  bartlett$p.value>=0.1){
            temp=aov(model,data)
            p[1]="ANOVA"
            p[2]=sprintf("%.3f",round(summary(temp)[[1]][1,4],3))
            p[3]=sprintf(paste0("%.",decimal.p,"f"),summary(temp)[[1]][1,5])
            #"Tukey" or "Dunnett" method
            temp=summary(glht(temp,linfct=mcp(group=mcp.test.method)))
            temp.s=temp$test$tstat
            temp.p=temp$test$pvalues
            p[4:(3+length(temp.s)+length(temp.p))]=sprintf("%.3f",round(c(temp.p,temp.s),3))
          } #if(min(p.shapiro)>=0.1 &  bartlett$p.value>=0.1)

          p=data.frame(t(p))
          r.name=NULL
          for(k in 1:(length(g.lev)-1)){
            r.name=c(r.name,paste(g.lev[k],"_vs_",g.lev[(k+1):length(g.lev)],sep=""))
          }
          rownames(p)=nc.name
          colnames(p)=c("method","statistic","p.value",
                        paste("p.",r.name,sep=""),paste("stat.",r.name,sep=""))

        } #if(mcp.stat)

        if(p[1]=="Kruskal-Wallis test" & !mean){
          q0=quantile(data[,i],na.rm=T)
          qr0=paste(sprintf("%.3f",round(q0[3],3)),"(",sprintf("%.3f",round(q0[2],3)),",",sprintf("%.3f",round(q0[4],3)),")",sep="")
          q=aggregate(data[,i],by=list(data[,g]),quantile,na.rm=TRUE)
          qr=paste(sprintf("%.3f",round(q[,-1][,3],3)),"(",sprintf("%.3f",round(q[,-1][,2],3)),",",sprintf("%.3f",round(q[,-1][,4],3)),")",sep="")
          qr=matrix(qr,byrow=T,nrow=1)
          des=data.frame(cbind(qr0,qr))
          names(des)=c("total",levels(data[,g]))
          rownames(des)=nc.name
        }
        if(p[1]=="ANOVA" | mean){
          m0=mean(data[,i],na.rm=T)
          s0=sd(data[,i],na.rm=T)
          ms0=paste(sprintf("%.3f",round(m0,3)),"+/-",sprintf("%.3f",round(s0,3)),sep="")
          m=aggregate(data[,i],by=list(data[,g]),mean,na.rm=TRUE)
          s=aggregate(data[,i],by=list(data[,g]),sd,na.rm=TRUE)
          ms=paste(as.matrix(sprintf("%.3f",round(m[,-1],3))),"+/-",as.matrix(sprintf("%.3f",round(s[,-1],3))),sep="")
          ms=matrix(ms,byrow=T,nrow=1)
          des=data.frame(cbind(ms0,ms))
          names(des)=c("total",levels(data[,g]))
          rownames(des)=nc.name
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

  ##########################################################
  ## function of factor variables from two or more groups ##
  ##########################################################

  fisher.des=function(nc.g=nc.g,nc.f=nc.f,dataInS=dataIn,ftestS=fisher){

    fisher.des.single=function(g=nc.g,i,data=dataInS,test=ftestS){
      data[,g]=as.factor(data[,g])
      data[,i]=as.factor(data[,i])
      t=table(data[,i])
      p=prop.table(t)
      p=sprintf(paste0("%.",decimal.prop,"f"),p*100)
      rst0=paste(t," (",p,")",sep="")
      t=table(data[,i],data[,g])
      p=prop.table(t,margin=table.margin)
      p=sprintf(paste0("%.",decimal.prop,"f"),p*100)
      rst=data.frame(matrix(paste(t," (",p,")",sep=""),ncol=length(levels(data[,g]))))

      if(test){

        fisher=try(fisher.test(t),silent=TRUE)
        if(inherits(fisher,"try-error")){
          fisher=fisher.test(t,simulate.p.value=TRUE)
        }else{
          fisher=fisher.test(t)
        }

        p=sprintf(paste0("%.",decimal.p,"f"),fisher$p.value)
        method=c("Fisher test",rep(NA,nrow(t)-1))
        statistic=rep(NA,nrow(t))
        p=c(p,rep(NA,nrow(t)-1))
        rst=cbind(rst0,rst,method,statistic,p)
        names(rst)=c("total",levels(data[,g]),"method","statistic","p.value")
        rownames(rst)=paste(names(data)[i],levels(data[,i]))

        if(length(levels(data[,g]))>2){
          pairwiseComb=combn(levels(data[,g]),2)
          results=matrix(NA,ncol=ncol(pairwiseComb),nrow=nrow(t))
          rownames(results)=rownames(rst)
          colnames(results)=apply(pairwiseComb,2,function(x){paste(x[1],"_vs_",x[2],sep="")})
          for(j in seq(ncol(pairwiseComb))){
            tempCol1=pairwiseComb[1,j]
            tempCol2=pairwiseComb[2,j]
            cols=c(tempCol1,tempCol2)
            tempMat=rbind(t[,cols])
            tempFisher=try(fisher.test(tempMat,alternative="two.sided"),silent=TRUE)
            if(inherits(tempFisher,"try-error")){
              tempFisher=fisher.test(tempMat,alternative="two.sided",simulate.p.value=TRUE)
            }else{
              tempFisher=fisher.test(tempMat,alternative="two.sided")
            }
            results[1,colnames(results)[j]]=tempFisher$p.value
          }#for(j in seq(ncol(pairwiseComb)))
          results[1,]=sprintf(paste0("%.",decimal.p,"f"),p.adjust(results[1,],method="fdr"))

          if(!mcp.stat){
            rst=cbind(rst,results)
          }

          if(mcp.stat){
            s.results=results
            s.results[1,]=NA
            colnames(s.results)=paste("stat.",colnames(results),sep="")
            colnames(results)=paste("p.",colnames(results),sep="")
            rst=cbind(rst,results,s.results)
          }

        } #if(length(levels(data[,g]))>2)

      } #if(test)

      if(!test){
        chisq=chisq.test(t)
        p=sprintf(paste0("%.",decimal.p,"f"),chisq$p.value)
        method=c("Chi-square test",rep(NA,nrow(t)-1))
        statistic=c(sprintf("%.3f",round(chisq$statistic,3)),rep(NA,nrow(t)-1))
        p=c(p,rep(NA,nrow(t)-1))
        rst=cbind(rst0,rst,method,statistic,p)
        names(rst)=c("total",levels(data[,g]),"method","statistic","p.value")
        rownames(rst)=paste(names(data)[i],levels(data[,i]))

        if(length(levels(data[,g]))>2){
          pairwiseComb=combn(levels(data[,g]),2)
          results=matrix(NA,ncol=ncol(pairwiseComb),nrow=nrow(t))
          s.results=matrix(NA,ncol=ncol(pairwiseComb),nrow=nrow(t))
          rownames(results)=rownames(rst)
          rownames(s.results)=rownames(rst)
          colnames(results)=apply(pairwiseComb,2,function(x){paste(x[1],"_vs_",x[2],sep="")})
          colnames(s.results)=apply(pairwiseComb,2,function(x){paste(x[1],"_vs_",x[2],sep="")})
          for(j in seq(ncol(pairwiseComb))){
            tempCol1 <- pairwiseComb[1,j]
            tempCol2 <- pairwiseComb[2,j]
            cols <- c(tempCol1, tempCol2)
            tempMat=rbind(t[,cols])
            tempChisq=chisq.test(tempMat)
            results[1,colnames(results)[j]] <- tempChisq$p.value
            s.results[1,colnames(s.results)[j]] <- tempChisq$statistic
          }#for(j in seq(ncol(pairwiseComb)))
          results[1,]=sprintf(paste0("%.",decimal.p,"f"),p.adjust(results[1,], method="fdr"))
          s.results[1,]=sprintf("%.3f",round(s.results[1,],3))

          if(!mcp.stat){
            rst=cbind(rst,results)
          }

          if(mcp.stat){
            colnames(s.results)=paste("stat.",colnames(s.results),sep="")
            colnames(results)=paste("p.",colnames(results),sep="")
            rst=cbind(rst,results,s.results)
          }

        }#if(length(levels(data[,g]))>2)

      }#if(!test)
      return(rst)
    } #fisher.des.single

    f.rst=fisher.des.single(g=nc.g,i=nc.f[1],data=dataInS,test=ftestS)

    if(length(nc.f)==1){f.rst=f.rst}else{
      for(jj in 2:length(nc.f)){
        f.rst=rbind(f.rst,fisher.des.single(g=nc.g,i=nc.f[jj],data=dataInS,test=ftestS))
      }
    }
    for(i in 1:ncol(f.rst)){f.rst[,i]=as.character(f.rst[,i])}
    return(f.rst)
  } #fisher.des=function()

  ##################################
  ## function for ordinal factors ##
  ##################################

  of.2g=function(nc.g,nc.of,dataInS){
    dataInS[,nc.g]=as.factor(dataInS[,nc.g])
    rst2=NULL
    for(i in 1:length(nc.of)){
      if(class(dataInS[,nc.of[i]])!="factor"){stop("Some variables in columns of nc.of are not factors!")}
      rnames=paste(names(dataInS)[nc.of[i]],levels(dataInS[,nc.of[i]]))
      t=table(dataInS[,nc.of[i]])
      p=prop.table(t)
      p=sprintf(paste0("%.",decimal.prop,"f"),p*100)
      rst0=paste(t," (",p,")",sep="")
      t=table(dataInS[,nc.of[i]],dataInS[,nc.g])
      p=prop.table(t,margin=table.margin)
      p=sprintf(paste0("%.",decimal.prop,"f"),p*100)
      rst=cbind(rst0,data.frame(matrix(paste(t," (",p,")",sep=""),ncol=length(levels(dataInS[,nc.g])))))

      temp=wilcox.test(as.numeric(dataInS[,nc.of[i]])~dataInS[,nc.g],exact=FALSE)
      s=sprintf("%.3f",round(temp$statistic,3))
      p=sprintf(paste0("%.",decimal.p,"f"),temp$p.value)
      rst=cbind(rst,data.frame(method=c("Wilcoxon test",rep(NA,nrow(rst)-1)),
                               statistic=c(s,rep(NA,nrow(rst)-1)),
                               p=c(p,rep(NA,nrow(rst)-1))))
      rownames(rst)=rnames
      rst2=rbind(rst2,rst)

    }
   names(rst2)=c("total",levels(dataInS[,nc.g]),"method","statistic","p.value")
   return(rst2)
  }

  of.3g=function(nc.g,nc.of,dataInS,mcp.stat){
    dataInS[,nc.g]=as.factor(dataInS[,nc.g])
    names(dataInS)[nc.g]="group"
    pairwiseComb=combn(levels(dataInS[,nc.g]),2)
    pairwise.name=apply(pairwiseComb,2,function(x){paste(x[1],"_vs_",x[2],sep="")})
    rst2=NULL
    for(i in 1:length(nc.of)){
      if(class(dataInS[,nc.of[i]])!="factor"){stop("Some variables in columns of nc.of are not factors!")}
      rnames=paste(names(dataInS)[nc.of[i]],levels(dataInS[,nc.of[i]]))
      t=table(dataInS[,nc.of[i]])
      p=prop.table(t)
      p=sprintf(paste0("%.",decimal.prop,"f"),p*100)
      rst0=paste(t," (",p,")",sep="")
      t=table(dataInS[,nc.of[i]],dataInS[,nc.g])
      p=prop.table(t,margin=table.margin)
      p=sprintf(paste0("%.",decimal.prop,"f"),p*100)
      rst=cbind(rst0,data.frame(matrix(paste(t," (",p,")",sep=""),ncol=length(levels(dataInS[,nc.g])))))

      dataInS[,nc.of[i]]=as.numeric(dataInS[,nc.of[i]])
      model=formula(paste(names(dataInS)[nc.of[i]],"~ group"))
      temp=kruskal.test(model,data=dataInS)
      s=sprintf("%.3f",round(temp$statistic,3))
      p=sprintf("%.3f",round(temp$p.value,3))
      rst=cbind(rst,data.frame(method=c("Kruskal-Wallis test",rep(NA,nrow(rst)-1)),
                               statistic=c(s,rep(NA,nrow(rst)-1)),
                               p=c(p,rep(NA,nrow(rst)-1))))

      temp=suppressWarnings(kwAllPairsNemenyiTest(x=dataInS[,nc.of[i]],g=dataInS[,nc.g],method="Chisq"))#Tukey)
      temp.s=as.numeric(temp$statistic)
      temp.s=temp.s[!is.na(temp.s)]
      temp.s=sprintf("%.3f",round(temp.s,3))
      temp.s2=matrix(rep(NA,nrow(rst)*length(temp.s)),ncol=length(temp.s),byrow=TRUE)
      temp.s2[1,]=temp.s
      temp.p=as.numeric(temp$p.value)
      temp.p=temp.p[!is.na(temp.p)]
      temp.p=sprintf(paste0("%.",decimal.p,"f"),temp.p)
      temp.p2=matrix(rep(NA,nrow(rst)*length(temp.p)),ncol=length(temp.p),byrow=TRUE)
      temp.p2[1,]=temp.p

      if(mcp.stat==TRUE){rst=cbind(rst,temp.p2,temp.s2)}
      if(mcp.stat==FALSE){rst=cbind(rst,temp.p2)}
      row.names(rst)=rnames
      rst2=rbind(rst2,rst)
    }
    if(mcp.stat==TRUE){
      names(rst2)=c("total",levels(dataInS[,nc.g]),"method","statistic","p.value",
                    paste0("p.",pairwise.name),
                    paste0("stat.",pairwise.name))
      }
    if(mcp.stat==FALSE){
      names(rst2)=c("total",levels(dataInS[,nc.g]),"method","statistic","p.value",pairwise.name)
      }
    return(rst2)
  }

  ###################
  ## main function ##
  ###################

  if(is.null(nc.g)){message("Please input the value of 'nc.g'!")}else{

    dataIn[,nc.g]=as.factor(dataIn[,nc.g])

    if(length(levels(dataIn[,nc.g]))==2){

      if(is.null(nc.n)){message("There is not any value of 'nc.n'!");n.rst=NULL}else{
        n.rst=t.wilcox.des(nc.g=nc.g,nc.n=nc.n,dataInS=dataIn,t2=t)
      }

      if(is.null(nc.f)){message("There is not any value of 'nc.f'!");f.rst=NULL}else{
        f.rst=fisher.des(nc.g=nc.g,nc.f=nc.f,dataInS=dataIn,ftestS=fisher)
      }

      if(is.null(nc.of)){message("There is not any value of 'nc.of'!");of.rst=NULL}else{
        of.rst=of.2g(nc.g=nc.g,nc.of=nc.of,dataInS=dataIn)
      }

    }#if(length(levels(dataIn[,nc.g]))==2)

    if(length(levels(dataIn[,nc.g]))>2){

      if(is.null(nc.n)){message("There is not any value of 'nc.n'!");n.rst=NULL}else{
        n.rst=aov.kru(g=nc.g,j=nc.n,data=dataIn,aov=aov)
      }

      if(is.null(nc.f)){message("There is not any value of 'nc.f'!");f.rst=NULL}else{
        f.rst=fisher.des(nc.g=nc.g,nc.f=nc.f,dataInS=dataIn,ftestS=fisher)
      }

      if(is.null(nc.of)){message("There is not any value of 'nc.of'!");of.rst=NULL}else{
        of.rst=of.3g(nc.g=nc.g,nc.of=nc.of,dataInS=dataIn,mcp.stat=mcp.stat)
      }

    }#if(length(levels(dataIn[,nc.g]))>2){

    if(is.null(f.rst) & is.null(of.rst)){return(n.rst)}
    if(is.null(n.rst) & is.null(of.rst)){return(f.rst)}
    if(is.null(n.rst) & is.null(f.rst)){return(of.rst)}
    if(!is.null(n.rst) & !is.null(f.rst) & is.null(of.rst)){rst=rbind(n.rst,f.rst);return(rst)}
    if(!is.null(n.rst) & !is.null(of.rst) & is.null(f.rst)){rst=rbind(n.rst,of.rst);return(rst)}
    if(!is.null(f.rst) & !is.null(of.rst) & is.null(n.rst)){rst=rbind(f.rst,of.rst);return(rst)}
    if(!is.null(n.rst) & !is.null(f.rst) & !is.null(of.rst)){rst=rbind(n.rst,f.rst,of.rst);return(rst)}

  }#main function

}#easyDes
