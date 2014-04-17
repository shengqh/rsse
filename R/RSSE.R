options("guiToolkit"="RGtk2")

est_pvalue_store<-function(X1X0Sum, n, phi, w=1) {
  a<-0:X1X0Sum
  if(w==1){
    anom<-dnbinom(a, mu=X1X0Sum/2, size=n/phi)
    zanom<-rev(anom)
  }else{
    anom<-dnbinom(a, mu=(w*X1X0Sum/(w+1)), size=n/phi)
    zanom<-dnbinom( X1X0Sum-a, mu=(X1X0Sum/(w+1)), size=n/phi)
  }
  azanom<-anom * zanom
  yy<-azanom/sum(azanom)
  return(yy)
}

est_pvalue_store_get_power3<-function(x1,x0,yy) {
  yy<-yy[[(x1+x0+1)]]
  if (x1>=x0) {
    pvalue<-2* sum(yy[(x1+1):length(yy)])
  } else {
    pvalue<-2*sum(yy[1:(x1+1)])
  }
  return(min(pvalue, 1))
}

est_pvalue<-function(x1,x0,n, phi, w=1) {
  yy<-est_pvalue_store(X1X0Sum=x1+x0,n=n, phi=phi, w=w)
  if (x1>=x0) {
    pvalue<-2* sum(yy[(x1+1):length(yy)])
  } else {
    pvalue<-2*sum(yy[1:(x1+1)])
  }
  return(min(pvalue, 1))
}

est_power3<-function(n, w=1, rho=2.0, mu0=5, phi_0=1, beta=0.2, alpha=0.05, error=0.001){
  mu0<-mu0
  mu1<-mu0*(rho*w)
  phi_1<-phi_0
  q0_u<-qnbinom(1-error, size=n/phi_0, mu=n*mu0)
  q0_l<-qnbinom(error, size=n/phi_0, mu=n*mu0)
  q1_u<-qnbinom(1-error, size=n/phi_1, mu=n*mu1)
  q1_l<-qnbinom(error, size=n/phi_1, mu=n*mu1)
  
  a<-0
  temp1<-dnbinom(q1_l:q1_u, mu=(n*mu1), size=n/phi_1)
  temp2<-dnbinom(q0_l:q0_u, mu=(n*mu0), size=n/phi_0)
  if (max(q0_u,q0_l,q1_u,q1_l)>=10000) { #Method2, doesn't store every pvalue but do it every time
    getPvalue<-function(x1,x0,yy=yyStore,...) {
      est_pvalue(x1,x0,n=n, phi=phi_0, w=w,...)
    }
  } else { #Method1, store every pvalue
    yMin<-min(q1_l,q1_u)+min(q0_l,q0_u)
    yMax<-max(q1_l,q1_u)+max(q0_l,q0_u)
    yyStore<-list()
    for (y in yMin:yMax) {
      yyStore[[y+1]]<-est_pvalue_store(y,n=n, phi=phi_0, w=w)
    }
    getPvalue<-function(x1,x0,yy=yyStore,...) {
      est_pvalue_store_get_power3(x1,x0,yy=yyStore,...)
    }
  }
  
  aNRow<-length(q1_l:q1_u)
  aNCol<-length(q0_l:q0_u)
  q0_l_loop<-q0_l
  q0_u_loop<-q0_u
  q1_l_loop<-q1_l
  q1_u_loop<-q1_u
  for (x in q0_l_loop:q0_u_loop) {
    for (y in q1_l_loop:q1_u_loop) {
      if (x>y) {break;}
      temp<-getPvalue(x1=y,x0=x)
      if (temp<=alpha) {
        a<-a+sum(temp1[(y-q1_l+1):aNRow]*temp2[(x-q0_l+1)])
        q1_l_loop<-y
        break;
      }
    }
  }
  q0_l_loop<-q0_l
  q0_u_loop<-q0_u
  q1_l_loop<-q1_l
  q1_u_loop<-q1_u
  for (y in q1_l_loop:q1_u_loop) {
    for (x in q0_l_loop:q0_u_loop) {
      if (x<y) {break;}
      temp<-getPvalue(x1=y,x0=x)
      if (temp<=alpha) {
        a<-a+sum(temp1[(y-q1_l+1)]*temp2[(x-q0_l+1):aNCol])
        q0_l_loop<-x
        break;
      }
    }
  }
  power<-a
  return(power-(1-beta))
}

##' sample_size
##' 
##' A function to estitamete the sample size for differential expression analysis of RNA-seq data.
##' 
##' @param m Total number of genes for testing.
##' @param m1 Expected number of prognostic genes.
##' @param power Power to detecte prognostic genes.
##' @param f FDR level
##' @param w Ratio of normalization factors between two groups.
##' @param rho minimum fold changes for prognostic genes between two groups.
##' @param mu0 Average read counts for prognostic genes.
##' @param phi_0 Dispersion for prognostic genes.
##' @param showMessage Logical. Display the message in the estimation process.
##' @importFrom ssanv uniroot.integer
##' @export
##' @examples 
##' \dontrun{#Input initial parameters. 
##' sample_size(m=10000, m1=100, power=0.8, f=0.1, w=1, rho=2, mu0=5, phi_0=1)
##' }
sample_size<-function(m=10000, m1=100, power=0.8, f=0.1, w=1, rho=2, mu0=5, phi_0=1,showMessage=F){
  r1<-m1 * power
  beta<-1-power
  step.power<-5
  
  alpha_star<-r1*f/((m-m1)*(1-f))
  z_alpha<-qnorm(1-alpha_star/2, lower.tail=T)
  z_beta<-qnorm(power, lower.tail=T)
  n_w<-( ( z_alpha + z_beta )^2* (1 + rho/w +phi_0*mu0*(1+rho^2)) )/ ( (rho-1)^2*mu0 )
  end.point<-round(n_w)+30
  
  start.point<-1

  p1<-est_power3(n=start.point,w=w, rho=rho, mu0=mu0, phi_0=phi_0, beta=beta, alpha=alpha_star)
  if (p1>0) {
    return(start.point)
  } else {
    if (p1<=0.8) {
      step.power<-7
    } else if (p1<=0.6) {
      step.power<-6
    }
    n_Exact<-ssanv::uniroot.integer(est_power3, c(start.point, end.point), w=w, rho=rho, 
                             mu0=mu0, phi_0=phi_0, beta=beta, alpha=alpha_star, pos.side=T,
                             step.up=T, step.power=step.power,print.steps=showMessage)$root  
    return(n_Exact)
  }
}

newlabel<-function(caption, fontlist, container){
  result<-glabel(caption, container=container)
  font(result)<-fontlist
  result
}

newedit<-function(defaultvalue, fontlist, container, width=10){
  result<-gedit(defaultvalue, width=width, container=container)
  font(result)<-fontlist
  result
}

##' rsse_ui
##' 
##' A function to display user friendly interface to do estitamete the sample size for differential expression analysis of RNA-seq data.
##' 
##' @export
##' @import gWidgets
##' @import gWidgetsRGtk2
##' @examples 
##' \dontrun{
##' rsse_ui()
##' }
rsse_ui<-function(){
  win <- gWidgets::gwindow("RNASeq Sample Size Estimation 0.98.3",width=1024, height=768, visible=FALSE)
  topgroup<-gWidgets::ggroup(container=win,horizontal=F)
  fontlist <- list(size=14)
  redfontlist <- list(size=14, color="red")
  bluefontlist <- list(size=14, color="blue")
  
  parameters <- gframe("Parameters", container=topgroup, horizontal=F)
  font(parameters)<-fontlist
  tbl <- glayout(container = parameters)
  
  tbl[1,1] <- (mLabel1<-newlabel("m:", fontlist, tbl))
  tbl[1,2] <- (mEdit<-newedit("10000", fontlist, tbl))
  tbl[1,3] <- (mLabel2<-newlabel("total number of genes for testing", fontlist, tbl))
  
  tbl[2,1] <- (m1Label1<-newlabel("m1:", fontlist, tbl))
  tbl[2,2] <- (m1Edit<-newedit("100", fontlist, tbl))
  tbl[2,3] <- (m1Label2<-newlabel("expected number of prognostic genes", fontlist, tbl))
  
  tbl[3,1] <- (powerLabel1<-newlabel("power:", fontlist, tbl))
  tbl[3,2] <- (powerEdit<-newedit("0.8", fontlist, tbl))
  tbl[3,3] <- (powerLabel2<-newlabel("power to detect prognostic genes", fontlist, tbl))
  
  tbl[4,1] <- (fLabel1<-newlabel("f:", fontlist, tbl))
  tbl[4,2] <- (fEdit<-newedit("0.1", fontlist, tbl))
  tbl[4,3] <- (fLabel2<-newlabel("FDR level", fontlist, tbl))
  
  tbl[5,1] <- (wLabel1<-newlabel("w:", fontlist, tbl))
  tbl[5,2] <- (wEdit<-newedit("1.0", fontlist, tbl))
  tbl[5,3] <- (wLabel2<-newlabel("ratio of normalization factors between two groups", fontlist, tbl))
  
  tbl[6,1] <- (rhoLabel1<-newlabel("rho:", fontlist, tbl))
  tbl[6,2] <- (rhoEdit<-newedit("2.0", fontlist, tbl))
  tbl[6,3] <- (rhoLabel2<-newlabel("minimum fold changes for prognostic genes in control group", fontlist, tbl))
  
  tbl[7,1] <- (mu0Label1<-newlabel("mu0:", fontlist, tbl))
  tbl[7,2] <- (mu0Edit<-newedit("5", fontlist, tbl))
  tbl[7,3] <- (mu0Label2<-newlabel("minimum average read counts among the prognostic genes in the control group", fontlist, tbl))
  
  tbl[8,1] <- (phi0Label1<-newlabel("phi0:", fontlist, tbl))
  tbl[8,2] <- (phi0Edit<-newedit("1", fontlist, tbl))
  tbl[8,3] <- (phi0Label2<-newlabel("dispersion for prognostic genes in control group", fontlist, tbl))
  
  samplesize <- gframe("Calculation Result", container=topgroup, horizontal=F)
  resultEdit<-gtext("Estimated sample size=", container=samplesize, expand=TRUE, font.attr=fontlist)
  enabled(resultEdit)<-FALSE
  
  description <- gframe("Description", container=topgroup, horizontal=F, expand=TRUE)
  dText<-gtext("", container=description, expand=TRUE, font.attr=fontlist)
  enabled(dText)<-FALSE
  
  doCalculation<-function(h, ...){
    svalue(resultEdit)<-"Estimated sample size=???"
    font(resultEdit)<-bluefontlist
    svalue(dText)<-""
    ret<-sample_size(as.numeric(svalue(mEdit)), 
                     as.numeric(svalue(m1Edit)), 
                     as.numeric(svalue(powerEdit)),
                     as.numeric(svalue(fEdit)),
                     as.numeric(svalue(wEdit)),
                     as.numeric(svalue(rhoEdit)),
                     as.numeric(svalue(mu0Edit)),
                     as.numeric(svalue(phi0Edit)))
    svalue(dText)<-paste0("We are planning a study to identify differential gene expression between two groups. Prior data indicate that ",
                          "the minimum average read counts among the prognostic genes in the control group is ", svalue(mu0Edit), 
                          ", the maximum dispersion is ", svalue(phi0Edit),
                          ", and the ratio of the geometric mean of normalization factors is ", svalue(wEdit),
                          ". Suppose that the total number genes for testing is ", svalue(mEdit),
                          " and the top ", svalue(m1Edit), 
                          " genes are prognostic. If the desired minimum fold changes is ", svalue(rhoEdit),
                          ", we will need to study ", ret,
                          " subjects to be able to reject the null hypothesis that the fold changes is 1 with probability (power) ", svalue(powerEdit),
                          " (", as.numeric(svalue(m1Edit)) * as.numeric(svalue(powerEdit)), "/", svalue(m1Edit),
                          ") using exact test. The FDR associated with this test of this null hypothesis is ", svalue(fEdit) ,".")
    svalue(resultEdit)<-paste0("Estimated sample size=", ret)
    font(resultEdit)<-redfontlist
  }
  
  button.group<-ggroup(container=topgroup)
  addSpring(button.group)
  calcbutton<-gbutton("Calculate", container=button.group, handler=doCalculation)
  copybutton<-gbutton("Copy description", container=button.group, handler = function(h, ...) writeClipboard(svalue(dText)) )
  closebutton<-gbutton("Close", container=button.group, handler = function(h,...) dispose(win))
  
  font(calcbutton)<-fontlist
  font(copybutton)<-fontlist
  font(closebutton)<-fontlist
  
  visible(win)<-TRUE
}

getvalues<-function(editFrom, name){
  v<-svalue(editFrom)
  result<-stringr::str_trim(unlist(strsplit(v, c(',',';'))))
  result<-as.double(result)
  if(any(is.na(result))){
    stop(paste("The values of ", name, " should be all numeric: ", v))
  }
  return(result)
}

##' rsse_batch_ui
##' 
##' A function to display user friendly interface to do estitamete the sample size for differential expression analysis of RNA-seq data at batch mode.
##' The result can be saved to file in csv format.
##' 
##' @export
##' @examples 
##' \dontrun{
##' rsse_batch_ui()
##' }
rsse_batch_ui<-function(){
  win <- gwindow("RNASeq Sample Size Estimation 0.98.3",width=1024, height=768, visible=FALSE)
  topgroup<-ggroup(container=win,horizontal=F)
  fontlist <- list(size=14)
  redfontlist <- list(size=14, color="red")
  bluefontlist <- list(size=14, color="blue")
  editwidth<-25
  
  parameters <- gframe("Parameters", container=topgroup, horizontal=F)
  font(parameters)<-fontlist
  tbl <- glayout(container = parameters)
  
  tbl[1,1] <- newlabel("m:", fontlist, tbl)
  tbl[1,2] <- (mEdit1<-newedit("10000", fontlist, tbl, width=editwidth))
  tbl[1,3] <- newlabel("total number of genes for testing", fontlist, tbl)
  
  tbl[2,1] <- newlabel("m1:", fontlist, tbl)
  tbl[2,2] <- (m1Edit1<-newedit("100", fontlist, tbl, width=editwidth))
  tbl[2,3] <- newlabel("expected number of prognostic genes", fontlist, tbl)
  
  tbl[3,1] <- newlabel("power:", fontlist, tbl)
  tbl[3,2] <- (powerEdit1<-newedit("0.8,0.9", fontlist, tbl, width=editwidth))
  tbl[3,3] <- newlabel("power to detect prognostic genes", fontlist, tbl)
  
  tbl[4,1] <- newlabel("f:", fontlist, tbl)
  tbl[4,2] <- (fEdit1<-newedit("0.1", fontlist, tbl, width=editwidth))
  tbl[4,3] <- newlabel("FDR level", fontlist, tbl)
  
  tbl[5,1] <- newlabel("w:", fontlist, tbl)
  tbl[5,2] <- (wEdit1<-newedit("1.0", fontlist, tbl, width=editwidth))
  tbl[5,3] <- newlabel("ratio of normalization factors between two groups", fontlist, tbl)
  
  tbl[6,1] <- newlabel("rho:", fontlist, tbl)
  tbl[6,2] <- (rhoEdit1<-newedit("2.0,3.0", fontlist, tbl, width=editwidth))
  tbl[6,3] <- newlabel("minimum fold changes for prognostic genes in control group", fontlist, tbl)
  
  tbl[7,1] <- newlabel("mu0:", fontlist, tbl)
  tbl[7,2] <- (mu0Edit1<-newedit("5,10", fontlist, tbl, width=editwidth))
  tbl[7,3] <- newlabel("average read counts for prognostic gene in control group", fontlist, tbl)
  
  tbl[8,1] <- newlabel("phi0:", fontlist, tbl)
  tbl[8,2] <- (phi0Edit1<-newedit("1", fontlist, tbl, width=editwidth))
  tbl[8,3] <- newlabel("dispersion for prognostic genes in control group", fontlist, tbl)
  
  samplesize <- gframe("Calculation Result", container=topgroup, horizontal=F, expand=TRUE)
  result<-as.matrix(data.frame(m=rep(0, 1),
                               m1=rep(0, 1),
                               power=rep(0, 1),
                               f=rep(0, 1),
                               w=rep(0, 1),
                               rho=rep(0, 1),
                               mu0=rep(0, 1),
                               phi0=rep(0, 1),
                               sampleSize=rep(NA, 1)))
  gt<-NA
  hasvalue<-FALSE
  
  doCalculation<-function(h, ...){
    ms<-getvalues(mEdit1, "m")
    m1s<-getvalues(m1Edit1, "m1")
    powers<-getvalues(powerEdit1, "power")
    fs<-getvalues(fEdit1, "f")
    ws<-getvalues(wEdit1, "w")
    rhos<-getvalues(rhoEdit1, "rho")
    mu0s<-getvalues(mu0Edit1, "mu0")
    phi0s<-getvalues(phi0Edit1, "phi0")
    
    totalcount<-length(ms) * length(m1s) * length(powers) * length(fs) * length(ws) * length(rhos) * length(mu0s) * length(phi0s)
    result<<-as.matrix(data.frame(m=rep(0, totalcount),
                                  m1=rep(0, totalcount),
                                  power=rep(0, totalcount),
                                  f=rep(0, totalcount),
                                  w=rep(0, totalcount),
                                  rho=rep(0, totalcount),
                                  mu0=rep(0, totalcount),
                                  phi0=rep(0, totalcount),
                                  sampleSize=rep(NA, totalcount)))
    
    index<-0
    for(m in ms){
      for(m1 in m1s){
        for(power in powers){
          for(f in fs){
            for(w in ws){
              for(rho in rhos){
                for(mu0 in mu0s){
                  for(phi0 in phi0s){
                    index = index+1
                    result[index,]<<-c(m, m1, power, f, w, rho, mu0, phi0, NA)
                  }
                }
              }
            }
          }
        }
      }
    }
    
    if(hasvalue){
      delete(samplesize, gt)
    }
    gt<<-gtable(result, container=samplesize, expand=TRUE)
    font(gt)<-fontlist
    hasvalue<<-TRUE
    
    enabled(calcbutton)<-FALSE
    enabled(copybutton)<-FALSE

    lapply(c(1:totalcount), function(x){
      result[x,9]<<-sample_size(as.numeric(result[x,1]),
                                as.numeric(result[x,2]),
                                as.numeric(result[x,3]),
                                as.numeric(result[x,4]),
                                as.numeric(result[x,5]),
                                as.numeric(result[x,6]),
                                as.numeric(result[x,7]),
                                as.numeric(result[x,8]))
      gt[x,9]<-as.numeric(result[x,9])
    })
    
    enabled(calcbutton)<-TRUE
    enabled(copybutton)<-TRUE
  }
  
  button.group<-ggroup(container=topgroup)
  addSpring(button.group)
  calcbutton<-gbutton("Calculate", container=button.group, handler=doCalculation)
  copybutton<-gbutton("Save to file", container=button.group, handler = function(h, ...) {
    fname <- gfile(gettext("Filename to save to"), type="save")
    if(nchar(fname)) {
      if(!file.exists(fname) || gconfirm(gettext("Overwrite file?"))){
        write.csv(result, fname)
      }
    }
  } )
  closebutton<-gbutton("Close", container=button.group, handler = function(h,...) dispose(win))
  
  font(calcbutton)<-fontlist
  font(copybutton)<-fontlist
  font(closebutton)<-fontlist
  
  visible(win)<-TRUE
}
