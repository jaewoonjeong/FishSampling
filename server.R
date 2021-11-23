function(input, output,session){
  observeEvent(input$btn1.1,{shinyalert(title = "Probability of corresponding to the region of density plot between the two dashed green vertical lines", type = "info")})
  observeEvent(input$btn1.2,{shinyalert(title = "The products of 'Number of sampled pen' and 'Number of sampled fish per pen' may be different to 'Number of total sampled fish' due to a rounding issue. The number of total sampled fish that actually used in simulations was the product, not the number of total fish to sample you supposed.", type = "info")})
  observeEvent(input$btn1.3,{shinyalert(title = "The black solid and green dashed lines represent abundance and higher and lower target limits, respectively", type = "info")})
  observeEvent(input$btn2.1,{shinyalert(title = "The probability corresponds to the right region to the red dashed line.", type = "info")})
  observeEvent(input$btn2.2,{shinyalert(title = "The products of 'Number of sampled pen' and 'Number of sampled fish per pen' may be different to 'Number of total sampled fish' due to a rounding issue. The number of total sampled fish that actually used in simulations was the product, not the number of total fish to sample you supposed.", type = "info")})
  observeEvent(input$btn2.3,{shinyalert(title = "The black solid line and red dashed line represent abundance and lice limit, respectively", type = "info")})
  observeEvent(input$btn3.1,{shinyalert(title = "The probability corresponds to the region of density plot to LEFT of the RED dashed vertical line, when Abundance2 is LOWER than Abundance1. The probability corresponds to the region of density plot to RIGHT of the RED dashed vertical line, when Abundance2 is HIGHER than Abundance1.", type = "info")})
  observeEvent(input$btn3.2,{shinyalert(title = "The products of 'Number of sampled pen' and 'Number of sampled fish per pen' may be different to 'Number of total sampled fish' due to a rounding issue. The number of total sampled fish that actually used in simulations was the product, not the number of total fish to sample you supposed.", type = "info")})
  observeEvent(input$btn3.3,{shinyalert(title = "Vertical lines represent the two mean abundance values", type = "info")})
  
  extract <- function(text) {text <- gsub(" ", "", text); split <- strsplit(text, ",", fixed = FALSE)[[1]]; as.numeric(split)}
  
  kappa=2.19 # dispersion parameter of negative binomial distribution
  fish=2000 # number of fish per pen
  
  #######################  SCENARIO 1
  mydata1<-reactive({if(input$clustering1=='no'){TotalSampleSize=as.numeric(input$TSS1)
  withProgress(message = 'Progress time',value = 0,{N<-10
  for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
    sam.mean=matrix(NA,length(TotalSampleSize),input$iteration1)
    for(i in 1:input$iteration1){for(s in 1:length(TotalSampleSize)){sam.mean[s,i]=mean(sample(rnbinom(fish*10,size=kappa,mu=input$ab1),TotalSampleSize[s]))}}}})
  results1<-list(s=sam.mean)} else {choice1=as.numeric(input$choice1)
  gap=seq(0,input$ab1,length.out=200)
  rho=numeric()
  g=1;  rho1=0
  while(rho1<input$icc1 & g<length(gap)){
    gg=seq(0,gap[g],length.out=round(input$pen1/2))
    if(input$pen1%%2==0){pen.AB1=c(rev(input$ab1-gg),input$ab1+gg)} else {pen.AB1=c(rev(input$ab1-gg),input$ab1+gg[-1])}
    PP=matrix(NA,fish,input$pen1)
    for(p in 1:input$pen1){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB1[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    rho1=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
    g=g+1}
  if(rho1==input$icc1){rho=rho1; pen.AB=pen.AB1} else{g=1;    rho2=0
  while(rho2<input$icc1 & g<length(gap)){
    if(input$pen1%%2==0){
      pen.AB2=c(rep(input$ab1-gap[g],round(input$pen1/2)),rep(input$ab1+gap[g],round(input$pen1/2)))} else {
        pen.AB2=c(rep(input$ab1-gap[g],round(input$pen1/2)),rep(input$ab1+gap[g],round(input$pen1/2))[-1])}
    PP=matrix(NA,fish,input$pen1)
    for(p in 1:input$pen1){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB2[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    rho2=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
    g=g+1}  
  rho=rho2; pen.AB=pen.AB2}
  
  TotalSampleSize=as.numeric(input$TSS1)
  SampleSizePerPen=matrix(NA,length(TotalSampleSize),length(choice1))
  for(s in 1:length(TotalSampleSize)){for(c in 1:length(choice1)){SampleSizePerPen[s,c]=round(TotalSampleSize[s]/choice1[c])}}
  withProgress(message = 'Progress time',value = 0,{N<-10
  for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
    sam.mean=matrix(NA,length(choice1)*length(TotalSampleSize)*input$iteration1); dim(sam.mean)=c(length(choice1),length(TotalSampleSize),input$iteration1)
    for(i in 1:input$iteration1){sam.pen=matrix(NA,input$pen1*length(choice1)*length(TotalSampleSize)); dim(sam.pen)=c(input$pen1,length(choice1),length(TotalSampleSize))
    for(s in 1:length(TotalSampleSize)){for(c in 1:length(choice1)){pp=sort(sample(1:input$pen1,choice1[c]))
    for(p in pp){sam.pen[p,c,s]=mean(sample(rnbinom(fish,size=kappa,mu=pen.AB[p]),SampleSizePerPen[s,c]))}}}
    sam.mean[,,i]=colMeans(sam.pen,na.rm=T)}}})
  results1<-list(b=sam.mean, c=SampleSizePerPen, d=pen.AB)}})
  
  #######################  OUTPUT 1B: 
  
  plo.1B.2<-eventReactive(input$goButton1, {results1<-mydata1()
  if(input$clustering1=='yes'){
    plot(x=1:length(results1$d),y=results1$d[sample(1:length(results1$d),length(results1$d))],las=1,xlab="Pen number",ylab="Abundance",main=c("Example abundance of each pen","with the supposed ICC"))
    abline(h=input$ab1,col='gray',lty=2)}})
  output$Plot.1B.2<-renderPlot({plo.1B.2()})
  
  plo.1B.3<-eventReactive(input$goButton1, {results1<-mydata1()
  if(input$clustering1=='no'){md=melt(results1$s) 
  md[,1]=as.numeric(input$TSS1)[md[,1]]; md[,1]=as.factor(md[,1]); md[,2]=as.factor(md[,2])
  colnames(md)=c("NumberOfTotalSampledFish","Iteration","Abundance")
  ggplot(md,aes(x=Abundance,y=NumberOfTotalSampledFish))+geom_density_ridges()+geom_vline(xintercept=input$ab1)+
    geom_vline(xintercept=c(input$u.m,input$l.m),colour=c("green"), linetype = "dashed")+ylab(ifelse(length(as.numeric(input$TSS1))==1,"Density","Total number of sampled fish"))+
    theme(text=element_text(size=18))+theme(axis.text.x=element_text(face="bold",size=18),axis.text.y = element_text(face="bold",size=18))
  } else {
    md=melt(results1$b)
    md[,1]=as.numeric(input$choice1)[md[,1]]; md[,1]=as.factor(md[,1]); 
    md[,2]=as.numeric(input$TSS1)[md[,2]]; md[,2]=as.factor(md[,2]); 
    md[,3]=as.factor(md[,3])
    colnames(md)=c("NumberOfSampledPens","NumberOfTotalSampledFish","Iteration","Abundance")
    ggplot(md, aes(x=Abundance,y=NumberOfSampledPens))+geom_density_ridges()+facet_wrap(~NumberOfTotalSampledFish)+theme(legend.position="bottom",legend.title=element_text("Number of sampled pens",size=18,face="bold"))+ggtitle("Total number of sampled fish")+ylab(ifelse(length(as.numeric(input$choice1))==1,"Density","Number of sampled pens"))+theme(plot.title=element_text(size=18),text=element_text(size=18))+theme(axis.text.x=element_text(face="bold",size=18),axis.text.y=element_text(face="bold",size=18))+geom_vline(xintercept=mean(extract(input$ab1)))+geom_vline(xintercept=c(input$u.m,input$l.m),colour=c("green"),linetype="dashed")}})
  output$Plot1<-renderPlot({plo.1B.3()})
  
  dataset.1<-eventReactive(input$goButton1,{results1<-mydata1()
  if(input$clustering1=='no'){
    CI=matrix(NA,length(as.numeric(input$TSS1)),1)
    for(s in 1:length(as.numeric(input$TSS1))){
      CI[s]=paste(format(round(mean(results1$s[s,]<input$u.m & results1$s[s,]>input$l.m),2),nsmall=2),"(",sum(results1$s[s,]<input$u.m & results1$s[s,]>input$l.m),"/",input$iteration1,")")}
    CI=matrix(CI[nrow(CI):1,],ncol=1)
    colnames(CI)=" "; rownames(CI)=rev(paste("Total Sample Size",as.numeric(input$TSS1)))   
    CI} else {
      CI=matrix(NA,length(as.numeric(input$TSS1)),length(as.numeric(input$choice1)))
      for(s in 1:length(as.numeric(input$TSS1))){for(t in 1:length(as.numeric(input$choice1))){
        CI[s,t]=paste(format(round(mean(results1$b[t,s,]<input$u.m & results1$b[t,s,]>input$l.m),2),nsmall=2),"(",sum(results1$b[t,s,]<input$u.m & results1$b[t,s,]>input$l.m),"/",input$iteration1,")")
      }}
      colnames(CI)=paste("Number of sampled pens",as.numeric(input$choice1))
      rownames(CI)=paste("Total number of sampled fish",as.numeric(input$TSS1))
      if(nrow(CI)==1 | ncol(CI)==1){CI=t(CI)} else {t(CI[,ncol(CI):1])}}
  })
  output$Table.1.1 <- renderTable({dataset.1()},rownames = TRUE)
  tt.1<-eventReactive(input$goButton1, {results1<-mydata1()
  SSperPen=format(results1$c)
  rownames(SSperPen)=paste("Total number of sampled fish",as.numeric(input$TSS1))
  colnames(SSperPen)=paste("Number of sampled pens",as.numeric(input$choice1))
  if(nrow(SSperPen)==1 | ncol(SSperPen)==1){SSperPen=t(SSperPen)} else {t(SSperPen[,ncol(SSperPen):1])}})
  output$Table.1.2<-renderTable({tt.1()},rownames = TRUE)
  
  ################################################## Scenario 2
  mydata2<-reactive({
    if(input$clustering2=='no'){pen.AB=c(input$ab1.2B,input$ab2.2B)
    TotalSampleSize=input$TSS2
    withProgress(message = 'Progress time',value = 0,{N<-10
    for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
      sam.mean=matrix(NA,length(TotalSampleSize)*input$iteration2*length(pen.AB)); dim(sam.mean)=c(length(TotalSampleSize),input$iteration2,length(pen.AB))
      for(s in 1:length(TotalSampleSize)){for(i in 1:input$iteration2){for(a in 1:length(pen.AB)){sam.mean[s,i,a]=mean(sample(rnbinom(fish*10,size=kappa,mu=pen.AB[a]),TotalSampleSize[s]))}}}}})
    results2<-list(h=sam.mean,z=sam.mean)} else {choice2=as.numeric(input$choice2)
    TotalSampleSize=as.numeric(input$TSS2) 
    SampleSizePerPen=matrix(NA,length(TotalSampleSize),length(choice2))
    for(s in 1:length(TotalSampleSize)){for(c in 1:length(choice2)){
      SampleSizePerPen[s,c]=round(TotalSampleSize[s]/choice2[c])
    }}
    gap=seq(0,input$ab1.2B,length.out=200)
    rho=numeric()
    g=1;  rho1=0
    while(rho1<input$icc1.2B & g<length(gap)){
      gg=seq(0,gap[g],length.out=round(input$pen2/2))
      if(input$pen2%%2==0){pen.AB1=c(rev(input$ab1.2B-gg),input$ab1.2B+gg)} else {
        pen.AB1=c(rev(input$ab1.2B-gg),input$ab1.2B+gg[-1])}
      PP=matrix(NA,fish,input$pen2)
      for(p in 1:input$pen2){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB1[p])}
      mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
      rho1=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
      g=g+1}
    if(rho1==input$icc1.2B){rho=rho1; pen.AB=pen.AB1} else{g=1;    rho2=0
    while(rho2<input$icc1.2B & g<length(gap)){
      if(input$pen2%%2==0){pen.AB2=c(rep(input$ab1.2B-gap[g],round(input$pen2/2)),rep(input$ab1.2B+gap[g],round(input$pen2/2)))} else {
        pen.AB2=c(rep(input$ab1.2B-gap[g],round(input$pen2/2)),rep(input$ab1.2B+gap[g],round(input$pen2/2))[-1])}
      PP=matrix(NA,fish,input$pen2)
      for(p in 1:input$pen2){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB2[p])}
      mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
      rho2=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
      g=g+1}  
    rho=rho2; pen.AB=pen.AB2}
    pen.AB11=pen.AB
    
    gap=seq(0,input$ab2.2B,length.out=200)
    rho=numeric()
    g=1;  rho1=0
    while(rho1<input$icc2.2B & g<length(gap)){
      gg=seq(0,gap[g],length.out=round(input$pen2/2))
      if(input$pen2%%2==0){pen.AB1=c(rev(input$ab2.2B-gg),input$ab2.2B+gg)} else {
        pen.AB1=c(rev(input$ab2.2B-gg),input$ab2.2B+gg[-1])}
      PP=matrix(NA,fish,input$pen2)
      for(p in 1:input$pen2){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB1[p])}
      mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
      rho1=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
      g=g+1}
    if(rho1==input$icc2.2B){rho=rho1; pen.AB=pen.AB1} else{g=1;    rho2=0
    while(rho2<input$icc2.2B & g<length(gap)){
      if(input$pen2%%2==0){pen.AB2=c(rep(input$ab2.2B-gap[g],round(input$pen2/2)),rep(input$ab2.2B+gap[g],round(input$pen2/2)))} else {
        pen.AB2=c(rep(input$ab2.2B-gap[g],round(input$pen2/2)),rep(input$ab2.2B+gap[g],round(input$pen2/2))[-1])}
      PP=matrix(NA,fish,input$pen2)
      for(p in 1:input$pen2){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB2[p])}
      mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
      rho2=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
      g=g+1}  
    rho=rho2; pen.AB=pen.AB2}
    pen.AB22=pen.AB
    
    pen.AB=matrix(NA,input$pen2,2); pen.AB[,1]=pen.AB11; pen.AB[,2]=pen.AB22
    withProgress(message = 'Progress time',value = 0,{N<-10
    for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
      sam.mean=matrix(NA,length(choice2),length(TotalSampleSize)*input$iteration2*2); dim(sam.mean)=c(length(choice2),length(TotalSampleSize),input$iteration2,2)
      for(a in 1:2){for(i in 1:input$iteration2){sam.pen=matrix(NA,input$pen2*length(choice2)*length(TotalSampleSize)); dim(sam.pen)=c(input$pen2,length(choice2),length(TotalSampleSize))
      for(s in 1:length(TotalSampleSize)){for(c in 1:length(choice2)){pp=sort(sample(1:input$pen2,choice2[c]))
      for(p in pp){sam.pen[p,c,s]=mean(sample(rnbinom(fish,size=kappa,mu=pen.AB[p,a]),SampleSizePerPen[s,c]))}}}
      sam.mean[,,i,a]=colMeans(sam.pen,na.rm=T)}}}}) # withProgress
    results2<- list(a=sam.mean,b=SampleSizePerPen,c=pen.AB11,d=pen.AB22)}})
  
  ##################### OUTPUT
  plo.2B.1<-eventReactive(input$goButton2, {results2<-mydata2()
  if(input$clustering2=='yes'){
    plot(x=1:length(results2$c),y=results2$c[sample(1:length(results2$c),length(results2$c))],las=1,xlab="Pen number",ylab="Abundance",main=c("Example abundance of each pen","with the supposed true abundance",sub="True Abundance 1"))
    abline(h=input$ab1.2B,col='gray',lty=2)
  }})
  output$Plot.2B.1<-renderPlot({plo.2B.1()})
  
  plo.2B.2<-eventReactive(input$goButton2, {results2<-mydata2()
  if(input$clustering2=='yes'){
    plot(x=1:length(results2$d),y=results2$d[sample(1:length(results2$d),length(results2$d))],las=1,xlab="Pen number",ylab="Abundance",main=c("Example abundance of each pen","with the supposed true abundance",sub="True Abundance 2"))
    abline(h=input$ab2.2B,col='gray',lty=2)
  }})
  output$Plot.2B.2<-renderPlot({plo.2B.2()})
  
  plo.2B.3<-eventReactive(input$goButton2, {results2<-mydata2()
  if(input$clustering2=='no'){
    dt=melt(results2$h); dt[,1]=as.numeric(input$TSS2)[dt[,1]]
    dt[,3]=ifelse(dt[,3]==1,"Abundance1","Abundance2")
    colnames(dt)=c("NumberOfTotalSampledFish","Iteration","ABUNDANCE","Abundance")
    ggplot(dt, aes(x=Abundance,colour=ABUNDANCE))+scale_color_manual(values=c('red','blue'))+geom_density(aes(group=ABUNDANCE),alpha=0.25)+facet_wrap(~NumberOfTotalSampledFish)+ggtitle("Total number of sampled fish")+theme(plot.title = element_text(size=18))+geom_vline(xintercept =c(input$ab1.2B,input$ab2.2B), linetype = "dashed", colour = rep(c(2,4),length(table(dt[,1]))))+ylab("Density")+theme(legend.position="bottom",text=element_text(size=18),axis.text.x=element_text(face="bold",size=18),axis.text.y=element_text(face="bold",size=18))+guides(col=guide_legend(""))
  } else {
    dt=melt(results2$a); dt[,1]=as.numeric(input$choice2)[dt[,1]]; dt[,2]=as.numeric(input$TSS2)[dt[,2]]; dt[,4]=ifelse(dt[,4]==1,"Abundance1","Abundance2")
    colnames(dt)=c("NumberOfSampledPens","NumberOfTotalSampledFish","Iteration","ABUNDANCE","Abundance")
    ggplot(dt,aes(x=Abundance,colour=ABUNDANCE))+scale_color_manual(values=c('red','blue'))+geom_density(aes(group=ABUNDANCE),alpha=0.25)+facet_grid(NumberOfSampledPens~NumberOfTotalSampledFish)+ggtitle("Total number of sampled fish")+theme(plot.title = element_text(size=18))+geom_vline(xintercept=c(input$ab1.2B,input$ab2.2B),linetype="dashed",colour=rep(c(2,4),length(table(dt[,2]))*length(table(dt[,1]))))+ylab("Density")+theme(legend.position="bottom",text=element_text(size=18),axis.text.x=element_text(face="bold",size=18),axis.text.y=element_text(face="bold",size=18))+guides(col=guide_legend(""))}})
  output$Plot2<-renderPlot({plo.2B.3()})
  
  dataset.2<-eventReactive(input$goButton2,{results2<-mydata2()
  if(input$clustering2=='no'){abc=matrix(NA,length(as.numeric(input$TSS2)),1)
  for(s in 1:length(as.numeric(input$TSS2))){abc[s]=ifelse(input$ab1.2B<input$ab2.2B, paste(format(round(mean(results2$h[s,,1]<results2$h[s,,2]),2),nsmall=2),"(",sum(results2$h[s,,1]<results2$h[s,,2]),"/",input$iteration2,")"),paste(round(mean(results2$h[s,,1]>results2$h[s,,2]),2),"(",sum(results2$h[s,,1]>results2$h[s,,2]),"/",input$iteration2,")"))}
  colnames(abc)=c(" "); rownames(abc)=paste("Total Sample Size",as.numeric(input$TSS2))
  abc} else {abc=matrix(NA,length(as.numeric(input$TSS2)),length(as.numeric(input$choice2)))
  for(s in 1:length(as.numeric(input$TSS2))){for(c in 1:length(as.numeric(input$choice2))){
    abc[s,c]=ifelse(input$ab1.2B<input$ab2.2B, paste(format(round(mean(results2$a[c,s,,1]<results2$a[c,s,,2]),2),nsmall=2),"(",sum(results2$a[c,s,,1]<results2$a[c,s,,2]),"/",input$iteration2,")"), paste(round(mean(results2$a[c,s,,1]>results2$a[c,s,,2]),2),"(",sum(results2$a[c,s,,1]>results2$a[c,s,,2]),"/",input$iteration2,")"))}}
  colnames(abc)=paste("Number of sampled pens",as.numeric(input$choice2)); rownames(abc)=paste("Total number of sampled fish",as.numeric(input$TSS2))
  t(abc)}})
  output$Table.2.1<-renderTable({dataset.2()},rownames = TRUE)
  
  tt.2<-eventReactive(input$goButton2, {results2<-mydata2()
  SSperPen=format(results2$b)
  colnames(SSperPen)=paste("Number of sampled pens",as.numeric(input$choice2))
  rownames(SSperPen)=paste("Total number of sampled fish",as.numeric(input$TSS2))
  if(nrow(SSperPen)==1 | ncol(SSperPen)==1){SSperPen=t(SSperPen)} else {t(SSperPen[,ncol(SSperPen):1])}
  })
  output$Table.2.2<-renderTable({tt.2()},rownames = TRUE)
  
  ############################## Scenario 3
  mydata3<-reactive({if(input$clustering3=='no'){TotalSampleSize=as.numeric(input$TSS3)
  withProgress(message = 'Progress time',value = 0,{N<-10
  for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
    sam.mean=matrix(NA,length(TotalSampleSize),input$iteration3)
    for(i in 1:input$iteration3){
      for(s in 1:length(TotalSampleSize)){
        sam.mean[s,i]=mean(sample(rnbinom(fish*10,size=kappa,mu=input$ab3),TotalSampleSize[s]))}}}})
  results3<-list(d=sam.mean)
  } else {
    choice3=as.numeric(input$choice3)
    gap=seq(0,input$ab3,length.out=200)
    rho=numeric()
    g=1;  rho1=0
    while(rho1<input$icc3 & g<length(gap)){
      gg=seq(0,gap[g],length.out=round(input$pen3/2))
      if(input$pen3%%2==0){pen.AB1=c(rev(input$ab3-gg),input$ab3+gg)} else {pen.AB1=c(rev(input$ab3-gg),input$ab3+gg[-1])}
      PP=matrix(NA,fish,input$pen3)
      for(p in 1:input$pen3){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB1[p])}
      mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
      rho1=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
      g=g+1}
    if(rho1==input$icc3){rho=rho1; pen.AB=pen.AB1} else{g=1;    rho2=0
    while(rho2<input$icc3 & g<length(gap)){
      if(input$pen3%%2==0){pen.AB2=c(rep(input$ab3-gap[g],round(input$pen3/2)),rep(input$ab3+gap[g],round(input$pen3/2)))} else {
        pen.AB2=c(rep(input$ab3-gap[g],round(input$pen3/2)),rep(input$ab3+gap[g],round(input$pen3/2))[-1])}
      PP=matrix(NA,fish,input$pen3)
      for(p in 1:input$pen3){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB2[p])}
      mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
      rho2=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
      g=g+1}  
    rho=rho2; pen.AB=pen.AB2}
    
    TotalSampleSize=as.numeric(input$TSS3)
    SampleSizePerPen=matrix(NA,length(TotalSampleSize),length(choice3))
    for(s in 1:length(TotalSampleSize)){for(c in 1:length(choice3)){SampleSizePerPen[s,c]=round(TotalSampleSize[s]/choice3[c])}}
    withProgress(message = 'Progress time',value = 0,{N<-10
    for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
      sam.mean=matrix(NA,length(TotalSampleSize)*length(choice3)*input$iteration3); dim(sam.mean)=c(length(TotalSampleSize),length(choice3),input$iteration3)
      for(s in 1:length(TotalSampleSize)){for(i in 1:input$iteration3){sam.pen=matrix(NA,input$pen3,length(choice3))
      for(c in 1:length(choice3)){pp=sort(sample(1:input$pen3,choice3[c]))
      for(p in pp){sam.pen[p,c]=mean(sample(rnbinom(fish,size=kappa,mu=pen.AB[p]),SampleSizePerPen[s,c]))}}
      sam.mean[s,,i]=colMeans(sam.pen,na.rm=T)}}}})
    results3<-list(a=sam.mean,b=SampleSizePerPen,c=pen.AB)}})
  
  #######################  OUTPUT 3B: 
  
  plo.3B.2<-eventReactive(input$goButton3, {results3<-mydata3()
  if(input$clustering3=='yes'){
    plot(x=1:length(results3$c),y=results3$c[sample(1:length(results3$c),length(results3$c))],las=1,xlab="Pen number",ylab="Abundance",main=c("Example abundance of each pen","with the supposed ICC"))
    abline(h=input$ab3,col='gray',lty=2)}})
  output$Plot.3B.2<-renderPlot({plo.3B.2()})
  
  plo.3B.1<-eventReactive(input$goButton3, {results3<-mydata3()
  if(input$clustering3=='no'){dt=melt(results3$d)
  dt[,1]=as.numeric(input$TSS3)[dt[,1]]; dt[,1]=as.factor(dt[,1]); dt[,2]=as.factor(dt[,2])
  colnames(dt)=c("NumberOfTotalSampledFish","Iteration","Abundance")
  ggplot(dt, aes(x=Abundance,y=NumberOfTotalSampledFish))+geom_density_ridges()+geom_vline(xintercept=input$ab3)+geom_vline(xintercept=input$th3,linetype="dashed",colour="red")+ylab(ifelse(length(as.numeric(input$TSS3))==1,"Density","Total number of sampled fish"))+theme(text=element_text(size=18))+theme(axis.text.x=element_text(face="bold",size=18),axis.text.y = element_text(face="bold",size=18))
  } else {dt=melt(results3$a)
  dt[,1]=as.numeric(input$TSS3)[dt[,1]]; dt[,2]=as.numeric(input$choice3)[dt[,2]]; dt[,2]=as.factor(dt[,2])
  colnames(dt)=c("NumberOfTotalSampledFish","NumberOfSampledPens","Iteration","Abundance")
  ggplot(dt, aes(x=Abundance, y=NumberOfSampledPens)) +geom_density_ridges()+facet_wrap(~NumberOfTotalSampledFish)+ggtitle("Total number of sampled fish")+theme(plot.title = element_text(size=18))+
    geom_vline(xintercept=mean(extract(input$ab3)))+geom_vline(xintercept=input$th3,linetype="dashed",colour="red")+
    ylab(ifelse(length(as.numeric(input$choice3))==1,"Density","Number of sampled pens"))+
    theme(legend.position="bottom",legend.title=element_text("Number of sampled pens",size=18,face="bold"))+theme(text=element_text(size=18))+
    theme(axis.text.x=element_text(face="bold",size=18),axis.text.y=element_text(face="bold",size=18))
  }})
  output$Plot3 <- renderPlot({plo.3B.1()})
  
  dataset.3<-eventReactive(input$goButton3,{results3<-mydata3()
  if(input$clustering3=='no'){cutline=matrix(NA,length(as.numeric(input$TSS3)),1)
  for(s in 1:length(as.numeric(input$TSS3))){
    cutline[s]=paste(format(1-round(mean(results3$d[s,]<input$th3),2),nsmall=2),"(",sum(results3$d[s,]>=input$th3),"/",input$iteration3,")")}
  cutline=matrix(cutline[nrow(cutline):1,],ncol=1)
  colnames(cutline)=" "; rownames(cutline)=rev(paste("Total number of sampled fish",as.numeric(input$TSS3)))
  cutline
  } else {
    cutline=matrix(NA,length(as.numeric(input$TSS3)),length(as.numeric(input$choice3)))
    for(c in 1:length(as.numeric(input$choice3))){
      for(s in 1:length(as.numeric(input$TSS3))){
        cutline[s,c]=paste(format(1-round(mean(results3$a[s,c,]<input$th3),2),nsmall=2),"(",sum(results3$a[s,c,]>=input$th3),"/",input$iteration3,")")}}
    colnames(cutline)=paste("Number of sampled pens",as.numeric(input$choice3)); rownames(cutline)=paste("Total number of sampled fish",as.numeric(input$TSS3))
    format(cutline,justify=c("right"))
    if(nrow(cutline)==1 | ncol(cutline)==1){cutline=t(cutline)} else {t(cutline[,ncol(cutline):1])}
    }
  })
  output$Table.3.1 <- renderTable({dataset.3()},rownames = TRUE)
  
  tt.3<-eventReactive(input$goButton3,{results3<-mydata3()
  SSperPen=format(results3$b)
  colnames(SSperPen)=paste("Number of sampled pens",as.numeric(input$choice3));
  rownames(SSperPen)=paste("Total number of sampled fish",as.numeric(input$TSS3))
  if(nrow(SSperPen)==1 | ncol(SSperPen)==1){SSperPen=t(SSperPen)} else {t(SSperPen[,ncol(SSperPen):1])}
  })
  output$Table.3.2<-renderTable({tt.3()},rownames = TRUE)
}
