server <-function(input, output,session){
  extract <- function(text) {text <- gsub(" ", "", text); split <- strsplit(text, ",", fixed = FALSE)[[1]]; as.numeric(split)}
  
  kappa=2.19 # dispersion parameter of negative binomial distribution
  fish=2000 # number of fish per pen
  
  #######################  SCENARIO 1B: mydata1B
  mydata.1B<-reactive({if(input$clustering1B=='no'){TotalSampleSize=as.numeric(input$TSS.1B)
  withProgress(message = 'Progress time',value = 0,{N<-10
  for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
    sam.mean=matrix(NA,length(TotalSampleSize),input$iteration.1B)
    for(i in 1:input$iteration.1B){for(s in 1:length(TotalSampleSize)){sam.mean[s,i]=mean(sample(rnbinom(fish*10,size=kappa,mu=input$ab.1B),TotalSampleSize[s]))}}}})
  results.1B<-list(s=sam.mean)} else {choice.1B=as.numeric(input$choice.1B)
  gap=seq(0,input$ab.1B,length.out=200)
  rho=numeric()
  g=1;  rho1=0
  while(rho1<input$icc.1B & g<length(gap)){
    gg=seq(0,gap[g],length.out=round(input$pen.1B/2))
    if(input$pen.1B%%2==0){pen.AB1=c(rev(input$ab.1B-gg),input$ab.1B+gg)} else {pen.AB1=c(rev(input$ab.1B-gg),input$ab.1B+gg[-1])}
    PP=matrix(NA,fish,input$pen.1B)
    for(p in 1:input$pen.1B){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB1[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    rho1=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
    g=g+1}
  if(rho1==input$icc.1B){rho=rho1; pen.AB=pen.AB1} else{g=1;    rho2=0
  while(rho2<input$icc.1B & g<length(gap)){
    if(input$pen.1B%%2==0){pen.AB2=c(rep(input$ab.1B-gap[g],round(input$pen.1B/2)),rep(input$ab.1B+gap[g],round(input$pen.1B/2)))} else {
      pen.AB2=c(rep(input$ab.1B-gap[g],round(input$pen.1B/2)),rep(input$ab.1B+gap[g],round(input$pen.1B/2))[-1])}
    PP=matrix(NA,fish,input$pen.1B)
    for(p in 1:input$pen.1B){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB2[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    rho2=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
    g=g+1}  
  rho=rho2; pen.AB=pen.AB2}
  
  TotalSampleSize=as.numeric(input$TSS.1B)
  SampleSizePerPen=matrix(NA,length(TotalSampleSize),length(choice.1B))
  for(s in 1:length(TotalSampleSize)){for(c in 1:length(choice.1B)){SampleSizePerPen[s,c]=round(TotalSampleSize[s]/choice.1B[c])}}
  withProgress(message = 'Progress time',value = 0,{N<-10
  for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
    sam.mean=matrix(NA,length(choice.1B)*length(TotalSampleSize)*input$iteration.1B); dim(sam.mean)=c(length(choice.1B),length(TotalSampleSize),input$iteration.1B)
    for(i in 1:input$iteration.1B){sam.pen=matrix(NA,input$pen.1B*length(choice.1B)*length(TotalSampleSize)); dim(sam.pen)=c(input$pen.1B,length(choice.1B),length(TotalSampleSize))
    for(s in 1:length(TotalSampleSize)){for(c in 1:length(choice.1B)){pp=sort(sample(1:input$pen.1B,choice.1B[c]))
    for(p in pp){sam.pen[p,c,s]=mean(sample(rnbinom(fish,size=kappa,mu=pen.AB[p]),SampleSizePerPen[s,c]))}}}
    sam.mean[,,i]=colMeans(sam.pen,na.rm=T)}}})
  results.1B<-list(b=sam.mean, c=SampleSizePerPen, d=pen.AB)}})
  
  #######################  OUTPUT 1B: 
  
  plo.1B.2<-eventReactive(input$goButton.1B, {results.1B<-mydata.1B()
  if(input$clustering1B=='yes'){
    plot(x=1:length(results.1B$d),y=results.1B$d[sample(1:length(results.1B$d),length(results.1B$d))],las=1,xlab="Pen number",ylab="Abundance",main=c("Example abundance of each pen","with the supposed ICC"))
    abline(h=input$ab.1B,col='gray',lty=2)}})
  output$Plot.1B.2<-renderPlot({plo.1B.2()})
  
  plo.1B.3<-eventReactive(input$goButton.1B, {results.1B<-mydata.1B()
  if(input$clustering1B=='no'){md=melt(results.1B$s) 
  md[,1]=as.numeric(input$TSS.1B)[md[,1]]; md[,1]=as.factor(md[,1]); md[,2]=as.factor(md[,2])
  colnames(md)=c("NumberOfTotalSampledFish","Iteration","Abundance")
  ggplot(md,aes(x=Abundance,y=NumberOfTotalSampledFish))+geom_density_ridges()+geom_vline(xintercept=input$ab.1B)+
    geom_vline(xintercept=c(input$u.m,input$l.m),colour=c("green"), linetype = "dashed")+ylab(ifelse(length(as.numeric(input$TSS.1B))==1,"Density","Total number of sampled fish"))+
    theme(text=element_text(size=18))+theme(axis.text.x=element_text(face="bold",size=18),axis.text.y = element_text(face="bold",size=18))
  } else {
    md=melt(results.1B$b)
    md[,1]=as.numeric(input$choice.1B)[md[,1]]; md[,1]=as.factor(md[,1]); 
    md[,2]=as.numeric(input$TSS.1B)[md[,2]]; md[,2]=as.factor(md[,2]); 
    md[,3]=as.factor(md[,3])
    colnames(md)=c("NumberOfSampledPens","NumberOfTotalSampledFish","Iteration","Abundance")
    ggplot(md, aes(x=Abundance,y=NumberOfSampledPens))+geom_density_ridges()+facet_wrap(~NumberOfTotalSampledFish)+theme(legend.position="bottom",legend.title=element_text("Number of sampled pens",size=18,face="bold"))+ggtitle("Total number of sampled fish")+ylab(ifelse(length(as.numeric(input$choice.1B))==1,"Density","Number of sampled pens"))+theme(plot.title=element_text(size=18),text=element_text(size=18))+theme(axis.text.x=element_text(face="bold",size=18),axis.text.y=element_text(face="bold",size=18))+geom_vline(xintercept=mean(extract(input$ab.1B)))+geom_vline(xintercept=c(input$u.m,input$l.m),colour=c("green"),linetype="dashed")}})
  output$Plot.1B.3<-renderPlot({plo.1B.3()})
  
  dataset.1<-eventReactive(input$goButton.1B,{results.1B<-mydata.1B()
  if(input$clustering1B=='no'){
    CI=matrix(NA,length(as.numeric(input$TSS.1B)),1)
    for(s in 1:length(as.numeric(input$TSS.1B))){CI[s]=paste(round(mean(results.1B$s[s,]<input$u.m & results.1B$s[s,]>input$l.m),2),"(",sum(results.1B$s[s,]<input$u.m & results.1B$s[s,]>input$l.m),"/",input$iteration.1B,")")}
    CI=matrix(CI[nrow(CI):1,],ncol=1)
    colnames(CI)="Probability"; rownames(CI)=rev(paste("Total Sample Size",as.numeric(input$TSS.1B)))   
    CI} else {CI=matrix(NA,length(as.numeric(input$TSS.1B)),length(as.numeric(input$choice.1B)))
    for(s in 1:length(as.numeric(input$TSS.1B))){for(t in 1:length(as.numeric(input$choice.1B))){
      CI[s,t]=paste(round(mean(results.1B$b[t,s,]<input$u.m & results.1B$b[t,s,]>input$l.m),2),"(",sum(results.1B$b[t,s,]<input$u.m & results.1B$b[t,s,]>input$l.m),"/",input$iteration.1B,")")
    }}
    colnames(CI)=paste("Number of sampled pens",as.numeric(input$choice.1B))
    rownames(CI)=paste("Total number of sampled fish",as.numeric(input$TSS.1B))
    format(CI,justify = c("right")); t(CI[,ncol(CI):1])
    }})
  output$Table.1B.1 <- renderTable({data.frame(dataset.1())},rownames = TRUE)
  
  tt.1<-eventReactive(input$goButton.1B, {results.1B<-mydata.1B()
  SSperPen=format(results.1B$c)
  colnames(SSperPen)=paste("Number of sampled pens",as.numeric(input$choice.1B))
  rownames(SSperPen)=paste("Total number of sampled fish",as.numeric(input$TSS.1B))
  format(SSperPen,justify = c("right")); t(SSperPen[,ncol(SSperPen):1])
  })
  output$TT.1<-renderTable({tt.1()},rownames = TRUE)
  
  ################################################## Scenario 2
  mydata.2B<-reactive({
    if(input$clustering2B=='no'){pen.AB=c(input$ab1.2B,input$ab2.2B)
    TotalSampleSize=input$TSS.2B
    withProgress(message = 'Progress time',value = 0,{N<-10
    for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
      sam.mean=matrix(NA,length(TotalSampleSize)*input$iteration.2B*length(pen.AB)); dim(sam.mean)=c(length(TotalSampleSize),input$iteration.2B,length(pen.AB))
      for(s in 1:length(TotalSampleSize)){for(i in 1:input$iteration.2B){for(a in 1:length(pen.AB)){sam.mean[s,i,a]=mean(sample(rnbinom(fish*10,size=kappa,mu=pen.AB[a]),TotalSampleSize[s]))}}}}})
    results.2B<-list(h=sam.mean,z=sam.mean)} else {choice.2B=as.numeric(input$choice.2B)
    TotalSampleSize=as.numeric(input$TSS.2B) 
    SampleSizePerPen=matrix(NA,length(TotalSampleSize),length(choice.2B))
    for(s in 1:length(TotalSampleSize)){for(c in 1:length(choice.2B)){
      SampleSizePerPen[s,c]=round(TotalSampleSize[s]/choice.2B[c])
    }}
    gap=seq(0,input$ab1.2B,length.out=200)
    rho=numeric()
    g=1;  rho1=0
    while(rho1<input$icc1.2B & g<length(gap)){
      gg=seq(0,gap[g],length.out=round(input$pen.2B/2))
      if(input$pen.2B%%2==0){pen.AB1=c(rev(input$ab1.2B-gg),input$ab1.2B+gg)} else {
        pen.AB1=c(rev(input$ab1.2B-gg),input$ab1.2B+gg[-1])}
      PP=matrix(NA,fish,input$pen.2B)
      for(p in 1:input$pen.2B){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB1[p])}
      mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
      rho1=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
      g=g+1}
    if(rho1==input$icc1.2B){rho=rho1; pen.AB=pen.AB1} else{g=1;    rho2=0
    while(rho2<input$icc1.2B & g<length(gap)){
      if(input$pen.2B%%2==0){pen.AB2=c(rep(input$ab1.2B-gap[g],round(input$pen.2B/2)),rep(input$ab1.2B+gap[g],round(input$pen.2B/2)))} else {
        pen.AB2=c(rep(input$ab1.2B-gap[g],round(input$pen.2B/2)),rep(input$ab1.2B+gap[g],round(input$pen.2B/2))[-1])}
      PP=matrix(NA,fish,input$pen.2B)
      for(p in 1:input$pen.2B){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB2[p])}
      mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
      rho2=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
      g=g+1}  
    rho=rho2; pen.AB=pen.AB2}
    pen.AB11=pen.AB
    
    gap=seq(0,input$ab2.2B,length.out=200)
    rho=numeric()
    g=1;  rho1=0
    while(rho1<input$icc2.2B & g<length(gap)){
      gg=seq(0,gap[g],length.out=round(input$pen.2B/2))
      if(input$pen.2B%%2==0){pen.AB1=c(rev(input$ab2.2B-gg),input$ab2.2B+gg)} else {
        pen.AB1=c(rev(input$ab2.2B-gg),input$ab2.2B+gg[-1])}
      PP=matrix(NA,fish,input$pen.2B)
      for(p in 1:input$pen.2B){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB1[p])}
      mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
      rho1=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
      g=g+1}
    if(rho1==input$icc2.2B){rho=rho1; pen.AB=pen.AB1} else{g=1;    rho2=0
    while(rho2<input$icc2.2B & g<length(gap)){
      if(input$pen.2B%%2==0){pen.AB2=c(rep(input$ab2.2B-gap[g],round(input$pen.2B/2)),rep(input$ab2.2B+gap[g],round(input$pen.2B/2)))} else {
        pen.AB2=c(rep(input$ab2.2B-gap[g],round(input$pen.2B/2)),rep(input$ab2.2B+gap[g],round(input$pen.2B/2))[-1])}
      PP=matrix(NA,fish,input$pen.2B)
      for(p in 1:input$pen.2B){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB2[p])}
      mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
      rho2=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
      g=g+1}  
    rho=rho2; pen.AB=pen.AB2}
    pen.AB22=pen.AB
    
    pen.AB=matrix(NA,input$pen.2B,2); pen.AB[,1]=pen.AB11; pen.AB[,2]=pen.AB22
    withProgress(message = 'Progress time',value = 0,{N<-10
    for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
      sam.mean=matrix(NA,length(choice.2B),length(TotalSampleSize)*input$iteration.2B*2); dim(sam.mean)=c(length(choice.2B),length(TotalSampleSize),input$iteration.2B,2)
      for(a in 1:2){for(i in 1:input$iteration.2B){sam.pen=matrix(NA,input$pen.2B*length(choice.2B)*length(TotalSampleSize)); dim(sam.pen)=c(input$pen.2B,length(choice.2B),length(TotalSampleSize))
      for(s in 1:length(TotalSampleSize)){for(c in 1:length(choice.2B)){pp=sort(sample(1:input$pen.2B,choice.2B[c]))
      for(p in pp){sam.pen[p,c,s]=mean(sample(rnbinom(fish,size=kappa,mu=pen.AB[p,a]),SampleSizePerPen[s,c]))}}}
      sam.mean[,,i,a]=colMeans(sam.pen,na.rm=T)}}}}) # withProgress
    results.2B<- list(a=sam.mean,b=SampleSizePerPen,c=pen.AB11,d=pen.AB22)}})
  
  ##################### OUTPUT
  plo.2B.1<-eventReactive(input$goButton.2B, {results.2B<-mydata.2B()
  if(input$clustering2B=='yes'){
    plot(x=1:length(results.2B$c),y=results.2B$c[sample(1:length(results.2B$c),length(results.2B$c))],las=1,xlab="Pen number",ylab="Abundance",main=c("Example abundance of each pen","with the supposed true abundance",sub="True Abundance 1"))
    abline(h=input$ab1.2B,col='gray',lty=2)
  }})
  output$Plot.2B.1<-renderPlot({plo.2B.1()})
  
  plo.2B.2<-eventReactive(input$goButton.2B, {results.2B<-mydata.2B()
  if(input$clustering2B=='yes'){
    plot(x=1:length(results.2B$d),y=results.2B$d[sample(1:length(results.2B$d),length(results.2B$d))],las=1,xlab="Pen number",ylab="Abundance",main=c("Example abundance of each pen","with the supposed true abundance",sub="True Abundance 2"))
    abline(h=input$ab2.2B,col='gray',lty=2)
  }})
  output$Plot.2B.2<-renderPlot({plo.2B.2()})
  
  plo.2B.3<-eventReactive(input$goButton.2B, {results.2B<-mydata.2B()
  if(input$clustering2B=='no'){
    dt=melt(results.2B$h); dt[,1]=as.numeric(input$TSS.2B)[dt[,1]]
    dt[,3]=ifelse(dt[,3]==1,"Abundance1","Abundance2")
    colnames(dt)=c("NumberOfTotalSampledFish","Iteration","ABUNDANCE","Abundance")
    ggplot(dt, aes(x=Abundance,colour=ABUNDANCE))+scale_color_manual(values=c('red','blue'))+geom_density(aes(group=ABUNDANCE),alpha=0.25)+facet_wrap(~NumberOfTotalSampledFish)+ggtitle("Total number of sampled fish")+theme(plot.title = element_text(size=18))+geom_vline(xintercept =c(input$ab1.2B,input$ab2.2B), linetype = "dashed", colour = rep(c(2,4),length(table(dt[,1]))))+ylab("Density")+theme(legend.position="bottom",text=element_text(size=18),axis.text.x=element_text(face="bold",size=18),axis.text.y=element_text(face="bold",size=18))+guides(col=guide_legend(""))
  } else {
    dt=melt(results.2B$a); dt[,1]=as.numeric(input$choice.2B)[dt[,1]]; dt[,2]=as.numeric(input$TSS.2B)[dt[,2]]; dt[,4]=ifelse(dt[,4]==1,"Abundance1","Abundance2")
    colnames(dt)=c("NumberOfSampledPens","NumberOfTotalSampledFish","Iteration","ABUNDANCE","Abundance")
    ggplot(dt,aes(x=Abundance,colour=ABUNDANCE))+scale_color_manual(values=c('red','blue'))+geom_density(aes(group=ABUNDANCE),alpha=0.25)+facet_grid(NumberOfSampledPens~NumberOfTotalSampledFish)+ggtitle("Total number of sampled fish")+theme(plot.title = element_text(size=18))+geom_vline(xintercept=c(input$ab1.2B,input$ab2.2B),linetype="dashed",colour=rep(c(2,4),length(table(dt[,2]))*length(table(dt[,1]))))+ylab("Density")+theme(legend.position="bottom",text=element_text(size=18),axis.text.x=element_text(face="bold",size=18),axis.text.y=element_text(face="bold",size=18))+guides(col=guide_legend(""))}})
  output$Plot.2B.3<-renderPlot({plo.2B.3()})
  
  dataset.2<-eventReactive(input$goButton.2B,{results.2B<-mydata.2B()
  if(input$clustering2B=='no'){abc=matrix(NA,length(as.numeric(input$TSS.2B)),1)
  for(s in 1:length(as.numeric(input$TSS.2B))){abc[s]=ifelse(input$ab1.2B<input$ab2.2B, paste(round(mean(results.2B$h[s,,1]<results.2B$h[s,,2]),2),"(",sum(results.2B$h[s,,1]<results.2B$h[s,,2]),"/",input$iteration.2B,")"),paste(round(mean(results.2B$h[s,,1]>results.2B$h[s,,2]),2),"(",sum(results.2B$h[s,,1]>results.2B$h[s,,2]),"/",input$iteration.2B,")"))}
  colnames(abc)=c("Probability"); rownames(abc)=paste("Total Sample Size",as.numeric(input$TSS.2B))
  abc} else {abc=matrix(NA,length(as.numeric(input$TSS.2B)),length(as.numeric(input$choice.2B)))
  for(s in 1:length(as.numeric(input$TSS.2B))){for(c in 1:length(as.numeric(input$choice.2B))){
    abc[s,c]=ifelse(input$ab1.2B<input$ab2.2B, paste(round(mean(results.2B$a[c,s,,1]<results.2B$a[c,s,,2]),2),"(",sum(results.2B$a[c,s,,1]<results.2B$a[c,s,,2]),"/",input$iteration.2B,")"), paste(round(mean(results.2B$a[c,s,,1]>results.2B$a[c,s,,2]),2),"(",sum(results.2B$a[c,s,,1]>results.2B$a[c,s,,2]),"/",input$iteration.2B,")"))}}
  colnames(abc)=paste("Number of sampled pens",as.numeric(input$choice.2B)); rownames(abc)=paste("Total number of sampled fish",as.numeric(input$TSS.2B))
  t(abc)}})
  output$Table.2B.1<-renderTable({data.frame(dataset.2())},rownames = TRUE)
  
  tt.2<-eventReactive(input$goButton.2B, {results.2B<-mydata.2B()
  SSperPen=format(results.2B$b)
  colnames(SSperPen)=paste("Number of sampled pens",as.numeric(input$choice.2B))
  rownames(SSperPen)=paste("Total number of sampled fish",as.numeric(input$TSS.2B))
  t(SSperPen)})
  output$TT.2<-renderTable({tt.2()},rownames = TRUE)
  
  ############################## Scenario 3
  mydata.3B<-reactive({if(input$clustering3B=='no'){TotalSampleSize=as.numeric(input$TSS.3B)
  withProgress(message = 'Progress time',value = 0,{N<-10
  for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
    sam.mean=matrix(NA,length(TotalSampleSize),input$iteration.3B)
    for(i in 1:input$iteration.3B){for(s in 1:length(TotalSampleSize)){sam.mean[s,i]=mean(sample(rnbinom(fish*10,size=kappa,mu=input$ab.3B),TotalSampleSize[s]))}}}})
  results.3B<-list(d=sam.mean) } else {choice.3B=as.numeric(input$choice.3B)
  gap=seq(0,input$ab.3B,length.out=200)
  rho=numeric()
  g=1;  rho1=0
  while(rho1<input$icc.3B & g<length(gap)){
    gg=seq(0,gap[g],length.out=round(input$pen.3B/2))
    if(input$pen.3B%%2==0){pen.AB1=c(rev(input$ab.3B-gg),input$ab.3B+gg)} else {pen.AB1=c(rev(input$ab.3B-gg),input$ab.3B+gg[-1])}
    PP=matrix(NA,fish,input$pen.3B)
    for(p in 1:input$pen.3B){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB1[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    rho1=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
    g=g+1}
  if(rho1==input$icc.3B){rho=rho1; pen.AB=pen.AB1} else{g=1;    rho2=0
  while(rho2<input$icc.3B & g<length(gap)){
    if(input$pen.3B%%2==0){pen.AB2=c(rep(input$ab.3B-gap[g],round(input$pen.3B/2)),rep(input$ab.3B+gap[g],round(input$pen.3B/2)))} else {
      pen.AB2=c(rep(input$ab.3B-gap[g],round(input$pen.3B/2)),rep(input$ab.3B+gap[g],round(input$pen.3B/2))[-1])}
    PP=matrix(NA,fish,input$pen.3B)
    for(p in 1:input$pen.3B){PP[,p]=rnbinom(fish,size=2.19,mu=pen.AB2[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    rho2=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)  
    g=g+1}  
  rho=rho2; pen.AB=pen.AB2}
  
  TotalSampleSize=as.numeric(input$TSS.3B)
  SampleSizePerPen=matrix(NA,length(TotalSampleSize),length(choice.3B))
  for(s in 1:length(TotalSampleSize)){for(c in 1:length(choice.3B)){SampleSizePerPen[s,c]=round(TotalSampleSize[s]/choice.3B[c])}}
  withProgress(message = 'Progress time',value = 0,{N<-10
  for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
    sam.mean=matrix(NA,length(TotalSampleSize)*length(choice.3B)*input$iteration.3B); dim(sam.mean)=c(length(TotalSampleSize),length(choice.3B),input$iteration.3B)
    for(s in 1:length(TotalSampleSize)){for(i in 1:input$iteration.3B){sam.pen=matrix(NA,input$pen.3B,length(choice.3B))
    for(c in 1:length(choice.3B)){pp=sort(sample(1:input$pen.3B,choice.3B[c]))
    for(p in pp){sam.pen[p,c]=mean(sample(rnbinom(fish,size=kappa,mu=pen.AB[p]),SampleSizePerPen[s,c]))}}
    sam.mean[s,,i]=colMeans(sam.pen,na.rm=T)}}}})
  results.3B<-list(a=sam.mean,b=SampleSizePerPen,c=pen.AB)}})
  
  #######################  OUTPUT 3B: 
  
  plo.3B.2<-eventReactive(input$goButton.3B, {results.3B<-mydata.3B()
  if(input$clustering3B=='yes'){
    plot(x=1:length(results.3B$c),y=results.3B$c[sample(1:length(results.3B$c),length(results.3B$c))],las=1,xlab="Pen number",ylab="Abundance",main=c("Example abundance of each pen","with the supposed ICC"))
    abline(h=input$ab.3B,col='gray',lty=2)}})
  output$Plot.3B.2<-renderPlot({plo.3B.2()})
  
  plo.3B.1<-eventReactive(input$goButton.3B, {results.3B<-mydata.3B()
  if(input$clustering3B=='no'){dt=melt(results.3B$d)
  dt[,1]=as.numeric(input$TSS.3B)[dt[,1]]; dt[,1]=as.factor(dt[,1]); dt[,2]=as.factor(dt[,2])
  colnames(dt)=c("NumberOfTotalSampledFish","Iteration","Abundance")
  ggplot(dt, aes(x=Abundance,y=NumberOfTotalSampledFish))+geom_density_ridges()+geom_vline(xintercept=input$ab.3B)+geom_vline(xintercept=input$th.3B,linetype="dashed",colour="red")+ylab(ifelse(length(as.numeric(input$TSS.3B))==1,"Density","Total number of sampled fish"))+theme(text=element_text(size=18))+theme(axis.text.x=element_text(face="bold",size=18),axis.text.y = element_text(face="bold",size=18))
  } else {dt=melt(results.3B$a)
  dt[,1]=as.numeric(input$TSS.3B)[dt[,1]]; dt[,2]=as.numeric(input$choice.3B)[dt[,2]]; dt[,2]=as.factor(dt[,2])
  colnames(dt)=c("NumberOfTotalSampledFish","NumberOfSampledPens","Iteration","Abundance")
  ggplot(dt, aes(x=Abundance, y=NumberOfSampledPens)) +geom_density_ridges()+facet_wrap(~NumberOfTotalSampledFish)+ggtitle("Total number of sampled fish")+theme(plot.title = element_text(size=18))+
    geom_vline(xintercept=mean(extract(input$ab.3B)))+geom_vline(xintercept=input$th.3B,linetype="dashed",colour="red")+
    ylab(ifelse(length(as.numeric(input$choice.3B))==1,"Density","Number of sampled pens"))+
    theme(legend.position="bottom",legend.title=element_text("Number of sampled pens",size=18,face="bold"))+theme(text=element_text(size=18))+
    theme(axis.text.x=element_text(face="bold",size=18),axis.text.y=element_text(face="bold",size=18))
  }})
  output$Plot.3B.1 <- renderPlot({plo.3B.1()})
  
  dataset.3<-eventReactive(input$goButton.3B,{results.3B<-mydata.3B()
  if(input$clustering3B=='no'){cutline=matrix(NA,length(as.numeric(input$TSS.3B)),1)
  for(s in 1:length(as.numeric(input$TSS.3B))){cutline[s]=ifelse(input$ab.3B>=input$th.3B,paste(round(mean(results.3B$d[s,]>=input$th.3B),2),"(",sum(results.3B$d[s,]>=input$th.3B),"/",input$iteration.3B,")"),paste(round(mean(results.3B$d[s,]<input$th.3B),2),"(",sum(results.3B$d[s,]<input$th.3B),"/",input$iteration.3B,")"))}
  cutline=matrix(cutline[nrow(cutline):1,],ncol=1)
  colnames(cutline)="Probability"; rownames(cutline)=rev(paste("Total number of sampled fish",as.numeric(input$TSS.3B)))
  cutline} else {
    cutline=matrix(NA,length(as.numeric(input$TSS.3B)),length(as.numeric(input$choice.3B)))
    for(c in 1:length(as.numeric(input$choice.3B))){for(s in 1:length(as.numeric(input$TSS.3B))){cutline[s,c]=ifelse(mean(extract(input$ab.3B))>=input$th.3B,paste(round(mean(results.3B$a[s,c,]>=input$th.3B),2),"(",sum(results.3B$a[s,c,]>=input$th.3B),"/",input$iteration.3B,")"),paste(round(mean(results.3B$a[s,c,]<input$th.3B),2),"(",sum(results.3B$a[s,c,]<input$th.3B),"/",input$iteration.3B,")"))}}
    colnames(cutline)=paste("Number of sampled pens",as.numeric(input$choice.3B)); rownames(cutline)=paste("Total number of sampled fish",as.numeric(input$TSS.3B))
    format(cutline,justify=c("right"))
    t(cutline[,ncol(cutline):1])}
  })
  output$Table.3B.1 <- renderTable({data.frame(dataset.3())},rownames = TRUE)
  
  tt.3<-eventReactive(input$goButton.3B,{results.3B<-mydata.3B()
  SSperPen=format(results.3B$b)
  colnames(SSperPen)=paste("Number of sampled pens",as.numeric(input$choice.3B));
  rownames(SSperPen)=paste("Total number of sampled fish",as.numeric(input$TSS.3B))
  t(SSperPen[,ncol(SSperPen):1])
  })
  output$TT.3<-renderTable({tt.3()},rownames = TRUE)
}
