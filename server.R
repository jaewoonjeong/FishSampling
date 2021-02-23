server <-function(input, output,session){
  extract <- function(text) {text <- gsub(" ", "", text); split <- strsplit(text, ",", fixed = FALSE)[[1]]; as.numeric(split)}
  
  kappa=2.19 # dispersion parameter of negative binomial distribution
  fish=2000 # number of fish per pen
  
  #######################  SCENARIO 1B: mydata1B
  kappa=2.19 # dispersion parameter of negative binomial distribution
  fish=1000 # number of fish per pen
  mydata.1B<-reactive({if(input$clustering1B=='no'){TotalSampleSize=as.numeric(input$TSS.1B)
  withProgress(message = 'Progress time',value = 0,{N<-10
  for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
    sam.mean=matrix(NA,length(TotalSampleSize),input$iteration.1B)
    for(i in 1:input$iteration.1B){for(s in 1:length(TotalSampleSize)){sam.mean[s,i]=mean(sample(rnbinom(fish*10,size=kappa,mu=input$ab.1B),TotalSampleSize[s]))}}}})
  results.1B<-list(s=sam.mean)} else {choice.1B=as.numeric(input$choice.1B)
  pen.AB=extract(input$penAB.1B)
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
  results.1B<-list(a=choice.1B, b=sam.mean, c=SampleSizePerPen, d=pen.AB)}})
  #######################  OUTPUT 1B: 
  output$ui.1B.1 <-renderUI({textInput("penAB.1B", "Abundance of Pens", value = NULL, placeholder = 'Use a comma between abundance')})
  output$Text.1B.1<-renderText({ABofPen=extract(input$penAB.1B)
  out<-round(sum(ABofPen)/length(ABofPen),2)
  ifelse(length(ABofPen)<1, "NULL", out)})
  output$Text.1B.2<-renderText({validate(need(input$penAB.1B, ""))
    pen.AB=extract(input$penAB.1B)
    PP=matrix(NA,1000,input$pen.1B)
    for(p in 1:input$pen.1B){PP[,p]=rnbinom(1000,size=2.19,mu=pen.AB[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    rho=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)
    rho})
  output$Plot.1B.2<-renderPlot({validate(need(input$penAB.1B, ""))
    pen.AB=extract(input$penAB.1B)
    plot(pen.AB, xlab="Pen number",ylab="Abundance",main="",ylim=c(0,max(pen.AB)),cex=1.5,las=1); abline(h=mean(pen.AB),lty=2)})
  
  plo.1B.3<-eventReactive(input$goButton.1B, {results.1B<-mydata.1B()
  if(input$clustering1B=='no'){  md=melt(results.1B$s) 
  md[,1]=as.numeric(input$TSS.1B)[md[,1]]; md[,1]=as.factor(md[,1]); md[,2]=as.factor(md[,2])
  colnames(md)=c("NumberOfTotalSampledFish","Iteration","Abundance")
  ggplot(md,aes(x=Abundance,y=NumberOfTotalSampledFish))+geom_density_ridges()+geom_vline(xintercept=input$ab.1B)+geom_vline(xintercept=c(input$u.m,input$l.m),colour=c("green"), linetype = "dashed")+ylab(ifelse(length(as.numeric(input$TSS.1B))==1,"Density","Total number of sampled fish"))+theme(text=element_text(size=18))+theme(axis.text.x=element_text(face="bold",size=18),axis.text.y = element_text(face="bold",size=18))
  } else {md=melt(results.1B$b)
  md[,1]=as.numeric(input$choice.1B)[md[,1]]; md[,1]=as.factor(md[,1]); md[,2]=as.numeric(input$TSS.1B)[md[,2]]; md[,2]=as.factor(md[,2]); md[,3]=as.factor(md[,3])
  colnames(md)=c("NumberOfSampledPens","NumberOfTotalSampledFish","Iteration","Abundance")
  ggplot(md, aes(x=Abundance,y=NumberOfSampledPens))+geom_density_ridges()+facet_wrap(~NumberOfTotalSampledFish)+theme(legend.position="bottom",legend.title=element_text("Number of sampled pens",size=18,face="bold"))+ggtitle("Total number of sampled fish")+ylab(ifelse(length(as.numeric(input$choice.1B))==1,"Density","Number of sampled pens"))+theme(plot.title = element_text(size=18), text=element_text(size=18))+theme(axis.text.x=element_text(face="bold",size=18),axis.text.y=element_text(face="bold",size=18))+geom_vline(xintercept=mean(extract(input$penAB.1B)))+geom_vline(xintercept=c(input$u.m,input$l.m),colour=c("green"),linetype="dashed")}})
  output$Plot.1B.3<-renderPlot({plo.1B.3()})
  
  dataset.1<-eventReactive(input$goButton.1B,{results.1B<-mydata.1B()
  if(input$clustering1B=='no'){
    CI=matrix(NA,length(as.numeric(input$TSS.1B)),1)
    for(s in 1:length(as.numeric(input$TSS.1B))){CI[s]=round(mean(results.1B$s[s,]<input$u.m & results.1B$s[s,]>input$l.m),2)}
    colnames(CI)="Probability"; rownames(CI)=paste("Total Sample Size",as.numeric(input$TSS.1B))   
    CI} else {CI=matrix(NA,length(as.numeric(input$TSS.1B)),length(as.numeric(input$choice.1B)))
    for(s in 1:length(as.numeric(input$TSS.1B))){for(t in 1:length(as.numeric(input$choice.1B))){CI[s,t]=round(mean(results.1B$b[t,s,]<input$u.m & results.1B$b[t,s,]>input$l.m),2)}}
    colnames(CI)=paste("Number of sampled pens",as.numeric(input$choice.1B))
    rownames(CI)=paste("Total number of sampled fish",as.numeric(input$TSS.1B))
    CI}})
  output$Table.1B.1 <- renderTable({data.frame(dataset.1())},rownames = TRUE)
  
  tt.1<-eventReactive(input$goButton.1B, {results.1B<-mydata.1B()
  SSperPen=format(results.1B$c)
  colnames(SSperPen)=c(paste("Number of sampled pens",as.numeric(input$choice.1B)[1]),as.numeric(input$choice.1B)[-1]);
  rownames(SSperPen)=paste("Total number of sampled fish",as.numeric(input$TSS.1B)); format(SSperPen,justify = c("right")); print(SSperPen)})
  output$TT.1<-renderTable({tt.1()},rownames = TRUE)
  
  ################################################## Scenario 2
  mydata.2B<-reactive({
    if(input$clustering2B=='no'){pen.AB=c(input$ab1.2B,input$ab2.2B)
    TotalSampleSize=input$TSS.2B
    withProgress(message = 'Progress time',value = 0,{N<-10
    for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
      sam.mean=matrix(NA,length(TotalSampleSize)*input$iteration.2B*length(pen.AB)); dim(sam.mean)=c(length(TotalSampleSize),input$iteration.2B,length(pen.AB))
      for(s in 1:length(TotalSampleSize)){for(i in 1:input$iteration.2B){for(a in 1:length(pen.AB)){sam.mean[s,i,a]=mean(sample(rnbinom(fish*10,size=kappa,mu=pen.AB[a]),TotalSampleSize[s]))}}}}})
    results.2B<-list(h=sam.mean,z=sam.mean)} else {
      choice.2B=as.numeric(input$choice.2B)
      TotalSampleSize=as.numeric(input$TSS.2B) 
      SampleSizePerPen=matrix(NA,length(TotalSampleSize),length(choice.2B))
      for(s in 1:length(TotalSampleSize)){for(c in 1:length(choice.2B)){
        SampleSizePerPen[s,c]=round(TotalSampleSize[s]/choice.2B[c])
      }}
      
      pen.AB1=extract(input$penAB1.2B); pen.AB2=extract(input$penAB2.2B)
      pen.AB=matrix(NA,input$pen.2B,2); pen.AB[,1]=pen.AB1; pen.AB[,2]=pen.AB2
      withProgress(message = 'Progress time',value = 0,{N<-10
      for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
        sam.mean=matrix(NA,length(choice.2B),length(TotalSampleSize)*input$iteration.2B*2); dim(sam.mean)=c(length(choice.2B),length(TotalSampleSize),input$iteration.2B,2)
        for(a in 1:2){for(i in 1:input$iteration.2B){sam.pen=matrix(NA,input$pen.2B*length(choice.2B)*length(TotalSampleSize)); dim(sam.pen)=c(input$pen.2B,length(choice.2B),length(TotalSampleSize))
        for(s in 1:length(TotalSampleSize)){for(c in 1:length(choice.2B)){pp=sort(sample(1:input$pen.2B,choice.2B[c]))
        for(p in pp){sam.pen[p,c,s]=mean(sample(rnbinom(fish,size=kappa,mu=pen.AB[p,a]),SampleSizePerPen[s,c]))}}}
        sam.mean[,,i,a]=colMeans(sam.pen,na.rm=T)}}}}) # withProgress
      results.2B<- list(a=sam.mean,b=SampleSizePerPen)}})
  
  output$ui.2B.1 <-renderUI({textInput("penAB1.2B", "Abundance1 of Pens", value = NULL, placeholder = 'Use a comma between abundance')})
  output$ui.2B.2 <-renderUI({textInput("penAB2.2B", "Abundance2 of Pens (Abundance of Interest)", value = NULL, placeholder = 'Use a comma between abundance')})
  output$Text.2B.1<-renderText({out<-round(sum(extract(input$penAB1.2B))/length(extract(input$penAB1.2B)),2); ifelse(length(extract(input$penAB1.2B))<1, "NULL", out)})
  output$Text.2B.2<-renderText({out<-round(sum(extract(input$penAB2.2B))/length(extract(input$penAB2.2B)),2); ifelse(length(extract(input$penAB2.2B))<1, "NULL", out)})
  output$Text.2B.3<-renderText({validate(need(input$penAB1.2B, ""))
    PP=matrix(NA,1000,input$pen.2B); for(p in 1:input$pen.2B){PP[,p]=rnbinom(1000,size=2.19,mu=extract(input$penAB1.2B)[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2]); round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)})
  output$Text.2B.4<-renderText({validate(need(input$penAB2.2B, ""))
    PP=matrix(NA,1000,input$pen.2B); for(p in 1:input$pen.2B){PP[,p]=rnbinom(1000,size=2.19,mu=extract(input$penAB2.2B)[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)})
  output$Plot.2B.1<-renderPlot({validate(need(input$penAB1.2B, ""))
    pen.AB=extract(input$penAB1.2B); plot(pen.AB,main="Abundance1", xlab="Pen number",ylab="Abundance",ylim=c(0,max(pen.AB)),cex=1.5,las=1); abline(h=mean(pen.AB),lty=2)})
  output$Plot.2B.2<-renderPlot({validate(need(input$penAB2.2B, ""))
    pen.AB=extract(input$penAB2.2B)
    plot(pen.AB,main="Abundance2", xlab="Pen number",ylab="Abundance",ylim=c(0,max(pen.AB)),cex=1.5,las=1); abline(h=mean(pen.AB),lty=2)})
  
  plo.2B.3<-eventReactive(input$goButton.2B, {results.2B<-mydata.2B()
  if(input$clustering2B=='no'){dt=melt(results.2B$h); dt[,1]=as.numeric(input$TSS.2B)[dt[,1]]
  dt[,3]=ifelse(dt[,3]==1,"Abundance1","Abundance2")
  colnames(dt)=c("NumberOfTotalSampledFish","Iteration","ABUNDANCE","Abundance")
  ggplot(dt, aes(x=Abundance,colour=ABUNDANCE))+scale_color_manual(values=c('red','blue'))+geom_density(aes(group=ABUNDANCE),alpha=0.25)+facet_wrap(~NumberOfTotalSampledFish)+ggtitle("Total number of sampled fish")+theme(plot.title = element_text(size=18))+geom_vline(xintercept =c(input$ab1.2B,input$ab2.2B), linetype = "dashed", colour = rep(c(2,4),length(table(dt[,1]))))+ylab("Density")+theme(legend.position="bottom",text=element_text(size=18),axis.text.x=element_text(face="bold",size=18),axis.text.y=element_text(face="bold",size=18))} else {dt=melt(results.2B$a); dt[,1]=as.numeric(input$choice.2B)[dt[,1]]; dt[,2]=as.numeric(input$TSS.2B)[dt[,2]]; dt[,4]=ifelse(dt[,4]==1,"Abundance1","Abundance2")
  colnames(dt)=c("NumberOfSampledPens","NumberOfTotalSampledFish","Iteration","ABUNDANCE","Abundance")
  ggplot(dt,aes(x=Abundance,colour=ABUNDANCE))+scale_color_manual(values=c('red','blue'))+geom_density(aes(group=ABUNDANCE),alpha=0.25)+facet_grid(NumberOfTotalSampledFish~NumberOfSampledPens)+ggtitle("Number of sampled pens")+theme(plot.title = element_text(size=18))+geom_vline(xintercept=c(mean(extract(input$penAB1.2B)),mean(extract(input$penAB2.2B))),linetype="dashed",colour=rep(c(2,4),length(table(dt[,2]))*length(table(dt[,1]))))+ylab("Density")+theme(legend.position="bottom",text=element_text(size=18),axis.text.x=element_text(face="bold",size=18),axis.text.y=element_text(face="bold",size=18))}})
  output$Plot.2B.3<-renderPlot({plo.2B.3()})
  
  dataset.2<-eventReactive(input$goButton.2B,{results.2B<-mydata.2B()
  if(input$clustering2B=='no'){abc=matrix(NA,length(as.numeric(input$TSS.2B)),1)
  for(s in 1:length(as.numeric(input$TSS.2B))){abc[s]=ifelse(input$ab1.2B<input$ab2.2B, mean(results.2B$h[s,,1]<results.2B$h[s,,2]), mean(results.2B$h[s,,1]>results.2B$h[s,,2]))}
  colnames(abc)=c("Probability"); rownames(abc)=paste("Total Sample Size",as.numeric(input$TSS.2B))
  abc} else {abc=matrix(NA,length(as.numeric(input$TSS.2B)),length(as.numeric(input$choice.2B)))
  for(s in 1:length(as.numeric(input$TSS.2B))){for(c in 1:length(as.numeric(input$choice.2B))){abc[s,c]=ifelse(mean(extract(input$penAB1.2B))<mean(extract(input$penAB2.2B)), mean(results.2B$a[c,s,,1]<results.2B$a[c,s,,2]), mean(results.2B$a[c,s,,1]>results.2B$a[c,s,,2]))}}
  colnames(abc)=paste("Number of sampled pens",as.numeric(input$choice.2B))
  rownames(abc)=paste("Total number of sampled fish",as.numeric(input$TSS.2B))
  abc}})
  output$Table.2B.1<-renderTable({data.frame(dataset.2())},rownames = TRUE)
  
  tt.2<-eventReactive(input$goButton.2B, {results.2B<-mydata.2B()
  SSperPen=format(results.2B$b)
  colnames(SSperPen)=paste("Number of sampled pens",as.numeric(input$choice.2B))
  rownames(SSperPen)=paste("Total number of sampled fish",as.numeric(input$TSS.2B))
  print(SSperPen)})
  output$TT.2<-renderTable({tt.2()},rownames = TRUE)
  
  ############################## Scenario 3
  mydata.3B<-reactive({if(input$clustering3B=='no'){
    TotalSampleSize=as.numeric(input$TSS.3B)
    withProgress(message = 'Progress time',value = 0,{N<-10
    for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
      sam.mean=matrix(NA,length(TotalSampleSize),input$iteration.3B)
      for(i in 1:input$iteration.3B){for(s in 1:length(TotalSampleSize)){sam.mean[s,i]=mean(sample(rnbinom(fish*10,size=kappa,mu=input$ab.3B),TotalSampleSize[s]))}}}})
    results.3B<-list(d=sam.mean) } else {
      choice.3B=as.numeric(input$choice.3B)
      pen.AB=extract(input$penAB.3B)
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
      results.3B<-list(a=sam.mean,b=SampleSizePerPen)}})
  #######################  OUTPUT 3B: 
  output$ui.3B.1 <-renderUI({textInput("penAB.3B", "Abundance of Pens", value = NULL, placeholder = 'Use a comma between abundance')})
  output$Text.3B.1<-renderText({ABofPen=extract(input$penAB.3B)
  out<-round(sum(ABofPen)/length(ABofPen),2)
  ifelse(length(ABofPen)<1, "NULL", out)})
  output$Text.3B.2<-renderText({validate(need(input$penAB.3B, ""))
    pen.AB=extract(input$penAB.3B)
    PP=matrix(NA,1000,input$pen.3B)
    for(p in 1:input$pen.3B){PP[,p]=rnbinom(1000,size=2.19,mu=pen.AB[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    rho=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)
    rho})
  output$Plot.3B.2<-renderPlot({validate(need(input$penAB.3B, ""))
    pen.AB=extract(input$penAB.3B)
    plot(pen.AB, xlab="Pen number",ylab="Abundance",main="",ylim=c(0,max(pen.AB)),cex=1.5,las=1); abline(h=mean(pen.AB),lty=2)})
  
  plo.3B.1<-eventReactive(input$goButton.3B, {results.3B<-mydata.3B()
  
  if(input$clustering3B=='no'){dt=melt(results.3B$d)
  dt[,1]=as.numeric(input$TSS.3B)[dt[,1]]; dt[,1]=as.factor(dt[,1]); dt[,2]=as.factor(dt[,2])
  colnames(dt)=c("NumberOfTotalSampledFish","Iteration","Abundance")
  ggplot(dt, aes(x=Abundance,y=NumberOfTotalSampledFish))+geom_density_ridges()+geom_vline(xintercept=input$ab.3B)+geom_vline(xintercept=input$th.3B,linetype="dashed",colour="red")+ylab(ifelse(length(as.numeric(input$TSS.3B))==1,"Density","Total number of sampled fish"))+theme(text=element_text(size=18))+theme(axis.text.x=element_text(face="bold",size=18),axis.text.y = element_text(face="bold",size=18))
  } else {dt=melt(results.3B$a)
  dt[,1]=as.numeric(input$TSS.3B)[dt[,1]]; dt[,2]=as.numeric(input$choice.3B)[dt[,2]]; dt[,2]=as.factor(dt[,2])
  colnames(dt)=c("NumberOfTotalSampledFish","NumberOfSampledPens","Iteration","Abundance")
  ggplot(dt, aes(x=Abundance, y=NumberOfSampledPens)) +geom_density_ridges()+facet_wrap(~NumberOfTotalSampledFish)+ggtitle("Number of sampled pens")+theme(plot.title = element_text(size=18))+
    geom_vline(xintercept=mean(extract(input$penAB.3B)))+geom_vline(xintercept=input$th.3B,linetype="dashed",colour="red")+
    ylab(ifelse(length(as.numeric(input$choice.3B))==1,"Density","Number of sampled pens"))+
    theme(legend.position="bottom",legend.title=element_text("Number of sampled pens",size=18,face="bold"))+theme(text=element_text(size=18))+
    theme(axis.text.x=element_text(face="bold",size=18),axis.text.y=element_text(face="bold",size=18))
  }})
  output$Plot.3B.1 <- renderPlot({plo.3B.1()})
  
  dataset.3<-eventReactive(input$goButton.3B,{results.3B<-mydata.3B()
  if(input$clustering3B=='no'){cutline=matrix(NA,length(as.numeric(input$TSS.3B)),1)
  for(s in 1:length(as.numeric(input$TSS.3B))){cutline[s]=ifelse(input$ab.3B>=input$th.3B,mean(results.3B$d[s,]>=input$th.3B),mean(results.3B$d[s,]<input$th.3B))}
  colnames(cutline)="Probability"; rownames(cutline)=paste("Total number of sampled fish",as.numeric(input$TSS.3B))
  cutline} else {cutline=matrix(NA,length(as.numeric(input$TSS.3B)),length(as.numeric(input$choice.3B)))
  for(c in 1:length(as.numeric(input$choice.3B))){for(s in 1:length(as.numeric(input$TSS.3B))){cutline[s,c]=ifelse(mean(extract(input$penAB.3B))>=input$th.3B,mean(results.3B$a[s,c,]>=input$th.3B),mean(results.3B$a[s,c,]<input$th.3B))}}
  colnames(cutline)=paste("Number of sampled pens",as.numeric(input$choice.3B)); rownames(cutline)=paste("Total number of sampled fish",as.numeric(input$TSS.3B))
  cutline}})
  output$Table.3B.1 <- renderTable({data.frame(dataset.3())},rownames = TRUE)
  tt.3<-eventReactive(input$goButton.3B,{results.3B<-mydata.3B()
  SSperPen=format(results.3B$b)
  colnames(SSperPen)=c(paste("Number of sampled pens",as.numeric(input$choice.3B)[1]),as.numeric(input$choice.3B)[-1]);
  rownames(SSperPen)=paste("Total number of sampled fish",as.numeric(input$TSS.3B)); print(SSperPen)})
  output$TT.3<-renderTable({tt.3()},rownames = TRUE)
}