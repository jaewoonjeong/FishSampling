server <-function(input, output,session){
  extract <- function(text) {text <- gsub(" ", "", text)
    split <- strsplit(text, ",", fixed = FALSE)[[1]]
    as.numeric(split)}
  
  kappa=2.19 # dispersion parameter of negative binomial distribution
  fish=1000 # number of fish per pen
  
  #######################  SCENARIO 1A
  cl.1A=c(0.8,0.9,0.95,0.99) # confidence level
  mydata.1A<-reactive({
    choice.1A=as.numeric(input$choice.1A)
    pen.sample.size=round(input$total.sample.size.1A/choice.1A) # sample size per pen
    if(input$clustering1A=='no'){pen.AB=rep(input$ab.1A,input$pen.1A)} else {pen.AB=extract(input$penAB.1A)}
    withProgress(message = 'Progress time',value = 0,{N<-10
      for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
        sam.mean=matrix(NA,length(choice.1A),input$iteration.1A)
        for(i in 1:input$iteration.1A){sam.pen=matrix(NA,input$pen.1A,length(choice.1A))
        for(c in 1:length(choice.1A)){pp=sort(sample(1:input$pen.1A,choice.1A[c]))
        for(p in pp){sam.pen[p,c]=mean(sample(rnbinom(fish,size=kappa,mu=pen.AB[p]),pen.sample.size[c]))}}
        sam.mean[,i]=colMeans(sam.pen,na.rm=T)}
        ssm=matrix(NA,length(choice.1A),input$iteration.1A)
        for(c in 1:length(choice.1A)){ssm[c,]=sort(sam.mean[c,],na.last = T)}
        range.ab=matrix(NA,length(cl.1A)*length(choice.1A)*2); dim(range.ab)=c(2,length(choice.1A),length(cl.1A))
        for(c in 1:length(choice.1A)){for(r in 1:length(cl.1A)){range.ab[,c,r]=ssm[c,][round(input$iteration.1A*c((1-cl.1A[r])/2,1-(1-cl.1A[r])/2))]}}}})
    results.1A<-list(a=round(range.ab,2), b=pen.sample.size, c=pen.AB, d=choice.1A)})

  output$ui.1A.1 <-renderUI({textInput("penAB.1A", "Abundance of Pens", value = NULL, placeholder = 'Use a comma between abundance')})
  output$Text.1A.1<-renderText({ABofPen=extract(input$penAB.1A)
    out<-round(sum(ABofPen)/length(ABofPen),2)
    ifelse(length(ABofPen)<1, "NULL", out)})
  output$Text.1A.2<-renderText({validate(need(input$penAB.1A, ""))
    pen.AB=extract(input$penAB.1A)
    PP=matrix(NA,1000,input$pen.1A)
    for(p in 1:input$pen.1A){PP[,p]=rnbinom(1000,size=2.19,mu=pen.AB[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    rho=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2); rho})
  output$Plot.1A.2<-renderPlot({validate(need(input$penAB.1A, ""))
    pen.AB=extract(input$penAB.1A)
    plot(pen.AB, xlab="Pen number",ylab="Abundance",main="",ylim=c(0,max(pen.AB)),cex=1.5,las=1); abline(h=mean(pen.AB),lty=2)})
  plo.1A.1<-eventReactive(input$goButton.1A,{results.1A<-mydata.1A()
    par(oma=c(0,0,0,0),mar=c(4,2,4,9),cex=1.2,pin=c(5.5,3))
    for(r in 1:length(cl.1A)){matplot(lwd=2,t(results.1A$a[,,r]),type="b",xlab="",ylab="",lty=r,col=r,pch=r,xlim=c(1,length(results.1A$d)),ylim=range(results.1A$a,na.rm=T),axes=F); ifelse(r==length(cl.1A),NA,par(new=T))}
    lines(x=1:length(input$choice.1A),y=rep(mean(results.1A$c),length(input$choice.1A)),lty=2)
    axis(1,1:length(results.1A$d),results.1A$d); axis(2,las=1); box(); axis(3,1:length(results.1A$d),results.1A$b); 
    mtext(cex=1.2,"Number of sampled pens",side=1,line=3); mtext(cex=1.2,"Abundance",2,3); mtext(cex=1.2,"Number of sampled fish per pen",3,3)
    legend(lwd=2,"right",border = "black",horiz = FALSE,xpd=TRUE,inset=c(-0.3,0),bty="o",legend=cl.1A,title="Confidence Level",pch=1:length(cl.1A),lty=1:length(cl.1A), col=1:length(cl.1A))})
  output$Plot.1A.1<-renderPlot({if(length(input$choice.1A)==1){NA} else {plo.1A.1()}})
  tab.1A.1<-eventReactive(input$goButton.1A, {results.1A<-mydata.1A()
    lm=results.1A$a[1,,]; um=matrix(NA,length(results.1A$d),length(cl.1A))
    for(i in 1:length(results.1A$d)){um[i,]=results.1A$a[2,(length(results.1A$d):1)[i],]}
    um=rev(um); margins=matrix(c(matrix(um),matrix(lm)),ncol=length(results.1A$d),byrow=TRUE,dimnames = list(c(),1:length(results.1A$d)))
    row.names(margins)=c(paste("Upper Margin",cl.1A[4]),paste("Upper Margin",cl.1A[3]),paste("Upper Margin",cl.1A[2]),paste("Upper Margin",cl.1A[1]),paste("Lower Margin",cl.1A[1]),paste("Lower Margin",cl.1A[2]),paste("Lower Margin",cl.1A[3]),paste("Lower Margin",cl.1A[4]))
    colnames(margins)=results.1A$d
    margins<-as.matrix(margins)})  
  #output$Table.1A.1 <- renderTable({tab.1A.1()},rownames = TRUE)
  
  #######################  SCENARIO 1B: mydata1B
  mydata.1B<-reactive({
    choice.1B=as.numeric(input$choice.1B)
    range.sample.size=input$range.1B[1]:input$range.1B[2]
    distance.sample.size=input$range.1B[2]-input$range.1B[1]+1
    margin=c(input$u.m,input$l.m)
    if(input$clustering1B=='no'){pen.AB=rep(input$ab.1B,input$pen.1B)} else {pen.AB=extract(input$penAB.1B)}
    withProgress(message = 'Progress time',value = 0,{N<-10
      for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
        sam.mean=matrix(NA,length(choice.1B)*distance.sample.size*input$iteration.1B); dim(sam.mean)=c(length(choice.1B),distance.sample.size,input$iteration.1B)
        for(i in 1:input$iteration.1B){sam.pen=matrix(NA,input$pen.1B*length(choice.1B)*distance.sample.size); dim(sam.pen)=c(input$pen.1B,length(choice.1B),distance.sample.size)
          for(s in 1:distance.sample.size){for(c in 1:length(choice.1B)){pp=sort(sample(1:input$pen.1B,choice.1B[c]))
            for(p in pp){sam.pen[p,c,s]=mean(sample(rnbinom(fish,size=kappa,mu=pen.AB[p]),range.sample.size[s]))}}}
          sam.mean[,,i]=colMeans(sam.pen,na.rm=T)}
        ss=numeric(); range.ab=matrix(NA,length(choice.1B)*distance.sample.size*2); dim(range.ab)=c(length(choice.1B),distance.sample.size,2) 
        for(c in 1:length(choice.1B)){for(s in 1:distance.sample.size){range.ab[c,s,]=sort(sam.mean[c,s,])[c(round((1-input$cl.1B)/2*input$iteration.1B),round((1-(1-input$cl.1B)/2)*input$iteration.1B))]}}
        for(c in 1:length(choice.1B)){ss[c]=range.sample.size[max(which(ifelse(range.ab[c,,1]>mean(pen.AB)-margin[1]&range.ab[c,,2]<mean(pen.AB)+margin[2],1,0)==0))]+1}}})
    ss=ifelse(ss==input$range.1B[2]+1,NA,ss)
    results.1B<-list(a=choice.1B, b=ss,  d=range.ab, f=pen.AB, j=distance.sample.size)})
  
  output$ui.1B.1 <-renderUI({textInput("penAB.1B", "Abundance of Pens", value = NULL, placeholder = 'Use a comma between abundance')})
  output$Text.1B.1<-renderText({ABofPen=extract(input$penAB.1B)
    out<-round(sum(ABofPen)/length(ABofPen),2)
    ifelse(length(ABofPen)<1, "NULL", out)})
  output$Text.1B.2<-renderText({validate(need(input$penAB.1B, ""))
    pen.AB=extract(input$penAB.1B)
    PP=matrix(NA,1000,input$pen.1B)
    for(p in 1:input$pen.1B){PP[,p]=rnbinom(1000,size=2.19,mu=pen.AB[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    rho=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2); rho})
  output$Plot.1B.2<-renderPlot({validate(need(input$penAB.1B, ""))
    pen.AB=extract(input$penAB.1B)
    plot(pen.AB, xlab="Pen number",ylab="Abundance",main="",ylim=c(0,max(pen.AB)),cex=1.5,las=1); abline(h=mean(pen.AB),lty=2)})
  tab.1B.1<-eventReactive(input$goButton.1B, {results.1B<-mydata.1B()
    matrix(as.character(c(results.1B$a,results.1B$b,results.1B$a*results.1B$b)),ncol=3,dimnames = list(c(),c("Number of Sampled Pens","Number of Sampled Fish per Pen","Total Number of Sampled Fish in a Farm")))})
  #output$Table.1B.1 <- renderTable({tab.1B.1()})
  plo.1B.1<-eventReactive(input$goButton.1B, {results.1B<-mydata.1B()
    par(oma=c(0,0,0,0),mar=c(4,2,4,6),cex=1.2,pin=c(4.85,3))
    for(c in 1:length(results.1B$a)){matplot(lwd=2,results.1B$d[c,,],xlim=c(1,results.1B$j),ylim=c(range(results.1B$d[,,])),xlab="Number of sampled fish per pen",ylab="Abundance",col=c,lty=c,type="l",axes=F); par(new=T)}
    axis(1,1:results.1B$j,input$range.1B[1]:input$range.1B[2]); axis(2,las=1); box()
    abline(h=mean(results.1B$f),lty=2); abline(h=mean(results.1B$f)-input$l.m,lty=3); abline(h=mean(results.1B$f)+input$u.m,lty=3)
    legend(lwd=2,"right",border = "black",horiz = FALSE,xpd=TRUE,inset=c(-0.473,0),bty="o",legend=results.1B$a,title="Number of sampled pens",lty=1:length(results.1B$a),col=1:length(results.1B$a))  })
  output$Plot.1B.1 <- renderPlot({plo.1B.1()})
  
  #######################  SCENARIO 2A: mydata2A
  cl.2A=c(0.8,0.9,0.95,0.99) # confidence leve
  mydata.2A<-reactive({choice.2A=as.numeric(input$choice.2A)
    pen.sample.size=round(input$total.sample.size.2A/choice.2A) # sample size per pen
    ab=c(input$ab1.2A,input$ab2.2A)
    if(input$clustering2A=='no'){pen.AB1=rep(input$ab1.2A,input$pen.2A)} else {pen.AB1=extract(input$penAB1.2A)}; if(input$clustering2A=='no'){pen.AB2=rep(input$ab2.2A,input$pen.2A)} else {pen.AB2=extract(input$penAB2.2A)}
    pen.AB=matrix(NA,input$pen.2A,2); pen.AB[,1]=pen.AB1; pen.AB[,2]=pen.AB2
    withProgress(message = 'Progress time',value = 0,{N<-10
      for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
        sam.mean=matrix(NA,length(choice.2A)*input$iteration.2A*2); dim(sam.mean)=c(length(choice.2A),input$iteration.2A,2)
        for(i in 1:input$iteration.2A){for(a in 1:2){sam.pen=matrix(NA,input$pen.2A,length(choice.2A))
        for(c in 1:length(choice.2A)){pp=sort(sample(1:input$pen.2A,choice.2A[c]))
        for(p in pp){sam.pen[p,c]=mean(sample(rnbinom(fish,size=kappa,mu=pen.AB[p,a]),pen.sample.size[c]),na.rm = T)
        sam.mean[,i,a]=colMeans(sam.pen,na.rm=T)}}}}
        ssm=matrix(NA,length(choice.2A)*input$iteration.2A*2); dim(ssm)=c(length(choice.2A),input$iteration.2A,2)
        for(a in 1:2){for(c in 1:length(choice.2A)){ssm[c,,a]=sort(sam.mean[c,,a],na.last = T)}}
        CV=matrix(NA,length(choice.2A),length(cl.2A)); s.power=matrix(NA,length(choice.2A),length(cl.2A))
        for(r in 1:length(cl.2A)){for(c in 1:length(choice.2A)){CV[c,r]=ifelse(ab[1]<ab[2], ssm[c,cl.2A[r]*input$iteration.2A,1], ssm[c,(1-cl.2A[r])*input$iteration.2A,1])
          s.power[c,r]=ifelse(ab[1]<ab[2], mean(ssm[c,cl.2A[r]*input$iteration.2A,1]<ssm[c,,2],na.rm = T), mean(ssm[c,(1-cl.2A[r])*input$iteration.2A,1]>ssm[c,,2],na.rm = T))}}}})
    results.2A<- list(a=s.power, b=round(CV,2), c=cl.2A, d=sam.mean, e=pen.sample.size, f=choice.2A)})
  
  output$ui.2A.1 <-renderUI({textInput("penAB1.2A", "Abundance 1 of Pens", value = NULL, placeholder = 'Use a comma between abundance')})
  output$ui.2A.2 <-renderUI({textInput("penAB2.2A", "Abundance 2 of Pens", value = NULL, placeholder = 'Use a comma between abundance')})
  output$Text.2A.1<-renderText({out<-round(sum(extract(input$penAB1.2A))/length(extract(input$penAB1.2A)),2)
    ifelse(length(extract(input$penAB1.2A))<1, "NULL", out)})
  output$Text.2A.2<-renderText({out<-round(sum(extract(input$penAB2.2A))/length(extract(input$penAB2.2A)),2)
    ifelse(length(extract(input$penAB2.2A))<1, "NULL", out)})
  output$Text.2A.3<-renderText({validate(need(input$penAB1.2A, ""))
    PP=matrix(NA,1000,input$pen.2A)
    for(p in 1:input$pen.2A){PP[,p]=rnbinom(1000,size=2.19,mu=extract(input$penAB1.2A)[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)})
  output$Text.2A.4<-renderText({validate(need(input$penAB2.2A, ""))
    PP=matrix(NA,1000,input$pen.2A); for(p in 1:input$pen.2A){PP[,p]=rnbinom(1000,size=2.19,mu=extract(input$penAB2.2A)[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2]); round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2)})
  output$Plot.2A.1<-renderPlot({validate(need(input$penAB1.2A, ""))
    pen.AB=extract(input$penAB1.2A); plot(pen.AB,main="Abundance1", xlab="Pen number",ylab="Abundance",ylim=c(0,max(pen.AB)),cex=1.5,las=1); abline(h=mean(pen.AB),lty=2)})
  output$Plot.2A.2<-renderPlot({validate(need(input$penAB2.2A, ""))
    pen.AB=extract(input$penAB2.2A); plot(pen.AB,main="Abundance2", xlab="Pen number",ylab="Abundance",ylim=c(0,max(pen.AB)),cex=1.5,las=1); abline(h=mean(pen.AB),lty=2)})
  tab.2A.1<-eventReactive(input$goButton.2A, {results.2A<-mydata.2A()
    df.1=matrix(NA,length(results.2A$f),length(results.2A$c)+1)
    df.1[,1]=results.2A$f; df.1[,-1]=results.2A$c
    colnames(df.1)=c("Number of sampled pens",paste(rep("Critical Value",length(results.2A$c)),results.2A$c)); df.1})  
  #output$Table.2A.1 <- renderTable({tab.2A.1()})
  tab.2A.2<-eventReactive(input$goButton.2A, {results.2A<-mydata.2A()
    df.2=matrix(NA,length(results.2A$f),length(results.2A$c)+1)
    df.2[,1]=results.2A$f; df.2[,-1]=results.2A$b
    colnames(df.2)=c("Number of sampled pens",paste(rep("Critical Value",length(results.2A$c)),results.2A$c)); df.2})  
  #output$Table.2A.2 <- renderTable({tab.2A.2()})
  plo.2A.3<-eventReactive(input$goButton.2A, {results.2A<-mydata.2A()
    par(oma=c(0,0,0,0),mar=c(4,2,4,9),cex=1.2,pin=c(6,2))
    matplot(lwd=2,results.2A$a,type="b",cex.lab=1.5, pch=1:length(results.2A$c),col=1:length(results.2A$c),lty=1:length(results.2A$c),xlab="Number of sampled pens",ylab="Statistical Power",axes=F)
    axis(1,1:length(results.2A$f),results.2A$f); axis(2,las=1,seq(0,1.4,0.2),c(seq(0,1,0.2),NA,NA)); axis(3,1:length(results.2A$f),results.2A$e); box(); 
    mtext(cex=1.5,"Number of sampled fish per pen",3,3)}) #,ylim=c(min(results.2A$a,na.rm = T),1)
  output$Plot.2A.3<-renderPlot({plo.2A.3()})
  plo.2A.4<-eventReactive(input$goButton.2A, {results.2A<-mydata.2A()
    par(oma=c(0,0,0,0),mar=c(4,2,4,9),cex=1.2,pin=c(6,2)) # ,ylim=c(min(results.2A$b,na.rm = T),max(results.2A$b,na.rm = T))
    matplot(lwd=2,results.2A$b,type="b",cex.lab=1.5,pch=1:length(results.2A$c),col=1:length(results.2A$c),lty=1:length(results.2A$c),xlab="Number of sampled pens",ylab="Critical Value",axes=F)
    axis(1,1:length(results.2A$f),results.2A$f); axis(2,las=1); axis(3,1:length(results.2A$f),results.2A$e); box(); mtext(cex=1.5,"Number of sampled fish per pen",3,3)
    legend(lwd=2,"right",border = "black",horiz = FALSE,xpd=TRUE,inset=c(-0.28,0),bty="o",title="Confidence Level",legend=results.2A$c,pch=1:length(results.2A$c),col=1:length(results.2A$c),lty=1:length(results.2A$c))})
  output$Plot.2A.4<-renderPlot({plo.2A.4()})
  
  #######################  SCENARIO 2B: mydata2B
  mydata.2B<-reactive({
    choice.2B=as.numeric(input$choice.2B)
    ab=c(input$ab1.2B,input$ab2.2B)
    range.sample.size=input$range.2B[1]:input$range.2B[2] # number of sampled fish per pen
    distance.sample.size=input$range.2B[2]-input$range.2B[1]+1
    if(input$clustering2B=='no'){pen.AB1=rep(input$ab1.2B,input$pen.2B)} else {pen.AB1=extract(input$penAB1.2B)}; if(input$clustering2B=='no'){pen.AB2=rep(input$ab2.2B,input$pen.2B)} else {pen.AB2=extract(input$penAB2.2B)}
    pen.AB=matrix(NA,input$pen.2B,2); pen.AB[,1]=pen.AB1; pen.AB[,2]=pen.AB2
    withProgress(message = 'Progress time',value = 0,{N<-10
    for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
      sam.mean=matrix(NA,length(choice.2B),distance.sample.size*input$iteration.2B*2); dim(sam.mean)=c(length(choice.2B),distance.sample.size,input$iteration.2B,2)
      for(a in 1:2){for(i in 1:input$iteration.2B){sam.pen=matrix(NA,input$pen.2B*length(choice.2B)*distance.sample.size); dim(sam.pen)=c(input$pen.2B,length(choice.2B),distance.sample.size)
      for(s in 1:distance.sample.size){for(c in 1:length(choice.2B)){pp=sort(sample(1:input$pen.2B,choice.2B[c]))
      for(p in pp){sam.pen[p,c,s]=mean(sample(rnbinom(fish,size=kappa,mu=pen.AB[p,a]),(input$range.2B[1]:input$range.2B[2])[s]))}}}
      sam.mean[,,i,a]=colMeans(sam.pen,na.rm=T)}}
      comp.cl=matrix(NA,length(choice.2B)*distance.sample.size*input$iteration.2B); dim(comp.cl)=c(length(choice.2B),distance.sample.size,input$iteration.2B)
      for(i in 1:input$iteration.2B){for(s in 1:distance.sample.size){for(c in 1:length(choice.2B)){comp.cl[c,s,i]=ifelse(input$ab.standard.2B=="Lower Abundance",sum(sam.mean[c,s,i,1]>sam.mean[c,s,,2]),sum(sam.mean[c,s,i,1]<sam.mean[c,s,,2]))}}}
      prob=matrix(NA,length(choice.2B),input$range.2B[2]-input$range.2B[1]+1); for(c in 1:length(choice.2B)){prob[c,]=1-rowMeans(comp.cl[c,,])/input$iteration.2B}
      ss=numeric(); for(c in 1:length(choice.2B)){sss=ifelse(mean(ifelse(rowMeans(comp.cl[c,,])/input$iteration.2B<(1-input$cl.2B),1,0))==0||mean(ifelse(rowMeans(comp.cl[c,,])/input$iteration.2B<(1-input$cl.2B),1,0))==1,NA,range.sample.size[1]+max(which(ifelse(rowMeans(comp.cl[c,,])/input$iteration.2B<(1-input$cl.2B),1,0)==0)))
      ss[c]=sss}
      for(c in 1:length(choice.2B)){ss[c]=ifelse(prob[c,input$range.2B[2]-input$range.2B[1]+1]<input$cl.2B,NA,ss[c])}
      power=numeric(); critical.value=numeric()
      if(input$ab.standard.2B=="Lower Abundance"){for(c in 1:length(choice.2B)){critical.value[c]=sort(sam.mean[c,ss[c]-input$range.2B[1]+1,,1])[input$iteration.2B*input$cl.2B]
      power[c]=sum(ifelse(sam.mean[c,which((range(input$range.2B)[1]:range(input$range.2B)[2])==ss[c]),,2]>critical.value[c],1,0))/input$iteration.2B}
      } else {
        for(c in 1:length(choice.2B)){critical.value[c]=sort(sam.mean[c,ss[c]-input$range.2B[1]+1,,1])[input$iteration.2B*(1-input$cl.2B)]
        power[c]=1-sum(ifelse(sam.mean[c,which((range(input$range.2B)[1]:range(input$range.2B)[2])==ss[c]),,2]>critical.value[c],1,0))/input$iteration.2B}}
      for(p in 1:length(choice.2B)){power[p]=ifelse(is.na(ss[p]),NA,power[p])}
      for(p in 1:length(choice.2B)){critical.value[p]=round(ifelse(is.na(ss[p]),NA,critical.value[p]),2)}}}) # withProgress
    results.2B<- list(a=ss, b=power, c=critical.value, d=comp.cl, e=distance.sample.size, f=sam.mean, g=choice.2B)}) # mydata2B
  
  output$ui.2B.1 <-renderUI({textInput("penAB1.2B", "Abundance 1 of Pens", value = NULL, placeholder = 'Use a comma between abundance')})
  output$ui.2B.2 <-renderUI({textInput("penAB2.2B", "Abundance 2 of Pens", value = NULL, placeholder = 'Use a comma between abundance')})
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
  tab.2B.1<-eventReactive(input$goButton.2B, {results.2B<-mydata.2B()
  matrix(as.character(c(results.2B$g,results.2B$a,results.2B$g*results.2B$a,round(results.2B$b,2),round(results.2B$c,2))),ncol=5,dimnames = list(c(),c("Number of Sampled Pens","Number of Sampled Fish per Pen","Total Number of Sampled Fish in a Farm","Power","Critical Value")))})
  #output$Table.2B.1<-renderTable({tab.2B.1()})
  plo.2B.3<-eventReactive(input$goButton.2B, {results.2B<-mydata.2B()
  par(oma=c(0,0,0,0),mar=c(4,2,4,9),cex=1.2,pin=c(4.7,3))
  for(c in 1:length(results.2B$g)){plot(lwd=2,1-rowMeans(results.2B$d[c,,])/input$iteration.2B,type="l",col=c,lty=c,xlim=c(1,results.2B$e),xlab="Sampled fish per pen",ylab="Confidence Level",axes = F); 
    par(new=T)}; axis(1,1:results.2B$e,input$range.2B[1]:input$range.2B[2]); axis(2,las=1); box(); par(new=F); abline(h=input$cl.2B,lty=2) 
  legend(lwd=2,"right",border = "black",horiz = FALSE,xpd=TRUE,inset=c(-0.8,0),bty="o",legend=results.2B$g,title="Number of sampled pens",col=1:length(results.2B$g),lty=1:length(results.2B$g))})
  output$Plot.2B.3<-renderPlot({plo.2B.3()})
 
  #######################  SCENARIO 3A: mydata3A
  mydata.3A<-reactive({choice.3A=as.numeric(input$choice.3A)
    pen.sample.size=round(input$total.sample.size.3A/choice.3A) # sample size per pen
    if(input$clustering3A=='no'){pen.AB=rep(input$ab.3A,input$pen.3A)} else {pen.AB=extract(input$penAB.3A)}
    withProgress(message = 'Progress time',value = 0,{N<-10
      for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
        sam.mean=matrix(NA,length(choice.3A),input$iteration.3A)
        for(i in 1:input$iteration.3A){sam.pen=matrix(NA,input$pen.3A,length(choice.3A))
        for(c in 1:length(choice.3A)){pp=sort(sample(1:input$pen.3A,choice.3A[c]))
        for(p in pp){sam.pen[p,c]=mean(sample(rnbinom(fish,size=kappa,mu=pen.AB[p]),pen.sample.size[c]))}}
        sam.mean[,i]=colMeans(sam.pen,na.rm=T)}
        ssm=matrix(NA,length(choice.3A),input$iteration.3A); for(c in 1:length(choice.3A)){ssm[c,]=sort(sam.mean[c,],na.last = T)}
        cutline=numeric(); for(p in 1:length(choice.3A)){cutline[p]=ifelse(mean(pen.AB)>=input$th.3A,mean(ssm[p,]>=input$th.3A),mean(ssm[p,]<input$th.3A))}}})
    results.3A<-list(a=cutline, b=ssm, c=pen.sample.size, d=choice.3A)})
 
  output$ui.3A.1 <-renderUI({textInput("penAB.3A", "Abundance of Pens", value = NULL, placeholder = 'Use a comma between abundance')})
  output$Text.3A.1<-renderText({ABofPen=extract(input$penAB.3A)
    out<-round(sum(ABofPen)/length(ABofPen),2)
    ifelse(length(ABofPen)<1, "NULL", out)})
  output$Text.3A.2<-renderText({validate(need(input$penAB.3A, ""))
    pen.AB=extract(input$penAB.3A)
    PP=matrix(NA,1000,input$pen.3A)
    for(p in 1:input$pen.3A){PP[,p]=rnbinom(1000,size=2.19,mu=pen.AB[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    rho=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2); rho})
  output$Plot.3A.2<-renderPlot({validate(need(input$penAB.3A, ""))
    pen.AB=extract(input$penAB.3A); plot(pen.AB, xlab="Pen number",ylab="Abundance",main="",ylim=c(0,max(pen.AB)),cex=1.5,las=1); abline(h=mean(pen.AB),lty=2)})
  plo.3A.1<-eventReactive(input$goButton.3A,{results.3A<-mydata.3A()
    par(oma=c(0,0,0,0),mar=c(4,2,4,9),cex=1.2,pin=c(5.5,2))
    plot(lwd=1.5,results.3A$a,main="",xlab="Number of sampled pens",ylab="Probability",type="b",axes=F)
    axis(1,1:length(results.3A$d),results.3A$d); axis(2,las=1) ; box(); axis(3,1:length(results.3A$d),results.3A$c); mtext(cex=1.5,"Number of sampled fish per pen",side=3,line=2)})
  output$Plot.3A.1<-renderPlot({plo.3A.1()})
  tab.3A.1<-eventReactive(input$goButton.3A, {results.3A<-mydata.3A()
  matrix(results.3A$a,nrow=1,dimnames = list(c(""),c(paste("Number of sampled pens: ",rep(results.3A$d)))))})  
  #output$Table.3A.1 <- renderTable({h4("Probability of correct decision"); tab.3A.1()})
  
  #######################  SCENARIO 3B: mydata3B
  mydata.3B<-reactive({
    choice.3B=as.numeric(input$choice.3B)
    if(input$clustering3B=='no'){pen.AB=rep(input$ab.3B,input$pen.3B)} else {pen.AB=extract(input$penAB.3B)}
    range.sample.size=input$range.3B[1]:input$range.3B[2]
    
    withProgress(message = 'Progress time',value = 0,{N<-10
      for(n in 1:N){incProgress(1/N,detail=paste("Doing simulation"))
        sam.mean=matrix(NA,length(range.sample.size)*length(choice.3B)*input$iteration.3B); dim(sam.mean)=c(length(range.sample.size),length(choice.3B),input$iteration.3B)
        for(s in 1:length(range.sample.size)){for(i in 1:input$iteration.3B){sam.pen=matrix(NA,input$pen.3B,length(choice.3B))
        for(c in 1:length(choice.3B)){pp=sort(sample(1:input$pen.3B,choice.3B[c]))
        for(p in pp){sam.pen[p,c]=mean(sample(rnbinom(fish,size=kappa,mu=pen.AB[p]),range.sample.size[s]))}}
        sam.mean[s,,i]=colMeans(sam.pen,na.rm=T)}}
        cutline=matrix(NA,length(range.sample.size),length(choice.3B))
        for(p in 1:length(choice.3B)){for(s in 1:length(range.sample.size)){ssm=sort(sam.mean[s,p,])
        cutline[s,p]=ifelse(input$ab.3B>=input$th.3B,mean(ssm>=input$th.3B),mean(ssm<input$th.3B))}}
        min.ss=numeric(); for(p in 1:length(choice.3B)){pass=ifelse(cutline[,p]>=input$cl.3B,1,0)
        min.ss[p]=ifelse(sum(pass)==length(pass),NA,ifelse(sum(pass)==0,NA,input$range.3B[1]+max(which(pass==0))))}}})
    results.3B<-list(a=cutline, b=min.ss, c=choice.3B)})

  output$ui.3B.1 <-renderUI({textInput("penAB.3B", "Abundance of Pens", value = NULL, placeholder = 'Use a comma between abundance')})
  output$Text.3B.1<-renderText({ABofPen=extract(input$penAB.3B)
    out<-round(sum(ABofPen)/length(ABofPen),2)
    ifelse(length(ABofPen)<1, "NULL", out)})
  output$Text.3B.2<-renderText({validate(need(input$penAB.3B, ""))
    pen.AB=extract(input$penAB.3B); PP=matrix(NA,1000,input$pen.3B)
    for(p in 1:input$pen.3B){PP[,p]=rnbinom(1000,size=2.19,mu=pen.AB[p])}
    mPP=melt(PP); colnames(mPP)=c("Order","PenID","Count"); mPP[,2]=as.factor(mPP[,2])
    rho=round(ICCest(PenID, Count, data=mPP, CI.type = "S")$ICC,2); rho})
  output$Plot.3B.2<-renderPlot({validate(need(input$penAB.3B, ""))
    pen.AB=extract(input$penAB.3B); plot(pen.AB, xlab="Pen number",ylab="Abundance",main="",ylim=c(0,max(pen.AB)),cex=1.5,las=1); abline(h=mean(pen.AB),lty=2)})
  tab.3B.1<-eventReactive(input$goButton.3B, {results.3B<-mydata.3B()
    matrix(as.character(results.3B$b),nrow=1,dimnames = list(c(""),c(paste("Number of sampled pens: ",rep(results.3B$c)))))})
  #output$Table.3B.1 <- renderTable({tab.3B.1()})
  plo.3B.1<-eventReactive(input$goButton.3B, {results.3B<-mydata.3B()
    par(oma=c(0,0,0,0),mar=c(4,2,4,9),cex=1.2,pin=c(5.5,3))
    matplot(lwd=2,results.3B$a,type="l",xlim=c(1,(input$range.3B[2]-input$range.3B[1]+1)),ylim=c(min(results.3B$a,na.rm=T),1),xlab="Number of sampled fish per pen",ylab="Probability",axes=F)
    axis(1,1:(input$range.3B[2]-input$range.3B[1]+1),input$range.3B[1]:input$range.3B[2]); axis(2,las=1); box(); abline(h=input$cl.3B, lty=2)
    legend(lwd=2,"right",border = "black",horiz = FALSE,xpd=TRUE,inset=c(-0.6,-0.6),bty="o",legend=results.3B$c,title="Number of sampled pens",lty=1:(length(results.3B$c)),col=1:(length(results.3B$c)))})
  output$Plot.3B.1 <- renderPlot({plo.3B.1()})
  }
