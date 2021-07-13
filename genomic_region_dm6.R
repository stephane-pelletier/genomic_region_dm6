library("shiny")
library(shinyWidgets)
library(DT)
library("Biostrings")  
library("BSgenome.Dmelanogaster.UCSC.dm6") 
library("GenomicRanges") 
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library ("ChIPpeakAnno")
library("org.Dm.eg.db")



ui <- 
fluidPage(
    fluidRow(
column(6,
h2("Genomic Range to Sequence")),
column(6,
h2("Sequence to Genomic Range")),
    	column(6, 
		selectInput("chr","Chromosome",choices=c("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY","chrM"),selected=c("chr2L"),multiple=F),
		#textOutput("dm6"),
		uiOutput("toCol"),
		uiOutput("toCol1")),
column(6,
       selectInput("chrnew","Chromosome",choices=c("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY","chrM"),selected=c("chr2L"),multiple=F),
       uiOutput("templateSequence"),
       br(uiOutput("reverse2PCR"))),
column(6,
       actionButton("submitGR2seq", "GR2Seq")),
column(6,
       actionButton("submitSeq2GR", "Seq2GR")),
column(6,
		br(verbatimTextOutput("text")),
		br(dataTableOutput(outputId = "view")),
		br(dataTableOutput(outputId = "view2")),
br(textOutput("runPCR")))

 	 )
)
server<-function(input, output,session){
	options(shiny.maxRequestSize=30*1024^2)

##GR2seq
	output$toCol<-renderUI({
		numericInput("start","Start",value=0)
 	  })
output$toCol1<-renderUI({
	numericInput("end","End",value=0)
 	   })
allgenes <-  eventReactive(input$submitGR2seq,{
	allGenes <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)})
GRan <-  reactive({
	if (is.na(input$start)==F){
GRan=GRanges(seqnames=as.factor(input$chr),
	ranges=IRanges(start=as.numeric(input$start),
		end=as.numeric(input$end)))
	}else{}})
seqe <-  eventReactive(input$submitGR2seq,{
	seqe=getSeq(Dmelanogaster, GRan()[1:length(GRan())])})
output$text<-renderText({
	if (is.na(input$start)==F){
	a<-c(as.vector(seqe()))
	}else{}})


extraCols_narrowPeak <-  eventReactive(input$submitGR2seq,{
	extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")})

annoData <-  eventReactive(input$submitGR2seq,{
	annoData <- annoGR(allgenes(), feature="gene")})

annotatedPeak <-  eventReactive(input$submitGR2seq,{
	annotatedPeak <- annotatePeakInBatch(GRan(), AnnotationData=annoData(),output="both")


annotatedPeak$symbol <- mapIds(org.Dm.eg.db, 
                     keys=annotatedPeak$"feature", 
                     column="SYMBOL", 
                     keytype="FLYBASE",
                     multiVals="first")


annotatedPeak$entrez <- mapIds(org.Dm.eg.db, 
                            keys=annotatedPeak$"feature", 
                            column="ENTREZID", 
                            keytype="FLYBASE",
                            multiVals="first")


annotatedPeak$name2 =   mapIds(org.Dm.eg.db,
                           keys=annotatedPeak$"feature", 
                           column="GENENAME",
                           keytype="FLYBASE",
                           multiVals="first")
annotatedPeak})


annotatedPeak3 <-  eventReactive(input$submitGR2seq,{
	annotatedPeak31 <- as.data.frame(annotatedPeak())
		annotatedPeak31})

output$view <-  renderDataTable({
	if (is.na(input$start)==F){
		annotatedPeak3()
			}else{}
  })

##Seq2GR

output$templateSequence<-renderUI({
	textInput("template", "Template")})

output$reverse2PCR<-renderUI({
	checkboxInput("template2", "Reverse complement", FALSE)})


template2use <-  eventReactive(input$submitSeq2GR,{
temp<-toupper(input$template)
temp<-DNAString(temp)
if (input$template2==TRUE){
temp <- reverseComplement(temp)}else{
temp<-DNAString(input$template)}})

mainChromosomes <-  eventReactive(input$submitSeq2GR,{
mainChromosomes <- c("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY","chrM")})
mainChrSeq <-  eventReactive(input$submitSeq2GR,{
mainChrSeq <- lapply(mainChromosomes(),function(x)BSgenome.Dmelanogaster.UCSC.dm6[[x]])
names(mainChrSeq) <- mainChromosomes()
mainChrSeq})

mainChrSeqSet <-  eventReactive(input$submitSeq2GR,{
mainChrSeqSet <- DNAStringSet(mainChrSeq())
mainChrSeqSet})

GRprocess <-  eventReactive(input$submitSeq2GR,{
fullseq<-mainChrSeqSet()[names(mainChrSeqSet())%in%input$chrnew]
fullseq})
GRprocess2 <-  eventReactive(input$submitSeq2GR,{
    mseqdna<-vmatchPattern(template2use(),GRprocess(),fixed=T,max.mismatch=0)
mseqdna})
matSeq <-  eventReactive(input$submitSeq2GR,{
mseqdna2<-as.matrix(as.data.frame(GRprocess2()))
#mseqdna2<-as.data.frame(GRprocess2())
mseqdna2})



output$view2 <-  renderDataTable({
matSeq()
  })

output$runPCR <- renderText({
})

}

shinyApp(ui, server)
