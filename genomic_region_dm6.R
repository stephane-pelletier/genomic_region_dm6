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
		selectInput("chr","Chromosome",choices=c("chr2L","chr2R","chr3L","chr3R","chr4","chrX"),selected=c("chr2L"),multiple=F),
		textOutput("dm6"),
		uiOutput("toCol"),
		uiOutput("toCol1"),
		textOutput("text"),
		dataTableOutput(outputId = "view"))
 	 )
)
server<-function(input, output,session){
	options(shiny.maxRequestSize=30*1024^2)
	output$toCol<-renderUI({
		numericInput("start","Start",value=input$end-25)
 	  })
output$toCol1<-renderUI({
	numericInput("end","End",value=input$start+25)
 	   })
allgenes <-  reactive({
	allGenes <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)})
GRan <-  reactive({
	if (is.na(input$start)==F){
GRan=GRanges(seqnames=as.factor(input$chr),
	ranges=IRanges(start=as.numeric(input$start),
		end=as.numeric(input$end)))
	}else{}})
seqe <-  reactive({
	seqe=getSeq(Dmelanogaster, GRan()[1:length(GRan())])})
output$text<-renderText({
	if (is.na(input$start)==F){
	a<-c(as.vector(seqe()))
	}else{}})


extraCols_narrowPeak <-  reactive({
	extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")})

annoData <-  reactive({
	annoData <- annoGR(allgenes(), feature="gene")})

annotatedPeak <-  reactive({
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


annotatedPeak3 <-  reactive({
	annotatedPeak31 <- as.data.frame(annotatedPeak())
		annotatedPeak31})

output$view <-  renderDataTable({
	if (is.na(input$start)==F){
		annotatedPeak3()
			}else{}
  })

}

shinyApp(ui, server)