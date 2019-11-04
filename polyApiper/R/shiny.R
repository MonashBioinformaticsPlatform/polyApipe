
oned_plot <- function(x, value, weight=NULL, cluster, value_label="", weight_label="") {
    show_weight <- !is.null(weight)
    if (is.null(weight))
        weight <- rep(1,length(x))
    
    df <- tibble(x=x,value=value,weight=weight,cluster=cluster) %>%
        arrange(weight, !is.na(value), abs(value), runif(n()))
    
    mid <- df %>%
        group_by(cluster) %>%
        summarize(x=mean(x),value=weighted.mean(value, weight))

    ggplot(df, aes(x=x, y=value, color=factor(cluster))) +
        geom_hline(yintercept=0) +
        geom_point(aes(size=weight)) +
        geom_label(aes(label=cluster), data=mid, show.legend=F) +
        scale_size_area(
            limits=c(0, max(weight)), max_size=log2(max(weight)+1), 
            guide=if (show_weight) "legend" else FALSE) +
        guides(color=FALSE, x=FALSE) +
        scale_x_continuous(breaks=c()) +
        labs(x="",y="", title=value_label, size=weight_label)
}


twod_plot <- function(x, y, value, weight=NULL, cluster, signed, value_label="", weight_label="") {
    show_weight <- !is.null(weight)
    if (is.null(weight))
        weight <- rep(1,length(x))
    
    df <- tibble(x=x,y=y,value=value,weight=weight,cluster=cluster) %>%
        arrange(weight, !is.na(value), abs(value), runif(n()))
    
    mid <- df %>%
        group_by(cluster) %>%
        summarize(
            x=mean(x),
            y=mean(y),
            value=weighted.mean(value, weight))

    plot <- ggplot(df, aes(x=x, y=y)) +
        geom_point(aes(size=weight, color=value)) +
        geom_label(aes(label=cluster), data=mid, show.legend=F) +
        coord_fixed() +
        scale_x_continuous(breaks=c()) +
        scale_y_continuous(breaks=c()) +
        scale_size_area(
            limits=c(0, max(weight)), max_size=log2(max(weight)+1), 
            guide=if (show_weight) "legend" else FALSE) +
        labs(x="",y="",title=value_label,color=value_label,size=weight_label)
    
    if (signed)
        plot <- plot + scale_color_gradient2(high="red",low="blue")
    else
        plot <- plot + scale_color_viridis_c()
    
    plot
}


ui <- function(request) {
    gene_tab <- tabPanel(
        "Gene",
        selectizeInput("gene", label="Gene", choices=c("Select a gene"=""), width="50%"),
        fluidRow(
            column(6,plotOutput("expr_1d")),
            column(6,plotOutput("expr_2d"))),
        fluidRow(
            column(6,plotOutput("shift_1d")),
            column(6,plotOutput("shift_2d"))))

    tabset <- navlistPanel(
        widths=c(2,10), well=FALSE,
        gene_tab)

    ui <- fluidPage(
        titlePanel("polyApipe shiny app"),
        tabset)
    
    ui
}


server <- 
        function(gene_expr, shift, peak_expr, layout_1d, layout_2d, clusters) 
        function(input, output, session) {    

    choices <- union(rownames(gene_expr), rownames(shift))
    names(choices) <- choices
    names(choices)[match(rownames(gene_expr), choices)] <- paste(rownames(gene_expr), rowData(gene_expr)$symbol)
    names(choices)[match(rownames(shift), choices)] <- paste(rownames(shift), rowData(shift)$symbol)
    
    stopifnot(identical(colnames(gene_expr), colnames(shift)))
    stopifnot(identical(colnames(gene_expr), colnames(peak_expr)))
   
    updateSelectizeInput(session, "gene", choices=c("Select a gene"="",choices), server=T)
    
    
    output$expr_1d <- renderPlot({
        req(input$gene %in% rownames(gene_expr))
        expr <- logcounts(gene_expr)[input$gene,]
        oned_plot(layout_1d, expr, cluster=clusters,
            value_label="log2 count")            
    })

    output$expr_2d <- renderPlot({
        req(input$gene %in% rownames(gene_expr))
        expr <- logcounts(gene_expr)[input$gene,]
        expr[expr==0] <- NA
        twod_plot(layout_2d[,1],layout_2d[,2], expr, cluster=clusters, signed=FALSE,
            value_label="log2 count")            
    })

    output$shift_1d <- renderPlot({
        req(input$gene %in% rownames(shift))
        value <- weitrix_x(shift)[input$gene,]
        weight <- weitrix_weights(shift)[input$gene,]
        oned_plot(layout_1d, value, weight=weight, cluster=clusters,
            value_label="shift", weight_label="UMI count")            
    })

    output$shift_2d <- renderPlot({
        req(input$gene %in% rownames(shift))
        value <- weitrix_x(shift)[input$gene,]
        weight <- weitrix_weights(shift)[input$gene,]
        twod_plot(layout_2d[,1], layout_2d[,2], value, weight=weight, cluster=clusters, signed=TRUE,
            value_label="shift", weight_label="UMI count")            
    })
}



#' @export
polyApipe_shiny <- function(
        gene_expr,
        shift,
        peak_expr,        
        layout_1d=NULL,
        layout_2d=NULL,
        clusters=NULL) {
    shinyApp(
        ui, 
        server(gene_expr=gene_expr, shift=shift, peak_expr=peak_expr, 
            layout_1d=layout_1d, layout_2d=layout_2d, clusters=clusters))
}



