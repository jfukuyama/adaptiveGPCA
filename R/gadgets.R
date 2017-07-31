#' Shiny gadget for adaptive gPCA
#'
#' Shiny gadget that shows the ordinations from an entire family of
#' gPCAs and returns a gPCA object with the one selected by the user.
#'
#' @param fullFamily The output from \code{\link{gpcaFullFamily}}
#' @param sample_data Optional data used for plotting the samples
#' @param sample_mapping An aesthetic mapping to be passed to
#' \code{\link[ggplot2]{ggplot}} for plotting the samples
#' @param sample_facet A \code{\link[ggplot2]{ggplot}} faceting
#' command used for faceting the samples.
#' @param var_data Optional data used for plotting the variables
#' @param var_mapping An aesthetic mapping to be passed to
#' \code{\link[ggplot2]{ggplot}} for plotting the variables
#' @param layout A vector of length 2. The first number gives the
#' number of columns (out of 12) for the sidebar, the second number
#' gives the number of columns (out of 12) for the sample plot in the
#' main panel.
#' @return This function will open a 'shiny' app in a browser
#' window. You can investigate the results for different values of
#' \eqn{r} with this app. Once you press the 'done' button, the app
#' will close and the function will return an R object containing the
#' results for the value of \eqn{r} (the regularization parameter)
#' that was chosen in the app. The returned object is a list
#' containing the variable loadings on the principal axes (\code{QV}),
#' the sample/row scores (\code{U}), and the fraction of the variance
#' explained by each of the axes (\code{vars}).
#' @import shiny
#' @importFrom ggplot2 aes ggplot geom_point aes_string
#' @examples
#' \dontrun{
#' data(AntibioticPhyloseq)
#' pp = processPhyloseq(AntibioticPhyloseq)
#' out.ff = gpcaFullFamily(pp$X, Q = pp$Q, D = pp$D, k = 2)
#' out.agpca = visualizeFullFamily(out.ff,
#'     sample_data = sample_data(AntibioticPhyloseq),
#'     sample_mapping = aes(x = Axis1, y = Axis2, color = condition),
#'     var_data = tax_table(AntibioticPhyloseq),
#'     var_mapping = aes(x = Axis1, y = Axis2, color = Phylum))
#' }
#' @export
visualizeFullFamily <-
    function(fullFamily, sample_data = NULL,
             sample_mapping = aes_string(x = "Axis1", y = "Axis2"),
             sample_facet = NULL, 
             var_data = NULL,
             var_mapping = aes_string(x = "Axis1", y = "Axis2"), layout = c(2,6)) {

        ui <- fluidPage(
            headerPanel("Visualization of adaptive gPCA"),
            sidebarPanel(width = layout[1],
                textOutput("rval"),
                sliderInput("r", "Select r index", value = 1, min = 1,
                            max = length(fullFamily[[1]]), step = 1),
                actionButton("down", "", icon = icon("chevron-left")),
                actionButton("up", "", icon = icon("chevron-right")),
                hr(),
                actionButton("done", "Done")),
            mainPanel(width = 12 - layout[1],
                fluidPage(
                    fluidRow(
                        column(layout[2], plotOutput("plot_samples")),
                        column(12 - layout[2], plotOutput("plot_species")))
                )
            )
        )
        

        server <- function(input, output, session) {
            output$rval = renderText({
                paste("r:", names(fullFamily$locations)[input$r])
            })
            observeEvent(input$up, {
                updateSliderInput(session, "r", value = input$r + 1)
            })
            observeEvent(input$down, {
                updateSliderInput(session, "r", value = input$r - 1)
            })

            output$plot_samples = renderPlot({
                if(!is.null(sample_data))
                    p = ggplot(data.frame(fullFamily$locations[[input$r]], sample_data),
                        sample_mapping) + geom_point()
                else
                    p = ggplot(data.frame(fullFamily$locations[[input$r]]), sample_mapping) +
                        geom_point()
                if(!is.null(sample_facet))
                    p = p + sample_facet
                p
            })
            output$plot_species = renderPlot({
                if(!is.null(var_data))
                    p = ggplot(data.frame(fullFamily$species[[input$r]], var_data),
                        var_mapping) + geom_point()
                else
                    p = ggplot(data.frame(fullFamily$species[[input$r]]), var_mapping) +
                        geom_point()
                p
            })


            observeEvent(input$done, {
                r = as.numeric(names(fullFamily$locations)[input$r])
                Qeig = fullFamily$Qeig
                evals = (rep(1 / (1 - r), ncol(Qeig$vectors)) +
                             r^(-1) * Qeig$values^(-1))^(-1)
                out.gpca = gpcaEvecs(fullFamily$X, k = ncol(fullFamily$locations[[1]]),
                    evecs = Qeig$vectors, evals = evals)
                out = list(V = out.gpca$V, U = out.gpca$U, QV = out.gpca$QV,
                           lambda = out.gpca$lambda, vars = out.gpca$vars,
                           r = r, evals = evals)
                class(out) = "adaptivegpca"
                stopApp(out)
            })
        }

  runGadget(ui, server)
}


#' Shiny gadget for tree/taxonomy inspection
#'
#' Shiny gadget that allows users to visualize the scores of the taxa
#' on the agpca axes, their positions on the phylogenetic tree, and
#' their taxonomic assignments.
#'
#' @param agpcafit An agpca object, created either by the function
#'     \code{\link{adaptivegpca}} or by
#'     \code{\link{visualizeFullFamily}}.
#' @param physeq A phyloseq object with a tree and a taxonomy table. 
#' @param axes The axes to plot, must be a vector of two whole
#'     numbers.
#' @param br.length Plot the tree with the branch lengths?
#' @param height The height, in pixels, of the plotting region.
#' @return The function will open a browser window showing the tree
#'     and the locations of the taxa on the selected agpca
#'     axes. "Brushing" over the plot will highlight the positions of
#'     the selected taxa on the tree and list their taxonomic
#'     assignments. Clicking the "done" button will exit the app and
#'     return a data frame containing the positions of the selected
#'     taxa on the agpca axes, the taxonomic assignments of the
#'     selected taxa, and their names.
#' @import shiny
#' @importFrom ggplot2 ggplot geom_point aes_string
#' @importFrom phyloseq plot_tree tax_table phy_tree
#' @examples
#' \dontrun{
#' data(AntibioticPhyloseq)
#' pp = processPhyloseq(AntibioticPhyloseq)
#' out.agpca = adaptivegpca(pp$X, pp$Q, k = 2)
#' treeInspect(out.agpca, AntibioticPhyloseq)
#' }
#' @export
inspectTaxonomy <- function(agpcafit, physeq,
                        axes = c(1,2), br.length = FALSE, height = 600) {
    check_axes(axes, agpcafit)
    check_phyloseq(agpcafit, physeq)
    axis.names = paste("Axis", axes, sep = "")
    axis.labels = paste("Axis ", axes, ": ", round(agpcafit$vars[axes] * 100, digits = 1), "%", sep = "")
    tt2 = data.frame(as(tax_table(physeq), "matrix"))
    ggdf = data.frame(agpcafit$QV, tt2)
    ggdf$otu = rownames(ggdf)
    tr2 = phy_tree(physeq)
    if(!br.length & !is.null(tr2$edge.length)) {
        tr2$edge.length = rep(1, length(tr2$edge.length))
    }
    ptree = plot_tree(tr2)

    ui = fluidPage(
        fluidRow(
            column(7,
                   h1("Taxonomy Inspection")),
            column(1,
                   div(style="float:right;margin-top:25px;",
                       actionButton("done", "Done")))),
        fluidRow(
            column(2,
                   fluidRow(plotOutput("treeplot", height = height))),
            column(6, 
                   fluidRow(plotOutput("plot1", height = height, brush = brushOpts(id = "plot1_brush"))))),
        fluidRow(
            column(8,
                   h3("Selected taxa"))),
        fluidRow(
            column(8,
                   tableOutput("brush_info")))
    )


    server = function(input, output) {
        output$plot1 = renderPlot({
            ggplot(ggdf, aes_string(axis.names[1], axis.names[2])) + geom_point() +
                xlab(axis.labels[1]) + ylab(axis.labels[2])
        })
        otus = reactive({
          brushedPoints(ggdf, input$plot1_brush)$otu  
        })
        output$treeplot = renderPlot({
            ptree +
                geom_point(aes_string(x = "xleft + 1", y = "y"), size = 3,
                           color = "#F98400", data = ptree$data[ptree$data$OTU %in% otus(),])
        })
        output$brush_info = renderTable({
            brushedPoints(ggdf, input$plot1_brush)[,-c(1,2,ncol(ggdf))]
        }, striped = TRUE)
        observe({
            if(input$done > 0){
                stopApp(brushedPoints(ggdf, input$plot1_brush))
            }
        })
    }

    runGadget(ui, server)
}


#' Check compatibility of agpca and phyloseq objects
#'
#' Check that the dimensions of the agpca object match the phyloseq
#' object and that the phyloseq object has a taxonomy table and a
#' phylogenetic tree.
#'
#' @param agpcafit An adaptivegpca object.
#' @param physeq A phyloseq object.
#' @importFrom phyloseq ntaxa access
#' @keywords internal
check_phyloseq <- function(agpcafit, physeq) {
    if(!inherits(agpcafit, "adaptivegpca")) {
        stop("agpcafit must be an object of class adaptivegpca")
    }
    if(!inherits(physeq, "phyloseq")) {
        stop("physeq must be an object of class phyloseq")
    }
    if(nrow(agpcafit$QV) != ntaxa(physeq)) {
        stop("Number of variables in agpcafit not equal to number of taxa in physeq")
    }
    if(is.null(access(physeq, "tax_table", errorIfNULL = FALSE))) {
        stop("physeq must have a tax_table element")
    }
    if(is.null(access(physeq, "phy_tree", errorIfNULL = FALSE))) {
        stop("physeq must have a phy_tree element")
    }
    return(TRUE)
    
}
