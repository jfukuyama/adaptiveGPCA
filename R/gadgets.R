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
                out = gpcaEvecs(fullFamily$X, k = ncol(fullFamily$locations[[1]]),
                    evecs = Qeig$vectors, evals = evals)
                out$r = r
                stopApp(out)
            })
        }

  runGadget(ui, server)
}
