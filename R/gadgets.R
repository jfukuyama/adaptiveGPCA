#' Shiny gadget for adaptive gPCA
#'
#' Shiny gadget that shows the ordinations from an entire family of
#' gPCAs and returns a gPCA object with the one selected
#'
#' @param fullFamily The output from gpcaFullFamily
#' @param sampleData Optional data used for plotting the samples
#' @param sample_mapping An aesthetic mapping to be passed to ggplot
#' for plotting the samples
#' @param varData Optional data used for plotting the variables
#' @param var_mapping An aesthetic mapping to be passed to ggplot for
#' plotting the variables
#' @param layout A vector of length 2. The first number gives the
#' number of columns (out of 12) for the sidebar, the second number
#' gives the number of columns (out of 12) for the sample plot in the
#' main panel.
#' @import shiny
#' @export
visualizeFullFamily <-
    function(fullFamily, sampleData = NULL, sample_mapping = aes(x = Axis1, y = Axis2),
             sample_facet = NULL, 
             varData = NULL, var_mapping = aes(x = Axis1, y = Axis2), layout = c(2,6)) {

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
                if(!is.null(sampleData))
                    p = ggplot(data.frame(fullFamily$locations[[input$r]], sampleData),
                        sample_mapping) + geom_point()
                else
                    p = ggplot(data.frame(fullFamily$locations[[input$r]]), sample_mapping) +
                        geom_point()
                if(!is.null(sample_facet))
                    p = p + sample_facet
                p
            })
            output$plot_species = renderPlot({
                if(!is.null(varData))
                    p = ggplot(data.frame(fullFamily$species[[input$r]], varData),
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


processFamily <- function(out.ff, rvec = (0:100)/100) {
    locationsfull = data.frame(Reduce(rbind, out.ff$locations))
    speciesfull = data.frame(Reduce(rbind, out.ff$species))
    locationsfull$r = rep(rvec, each = nrow(out.ff$locations[[1]]))
    speciesfull$r = rep(rvec, each = nrow(out.ff$species[[1]]))
    p1 = ggplot() +
        geom_point(aes(x = Axis1, y = Axis2, showSelected = r), data = locationsfull)
    p2 = ggplot() +
        geom_point(aes(x = Axis1, y = Axis2, showSelected = r), data = speciesfull)
    selecting = ggplot() +
        make_tallrect(locationsfull, "r") + 
        geom_point(aes(x = r, y = 1, clickSelects = r), data = locationsfull)

    animint2dir(list(p1 = p1, p2 = p2, selecting = selecting,
                     time = list(variable = "r", ms = 300)),
                out.dir = "~/simple", open.browser = FALSE)
    servr::httd("~/simple")
}
