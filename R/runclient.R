#' Client for Federated IRT Model Estimation
#' @note This shiny app should be used together with server version. Run server before run client, and get the correct address from server interface to initialize the estimating process.
#'
#' @description Initializes a client interface for the federated learning estimation of Item Response Theory (IRT) model parameters, connecting to a central server to participate in collaborative parameter estimation. It is essential to start the server prior to the client to ensure the client can establish a successful connection, otherwise an error will occur.
#'
#' @details The client interface, built with Shiny, provides an interactive platform that enables users to upload response matrix data in CSV format, connect to a central server, and receive the estimation results once the computation is complete. The client sends computed local statistics or partial results to the server, which then aggregates information from all clients to update the global IRT model parameters. Users can input the server's IP address and port number, reconnect if needed, and visualize the computed item and ability parameters through plots and tables displayed in the interface.
#'
#' The client is capable of uploading data, processing it locally to compute log-likelihood or gradient information, and sending these details to the server based on HTTP POST requests. The client also includes functionality to handle responses from the server, either to signal the status of the connection or to receive and display results of the federated estimation process. Through this interactive client-server architecture, the federated IRT model estimation becomes a seamless process, allowing participants to contribute computational resources while preserving data privacy within their local environments.
#'
#' Additional client functions include local IP retrieval for network communication, server connection initiation, response data processing, and result visualization. Interactive components built in Shiny enable a smooth user experience and real-time updates, making the client an integral part of the federated IRT model estimation framework.
#'
#' @return shows the discriminations and difficulties of each item and plot them. Also displays each students' abilities.
#'

#' @importFrom purrr map
#' @importFrom pracma quadl
#' @importFrom stats optim
#' @importFrom shiny shinyApp fluidPage titlePanel sidebarLayout sidebarPanel mainPanel textInput actionButton verbatimTextOutput fileInput dataTableOutput plotOutput reactiveValues renderPrint
#' @importFrom httr POST content
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal labs
#' @importFrom shinyjs useShinyjs
#' @importFrom callr r_bg r
#' @importFrom shiny runApp
#' @importFrom DT datatable
#'
#' @export

runclient = function(){

  appDir <- system.file("shiny", package = "FedIRT")
  clientAppFile <- file.path(appDir, "client.R")

  if (file.exists(clientAppFile)) {
    shiny::runApp(clientAppFile)
  } else {
    stop("Unable to launch the client app - file client.R not found in the FedIRT package")
  }
}
