# main application file

library(shiny)
library(shinydashboard)
library(fresh)
library(ggplot2)
library(dplyr)
library(pwr)
library(DT)

# source all ui and server modules
source("R/ui_inputs.R")
source("R/ui_outputs.R")
source("R/ui_description.R")
source("R/server_inputs.R")
source("R/server_outputs.R")

# ================================= 
# 1. UI
# ================================= 
ui <- dashboardPage(
  skin = "blue", 
  dashboardHeader(title = "Trial Comparison"),
  dashboardSidebar(
    sidebarMenu(id = "tabs",
                tags$li(class = "header", "TRIAL COMPARISON DESIGN"),
                menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
                menuItem("Github", icon = icon("github"), href = "https://github.com/peiyuliu-biostats/anpp-trial-comparison"),
                menuItem("Model Description", tabName = "model_description", icon = icon("book-open")),
                menuItem("shinyapps.io", icon = icon("cloud-upload-alt"), href = "https://peiyuliu.shinyapps.io/anpp-trial-comparison/"),
                menuItem("Author", tabName = "author", icon = icon("user-circle"))
    )
  ),
  dashboardBody(
    # theme 
    tags$head(tags$style(HTML("
  /* --- General & Header --- */
  .skin-blue .main-header .navbar,
  .skin-blue .main-header .logo {
    background-color: #6c757d !important; /* Top header bar: grey */
  }

  /* --- Sidebar --- */
  .skin-blue .main-sidebar {
    background-color: #f8f9fa !important; /* Sidebar background: light grey */
  }
  .skin-blue .sidebar-menu > li > a {
    color: #333333 !important; /* Sidebar text: black */
  }
  .skin-blue .sidebar-menu > li.header {
     color: #6f42c1 !important; /* Sidebar header text: purple */
     background: #f8f9fa !important;
     font-weight: bold;
  }
  /* Combined rule for active and hovered sidebar items */
  .skin-blue .sidebar-menu > li:hover > a,
  .skin-blue .sidebar-menu > li.active > a {
    background-color: #ffffff !important; /* Hover/active background: white */
  }
  /* Specific style for active item border */
  .skin-blue .sidebar-menu > li.active > a {
    border-left-color: #6f42c1 !important; /* Active item left border: purple */
  }

  /* --- Box Styles --- */
  /* Primary Box (Purple) */
  .skin-blue .box.box-solid.box-primary {
    border-color: #6f42c1 !important;
  }
  .skin-blue .box.box-solid.box-primary > .box-header {
    color: #ffffff !important;
    background-color: #6f42c1 !important;
  }

  /* Combined rules for Info & Success boxes */
  .skin-blue .box.box-solid.box-info,
  .skin-blue .box.box-solid.box-success {
    background-color: #ffffff !important; /* Box body for both: white */
  }
  .skin-blue .box.box-solid.box-info > .box-header,
  .skin-blue .box.box-solid.box-success > .box-header {
    color: #333333 !important; /* Header text for both: black */
  }

  /* Specific rules for Info Box (Light Purple) */
  .skin-blue .box.box-solid.box-info {
    border: 1px solid #d3bce3 !important;
  }
  .skin-blue .box.box-solid.box-info > .box-header {
    background-color: #f2e8f9 !important;
  }

  /* Specific rules for Success Box (Light Orange) */
  .skin-blue .box.box-solid.box-success {
    border: 1px solid #f7d3ba !important;
  }
  .skin-blue .box.box-solid.box-success > .box-header {
    background-color: #fce5d4 !important;
  }
"))),
    
    # assemble the body using the UI modules
    tabItems(
      # main dashboard boody tab, from input and output modules
      tabItem(tabName = "dashboard",
              uimod_content_inputs(),
              uimod_content_outputs()
      ),
      
      # static content tabs, from the description module
      uimod_content_doc()
    )
    
  )
)


# ================= 
# 2. Server
# ================= 
server <- function(input, output, session) {
  
  # 1. inputs (reactive state & core computings)
  reactives <- servermod_inputs_computing(input)
  
  
  # 2. outputs (Rendering)
  servermod_output_render(output, reactives)
  
  # 3. update UI (no needed now)
  
  # 4. UX enhancements (no needed now)
  
  # 5. performance optimization (Handled within modules)
  # `debounce` is used in server_inputs.R to improve performance with numeric inputs.
}

# ============ 
# 3. run  app
# ============= 
shinyApp(ui = ui, server = server)