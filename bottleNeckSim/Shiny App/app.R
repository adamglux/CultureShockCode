library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bslib)

# Load real datasets
# popSim_without_sorted <- readRDS("reg_sorted.RDS")
# popSim_without_unsorted <- readRDS("reg_unsorted.RDS")
# popSim_with_sorted <- readRDS("prev_sorted.RDS")
# popSim_with_unsorted <- readRDS("prev_unsorted.RDS")
# 
# popSim_with_unsorted <- popSim_with_unsorted %>%
#   filter(`80` != max(`80`))
# 
# popSim_with_sorted <- popSim_with_sorted %>%
#   filter(`80` != max(`80`))

popSim_withPrev <- readRDS("withPrev_data.RDS")
popSim_withoutPrev <- readRDS("withoutPrev_data.RDS")
##

## FG loss data 
FG_losswithPrev <- readRDS("AllFG_withprev.RDS")
FG_losswithoutPrev <- readRDS("AllFG_noprev.RDS")



##################
##################
## simulation priors
##################

pop_data <- tibble::tibble(
  Generation = 75:82,
  Year = c("1830","1855","1880","1905","1930","1955","1980","2005"),
  Phase = c("Before whaling", "Whaling starts", "Whaling continues",
            "Whaling continues", "Bottleneck", "Recovery starts",
            "Hopeful recovery", "Current population"),
  Distribution = c("30,000 – 40,000",
                 "100 – 5,000",
                 "36 – 750",
                 "36 – 750",
                 "36 – 750",
                 "36 – 750",
                 "36 – 750",
                 "500 – 2,500")
)



ui <- fluidPage(
  titlePanel("Population Simulation Viewer"),
  
  
  
  sidebarLayout(
    sidebarPanel(

      radioButtons("fg_prevalence", "Include FG Prevalence?",
                   choices = c("Equally Prevelant Foraging Ground" = "without",
                               "Dominant Foraging Ground" = "with")),
      helpText("Choose whether or not the simulations include FG prevalence effects. Note: these are two distinct datasets from two different runs of simulations. Manipulations below are subsets of these data."),


      # radioButtons("fg_sort", "Foraging Grounds:",
      #              choices = c("Sorted by FG size in Gen 82" = "sorted",
      #                          "Unsorted by FG size in Gen 82" = "unsorted")),
      # helpText("Select whether foraging grounds are sorted by size or left unsorted."),

      radioButtons("success_type", "Success based on:",
                   choices = c("Haplotype Diversity" = "hapdiv",
                               "Foraging Ground Loss" = "fg_loss")),
      helpText("Choose how to define successful simulations. Diversity takes simulations with haplotype diversity in gen 82 that is between the upper and lower limits (set below), and with number of haplotypes set below. FG loss is measured by sims with at least one site that is <1% of the sample population."),

      
      ### slider controls
      
      sliderInput("hap", "Lower and Upper Number of Haplotypes:",
                  min = 0, max = 20,
                  value = c(7,13)),
      helpText("Sets the lower and upper limits to the number of haplotypes in gen 82"),
      
      
      sliderInput("sd", "Set the 'success' boundaries within the standard error of the diversity given below:",
                  min = 0, max = 3, step = 0.25,
                  value = 2),
      helpText("Sets the parameters of 'success' for haplotype diversity based on the standard deviation of the set of simulations."),
      
      ## show loss plot 
      #checkboxInput("checkbox", "Show FG-losses Plot", FALSE), 
      #verbatimTextOutput("value"),
      #helpText("Choose whether to show foraging ground losses plot (may take a few moments to load on some computers)"),
      
      ## plot type 
      
      radioButtons("plot_type", "Plot scale:",
                   choices = c("Mean Population" = "raw",
                               "Log Mean Population" = "log")),
      helpText("Choose whether to plot raw or log-transformed mean population size."),
      
      ## set seed 
      numericInput( 
        "numeric", 
        "set seed", 
        value = 123, 
        min = 1, 
        max = 1e7 
      ),
      helpText("Random examples are extracted using the seed provided here.")

    ),
    
    
    navset_card_underline(
      
      nav_panel("About", 
                p("This is a tool built for manipulation and understanding of the data provided by simulations carried out using simuPOP. In the control-panel to the left, users can select whether to explore simulations carried out with or without initial foraging ground prevalence; and whether to set ‘success’ criteria based on haplotype diversity or by total foraging ground loss."),
                p("Numbers of haplotypes can also be reset to change the number of successful outcomes accordingly. Adjusting the standard error allows users to narrow the acceptance range based on the variation of the empirical diversity (see ‘success’measures below)."),
                p("Changing the plot scale allows users to see a log-transformed view in some plots. And finally, the ‘set seed’ option will allow users to change the random seed governing the random selection of instances in the Random examples tab. "),
                p("The tabs above will allow users to explore summary tables and plots as well as random examples of cases based on their ‘success’ criteria provide."),
                p("Under the ‘Summaries and Tables’ tab, is a summary table for the entire dataset (table 1.), showing n-number of instances, average number of haplotypes, average diversity, average foraging ground frequency, and proportion of the total. The next table (table 2.) shows a breakdown of foraging ground losses by ‘success’ criteria. Average population table shows the average population sizes by generation. This is further illustrated in the plot below. Population densities show population size distributions across all simulations. The average haplotype plot shows declining haplotype numbers on average across all generations of interest. Violin plots show distributions of haplotype diversity, and also foraging ground loss. The ‘FG Loss’ plot shows the number of individuals lost in the smallest feeding ground, for each simulation."),
                
                h3("'Success' measures"),
                p("Empirical diversity is based on genetic monitoring at the key wintering calving and socialising ground of Port Ross Maungahuka Auckland Islands. This sample includes a total haplotype number of 10 unique haplotypes, with a diversity of 0.687, using Nei (1987) pg. 180, and variance of 0.00135. We are able to adjust the 'success' threshold to within a set of upper and lower limits based on the empirical standard error, as well as by number of haplotypes in the sample."),
                textOutput("sd_limits"),
                h3("Population Size Distributions for Simulations"),
                h4("(based on Jackson et al. 2016)"),
                tableOutput("popTable"),
                p("The table above shows the prior distributions of population sizes used in the whale simulation. At generation 76 the population enters into a bottleneck, but is allowed to grow to 5,000 individuals. However, from generation 77 through 81, the population is forced into a tighter bottleneck, for which the parameters stay consistant. The simulation takes a random sample from the upper and lower boundary provided."),
    
                
                ###### references
                ######
                
                h3("References"),
                tags$ul(
                  tags$li("Nei, M. (1989). ",
                          em("Molecular Evolutionary Genetics"),
                          ". Columbia University Press. ISBN: 978-0-231-06321-0."),
                  
                  tags$li("Jackson, J. A., Carroll, E. L., Smith, T. D., Zerbini, A. N., Patenaude, N. J., & Baker, C. S. (2016). 
             An integrated approach to historical population assessment of the great whales: case of the New Zealand southern right whale. ",
                          em("Royal Society Open Science, 3"),
                          ", 150669. https://doi.org/10.1098/rsos.150669"),
                  
                  tags$li("Peng, B., & Kimmel, M. (2005). ",
                          em("simuPOP: a forward-time population genetics simulation environment"),
                          ". ",
                          em("Bioinformatics, 21(18)"),
                          ", 3686–3687. https://doi.org/10.1093/bioinformatics/bti584")
                )
      ),
      
      nav_panel("Summaries and Tables",
                #textOutput("sd_limits"),
                h3("Table 1. Simulations Summary Table"),
                tableOutput("summary_table"),
                tags$p("Table 1. depicts total number of simulations categorised into groups as “success” and “non-success” based on the criteria selected to the left. the average numbers of haplotypes, diversity, and foraging ground frequency are calculated within those groups. The proportion shows the amount of the total simulations belong to each group."),
                h3("Table 2. All Foraging Ground Declines to <1% by Success-group"),
                tableOutput("fg_table"),
                tags$p("Table 2. shows FG loss information based on empirical-data groupings. Subgroups of non-success are based on diversity thresholds (see above). Number of FG (n.FG) shows number of occurences, which may not sum to the total (n) if multiple occurences of loss (<= 1%) happen simutenously."),
                h3("Table 3. Avg. Population sizes (rounded to nearest whole num)"),
                tableOutput("pop_table"),
                
                h3("Table 4. Avg Haplotype Diversity by Year"),
                tableOutput("div_table"),
                
                # h2("Avg. Population sizes by generation for FG loss"),
                # tableOutput("pop_table2"),
            
                h3("Avg. Population Size by Generation"),
                plotOutput("pop_plot"),
                textOutput("plot_caption"),
                h2("Population Densities"),
                plotOutput("pop_density_plot"),
                h2("Avg Haplotypes Diversity of by Year"),
                plotOutput("div_var_plot"),
                
                
                
                
                h2("Avg Number of Haplotypes by Generation"),
                plotOutput("happlot"),
                h2("Average Haplotype Diversity"),
                plotOutput("hapdiv_violin"),
                h2("Foraging Ground Losses"),
                plotOutput("violin"),
                
                ## fg losses
                h2("FG Loss"),
                plotOutput("fg_loss_plot")
                
                
                ),
      
      nav_panel(
        "Random Examples",
        
        h2("Random Sample of Data"),
        tableOutput("sample_table"),
        
        h2("Number of Haplotypes"),
        plotOutput("sample_haps"), 
        
        h2("Foraging Ground trajectories"),
        plotOutput("sample_FGs")
        
      )
    )
  )
  
)

server <- function(input, output, session) {

  ##################
  ##################
  ## capture data
  ##################
  
  selected_data <- reactive({
    if (input$fg_prevalence == "without") {
      data <- popSim_withoutPrev
    } else {
      data <- popSim_withPrev
    }
    
    data <- data %>%
      filter(!is.na(HapDiv)) 
    
    if (input$success_type == "hapdiv") {
      data <- data %>%
        mutate(Group = if_else(HapDiv >= hapdiv_range()["lower"] & HapDiv <= hapdiv_range()["upper"] & hapCur >= min(input$hap) & hapCur <= max(input$hap), 
                               "success", "non-success")) 
      
    } else if (input$success_type == "fg_loss") {
        data <- data %>%
          mutate(
            Group = if_else(
              ratio <= 0.01,
              #ratioFG1 <= 0.01 | ratioFG2 <= 0.01 | ratioFG3 <= 0.01,
              "success", "non-success"))
    }
    
    data$Group <- factor(data$Group, levels = c("success", "non-success"))
    data <- data %>%
      filter(!is.na(Group)) %>%
      mutate(sim_id = row_number())  
      
  })
  
  
  ##################
  ##################
  ## set success threshold
  ##################
  
  hapdiv_range <- reactive({
    mean_hapdiv <- 0.6866629
    sd_hapdiv <- sqrt(0.001348966)
    lower <- round(mean_hapdiv - (input$sd * sd_hapdiv), 3)
    upper <- round(mean_hapdiv + (input$sd * sd_hapdiv), 3)
    c(lower = lower, upper = upper)
  })

  pop_summary <- reactive({
    
    sim_data <- selected_data()

    pop_long <- sim_data %>%
      select(sim_id, Group, `76`:`82`) %>%
      pivot_longer(
        cols = `76`:`82`,
        names_to = "generation",
        values_to = "pop_size"
      ) %>%
      mutate(
        generation = as.integer(generation),
        log_pop = log(pop_size)
      ) %>%
      filter(is.finite(log_pop))

    pop_long <- pop_long %>%
      group_by(Group, generation) %>%
      #summarise(mean_pop = mean(pop_size, na.rm = TRUE), .groups = "drop")
      summarise(
        mean_pop = mean(pop_size, na.rm = TRUE),
        lower = quantile(pop_size, 0.25, na.rm = TRUE),
        upper = quantile(pop_size, 0.75, na.rm = TRUE),
        .groups = "drop"
      )
  
  })
  
  
  ####################
  ####### text output
  ####################
  
  
  output$sd_limits <- renderText({
    sim_data <- selected_data()
    
    if (nrow(sim_data) == 0) return("No data to calculate limits.")
    
    #mean_hapdiv <- mean(filtered$HapDiv, na.rm = TRUE)
    sd_hapdiv <- sqrt(0.001348966)
    lower <- round(0.6866629 - (input$sd * sd_hapdiv), 3)
    upper <- round(0.6866629 + (input$sd * sd_hapdiv), 3)
    
    haplower <- min(input$hap)
    hapupper <- max(input$hap)
    
    paste0("In this document, 'successful' outcomes based on 'Haplotype Diversity' will have both haplotype diversity, measured in generation 82, between values of: Lower = ", lower, ", Upper = ", upper, ", and number of haplotypes that are >= ", haplower, " and <=", hapupper )
  })
  

  output$popTable <- renderTable({
    pop_data
  })


  ####################
  ####### plots
  ####################
  
  output$pop_plot <- renderPlot({
    sim_data <- pop_summary()
    
    raw_mode <- input$plot_type == "raw"
    
    # Define the trending hard limits (as line segments)
    hard_limits_df <- data.frame(
      generation = c(76, 77, 81, 82),
      lower = c(100, 36, 36, 500),
      upper = c(5000, 750, 750, 2500)
    )
    
    # Log-transform if needed
    if (!raw_mode) {
      hard_limits_df$lower <- log(hard_limits_df$lower)
      hard_limits_df$upper <- log(hard_limits_df$upper)
      
      sim_data <- sim_data %>%
        mutate(
          mean_pop = log(mean_pop),
          lower = log(lower),
          upper = log(upper)
        )
    }
    
    p <- ggplot(sim_data, aes(x = generation, color = Group, fill = Group)) +
      # Confidence ribbon
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color = NA) +
      
      
      
      # Add hard limit dashed lines (upper and lower trends)
      geom_line(data = hard_limits_df, aes(x = generation, y = lower),
                inherit.aes = FALSE, linetype = "dashed", color = "grey", size = 0.4) +
      geom_line(data = hard_limits_df, aes(x = generation, y = upper),
                inherit.aes = FALSE, linetype = "dashed", color = "grey", size = 0.4) +
      
      # Solid thin lines for CI bounds
      geom_line(aes(y = lower), linetype = "solid", size = 0.2, show.legend = FALSE) +
      geom_line(aes(y = upper), linetype = "solid", size = 0.2, show.legend = FALSE) +
      
      # Mean population line (thicker)
      geom_line(aes(y = mean_pop), size = 1) +
      
      theme_minimal() +
      labs(
        title = if (raw_mode) "Average Population Sizes by Generation" else "Average log(Population) Sizes by Generation",
        x = "Generation",
        y = if (raw_mode) "Mean Population Size" else "log(Mean) Population Size",
      )
    
    p 
    
  })
  
  

  output$plot_caption <- renderText({
    "Shaded areas are 25th–75th percentile ranges; grey dashed lines show hard upper and lower limits."
  })



  output$pop_density_plot <- renderPlot({
    sim_data <- selected_data()

  
    ###
    #sim_data$Group <- factor(sim_data$Group, levels = c("success", "non-success"))

    pop_long <- sim_data %>%
      select(sim_id, Group, `76`:`82`) %>%
      pivot_longer(cols = `76`:`82`, names_to = "generation", values_to = "pop_size") %>%
      mutate(
        generation = as.integer(generation),
        log_pop = log(pop_size)
      ) %>%
      filter(is.finite(log_pop))

    # Plot base
    if (input$plot_type == "raw") {
      p2 <- ggplot(pop_long, aes(x = pop_size, fill = Group)) +
        geom_density(alpha = 0.4) +
        facet_wrap(~ generation, scales = "free") +
        # scale_x_continuous(
        #   limits = c(0, 1000),
        #   breaks = c(0, 250, 500, 750, 1000)
        # ) +
        labs(
          title = "Population Distributions per Generation by Group",
          x = "Population Size",
          fill = "Group"
        )
    } else {
      p2 <- ggplot(pop_long, aes(x = log_pop, fill = Group)) +
        geom_density(alpha = 0.4) +
        facet_wrap(~ generation, scales = "free") +
        # scale_x_continuous(
        #   limits = c(log(1), log(1000)),
        #   breaks = log(c(1, 10, 100, 1000)),
        #   labels = c("1", "10", "100", "1000")
        # ) +
        labs(
          title = "log(Population) Distributions per Generation by Group",
          x = "log(Population Size)",
          fill = "Group"
        )
    }

    p2 +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  })
  
  
  output$fg_loss_plot <- renderPlot({
    
    # if (input$checkbox) {
    #   sim_data <- selected_data() %>%
    #     mutate(
    #       run = row_number(),
    #       diff = smallestFG1 - (509/3),
    #       bar_color = case_when(
    #         Group == "success" ~ "darkblue",
    #         Group == "non-success" ~ "darkred",
    #         TRUE ~ "lightgrey"
    #       )
    #     )
    #   
    #   ggplot(sim_data, aes(x = run, y = diff, fill = bar_color)) +
    #     geom_col() +
    #     geom_hline(yintercept = 0, colour = "black") +
    #     labs(x = "Run", y = "Number of Individuals Lost in Smallest FG") +
    #     scale_fill_identity() + 
    #     scale_x_continuous(
    #       breaks = c(0, 10000, 20000, 30000)
    #     ) +
    #     facet_wrap(~ Group, ncol = 1) +  # One column, stacked vertically
    #     theme_minimal() 
    # }
    
      sim_data <- selected_data() %>%
        mutate(
          run = row_number(),
          diff = smallestFG1 - (509/3),
          bar_color = case_when(
            Group == "success" ~  "#FBADA8",
            Group == "non-success" ~ "#67D9DD",
            TRUE ~ "lightgrey"
          )
        )
      
      ggplot(sim_data, aes(x = diff, fill = bar_color)) +
        geom_histogram(binwidth = 5, color = "black", alpha = 0.8) +  # adjust binwidth
        geom_vline(xintercept = 0, colour = "black") +
        labs(x = "Number of Individuals Lost in Smallest FG", 
             y = "Count of Runs") +
        scale_fill_identity() +
        facet_wrap(~ Group, ncol = 1, scales = "free_y") +   # still separate by Group
        theme_minimal()
    
  })
  

  plot_div_var <- reactive({
    sim_data <- selected_data()

    div_cols <- paste0("div", 75:81)
    var_cols <- paste0("var", 75:81)

    # map gens -> calendar years (match your table)
    gen_years <- c(`div75`=1830, `div76`=1855, `div77`=1880, `div78`=1905,
                   `div79`=1930, `div80`=1955, `div81`=1980)
    var_years <- c(`var75`=1830, `var76`=1855, `var77`=1880, `var78`=1905,
                   `var79`=1930, `var80`=1955, `var81`=1980)

    # summarise mean div and var by Group × generation
    div_long <- sim_data %>%
      filter(Group %in% c("success","non-success"), !is.na(HapDiv)) %>%
      group_by(Group) %>%
      summarise(across(all_of(div_cols), ~ mean(.x, na.rm = TRUE)), .groups="drop") %>%
      pivot_longer(-Group, names_to="gen", values_to="div") %>%
      mutate(Year = gen_years[gen])

    var_long <- sim_data %>%
      filter(Group %in% c("success","non-success"), !is.na(HapVar)) %>%
      group_by(Group) %>%
      summarise(across(all_of(var_cols), ~ mean(.x, na.rm = TRUE)), .groups="drop") %>%
      pivot_longer(-Group, names_to="gen", values_to="var") %>%
      mutate(Year = var_years[gen])

    left_join(div_long, var_long, by=c("Group","Year")) %>%
      mutate(
        ymin = div - sqrt(pmax(var, 0)),   # treat var as σ² → ribbon = ±σ
        ymax = div + sqrt(pmax(var, 0))
      )
  })

  
  
  output$div_var_plot <- renderPlot({
    df <- plot_div_var()
    
    ggplot(df, aes(x = Year, y = div, group = Group, color = Group, fill = Group)) +
      # variance ribbons (lighter)
      geom_ribbon(aes(ymin = ymin, ymax = ymax),
                  alpha = 0.25, colour = NA) +
      # diversity lines
      geom_line(linewidth = 1.2) +
      geom_point(size = 2) +
      
      scale_color_manual(values = c(
        "success"     = "#FBADA8",
        "non-success" = "#67D9DD"
      )) +
      scale_fill_manual(values = c(
        "success"     = "#FBADA8",
        "non-success" = "#67D9DD"
      )) +
      
      labs(x = "Year", y = "Haplotype diversity",
           title = "Diversity with variance envelopes") +
      theme_minimal(base_size = 12) +
      theme(legend.title = element_blank())
  })
  




  output$happlot <- renderPlot({
    sim_data <- selected_data()


    # Reshape haplotype columns
    hap_long <- sim_data %>%
      select(sim_id, Group, hap75, hap76, hap77, hap78, hap79, hap80, hap81, hapCur) %>%
      rename(hap82 = hapCur) %>%
      pivot_longer(
        cols = hap75:hap82,
        names_to = "generation",
        values_to = "num_haplotypes"
      ) %>%
      mutate(generation = as.integer(gsub("hap", "", generation)))

    # Calculate mean number of haplotypes per generation per group
    hap_summary <- hap_long %>%
      group_by(Group, generation) %>%
      summarise(mean_hap = mean(num_haplotypes, na.rm = TRUE), .groups = "drop")

    # Plot
    ggplot(hap_summary, aes(x = generation, y = mean_hap, color = Group)) +
      geom_line(size = 1.2) +
      labs(
        title = "Average Number of Haplotypes by Generation",
        x = "Generation",
        y = "Mean Number of Haplotypes",
        color = "Group"
      ) +
      theme_minimal()
  })



  output$hapdiv_violin <- renderPlot({

      sim_data <- selected_data()

      # Make the plot
      ggplot(sim_data, aes(x = Group, y = HapDiv, fill = Group)) +
        geom_violin(trim = FALSE, alpha = 0.6) +
        geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
        #facet_wrap(~ panel, ncol = 2) +
        geom_hline(yintercept = 0.69, linetype = "solid", color = "grey", size = 1) +
        annotate("text", x = 0.5, y = 0.72, label = "Empirical diversity = 0.69", color = "grey", size = 4, hjust = 0) +
        labs(
          title = "Haplotype Diversity (HapDiv) Across Prevalence and Success Types",
          x = "",
          y = "Haplotype Diversity"
        ) +
        #scale_fill_manual(values = c("success" = "#4CAF50", "non-success" = "#F44336")) +
        theme_minimal() +
        theme(
          legend.position = "none",
          strip.text = element_text(size = 12, face = "bold")
        )


    })


  output$violin <- renderPlot({
    #req(input$fg_sort == "unsorted")  # Only show if unsorted

    sim_data <- selected_data()

    # Pivot FG1-FG3 at generation 82 to long format
    fg_long <- sim_data %>%
      select(sim_id, Group, FG1, FG2, FG3) %>%
      pivot_longer(cols = FG1:FG3, names_to = "FG", values_to = "size")

    # Make the violin plot
    ggplot(fg_long, aes(x = Group, y = size, fill = Group)) +
      geom_violin(trim = FALSE, scale = "width", alpha = 0.6) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
      facet_wrap(~FG, scales = "free_y") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +

      labs(
        title = "FG Sizes at Generation 82 by Success Group",
        x = "",
        y = "Size",
        fill = "Group"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
  })





  ####################
  ####### summary tables
  ####################
  
  
  
  ####################
  ####### overall 


  output$summary_table <- renderTable({
    sim_data <- selected_data()  # Get the currently selected dataset (based on prevalence and sorting inputs)
    
    
    
      # Grouped summary (success and non-success)
      sandf <- sim_data %>%
        filter(Group %in% c("success", "non-success"), !is.na(HapDiv)) %>%
        group_by(Group) %>%
        summarise(
          n = n(),
          hapCur = mean(hapCur, na.rm = TRUE),
          HapDiv = mean(HapDiv, na.rm = TRUE),
          ratio = mean(ratio, na.rm = TRUE),
          #`n. FG1` = mean(FG1, na.rm = TRUE),
          #`n. FG2` = mean(FG2, na.rm = TRUE),
          #`n. FG3` = mean(FG3, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(pertotal = n / nrow(sim_data)) %>%
        rename(
          `avg. n.hap` = hapCur,
          `avg. diversity` = HapDiv,
          `avg. FG freq.` = ratio,
          proportion = pertotal
        )
      
      # Total summary
      tot <- sim_data %>%
        filter(Group %in% c("success", "non-success"), !is.na(HapDiv)) %>%
        summarise(
          n = n(),
          hapCur = mean(hapCur, na.rm = TRUE),
          HapDiv = mean(HapDiv, na.rm = TRUE),
          ratio = mean(ratio, na.rm = TRUE)
      
        ) %>%
        mutate(
          proportion = n / nrow(sim_data),
          Group = "Total"
        ) %>%
        rename(
          `avg. n.hap` = hapCur,
          `avg. diversity` = HapDiv,
          `avg. FG freq.` = ratio
        ) %>%
        select(Group, `avg. n.hap`, `avg. diversity`, `avg. FG freq.`, n, proportion)
      
      # Combine
      final_summary <- bind_rows(sandf, tot) %>%
        mutate(Group = factor(Group, levels = c("success", "non-success", "Total"))) %>%
        arrange(Group)
      
      final_summary
  })
  

  
  ####################
  ####### FG loss by success type
  
  output$fg_table <- renderTable({
    
    sim_data <- selected_data()
    
    ## feeding ground data
    group_summary <- sim_data %>%
      filter(Group %in% c("success", "non-success"), !is.na(HapDiv)) %>%
      filter(ratio <= 0.01) %>%
      group_by(Group) %>%
      summarise(
        n = n(),
        `avg. n.hap` = mean(hapCur, na.rm = TRUE),
        `avg. diversity` = mean(HapDiv, na.rm = TRUE),
        `n. FG1` = sum(ratioFG1 <= 0.01, na.rm = TRUE),
        `n. FG2` = sum(ratioFG2 <= 0.01, na.rm = TRUE),
        `n. FG3` = sum(ratioFG3 <= 0.01, na.rm = TRUE),
        .groups = "drop"
      )
    
    ## non-success subgroups
    non_success_summary <- sim_data %>%
      filter(Group == "non-success", !is.na(HapDiv)) %>%
      mutate(
        Group = case_when(
          HapDiv < hapdiv_range()["lower"] ~ "non-success_low",
          HapDiv > hapdiv_range()["upper"] ~ "non-success_high",
          TRUE ~ NA_character_  # ignore mid-range
        )
      ) %>%
      filter(!is.na(Group), ratio <= 0.01) %>%
      group_by(Group) %>%
      summarise(
        n = n(),
        `avg. n.hap` = mean(hapCur, na.rm = TRUE),
        `avg. diversity` = mean(HapDiv, na.rm = TRUE),
        `n. FG1` = sum(ratioFG1 <= 0.01, na.rm = TRUE),
        `n. FG2` = sum(ratioFG2 <= 0.01, na.rm = TRUE),
        `n. FG3` = sum(ratioFG3 <= 0.01, na.rm = TRUE),
        .groups = "drop"
      )
    
    ## Total summary
    total_summary <- sim_data %>%
      filter(Group %in% c("success", "non-success"), !is.na(HapDiv), ratio <= 0.01) %>%
      summarise(
        n = n(),
        `avg. n.hap` = mean(hapCur, na.rm = TRUE),
        `avg. diversity` = mean(HapDiv, na.rm = TRUE),
        `n. FG1` = sum(ratioFG1 <= 0.01, na.rm = TRUE),
        `n. FG2` = sum(ratioFG2 <= 0.01, na.rm = TRUE),
        `n. FG3` = sum(ratioFG3 <= 0.01, na.rm = TRUE)
      ) %>%
      mutate(Group = "Total") %>%
      select(Group, everything())
    
    ## Combine and arrange all
    final_summary <- bind_rows(
      group_summary,
      non_success_summary,
      total_summary
    ) %>%
      mutate(Group = factor(Group, levels = c(
        "success", "non-success", "non-success_low", "non-success_high", "Total"
      ))) %>%
      arrange(Group)
    
    final_summary
    
    
  })
    
  
  
  ####################
  ####### Population sizes hap diversity
  
  output$pop_table <- renderTable({
    
    sim_data <- selected_data() 
    
    ## feeding ground data
    
    group_summary <- sim_data %>%
      filter(Group %in% c("success", "non-success"), !is.na(HapDiv)) %>%
      group_by(Group) %>%
      summarise(
        across(`75`:`82`, ~ round(mean(.)), .names = "{.col}"),
        .groups = "drop"
      ) %>%
      mutate(Group = factor(Group, levels = c("success", "non-success"))) %>%
      arrange(Group)
    
    total_summary <- sim_data %>%
      filter(Group %in% c("success", "non-success"), !is.na(HapDiv)) %>%
      summarise(
        across(`75`:`82`, ~ round(mean(.)), .names = "{.col}")
      ) %>%
      mutate(Group = "Total") %>%
      select(Group, everything())
    
    final_summary <- bind_rows(group_summary, total_summary) %>%
      mutate(Group = factor(Group, levels = c("success", "non-success", "Total"))) %>%
      arrange(Group)
    
    final_summary <- final_summary %>%
      mutate(across(`75`:`82`, as.integer))
    
    names(final_summary) <- c("1830", "1855", "1880", "1905", "1930", "1955", "1980", "2005", "Current")
    
    
    final_summary
    
    
    
  })
  
  
  ####################
  ####### diversity by generation
  
  output$div_table <- renderTable({

    sim_data <- selected_data()

    div_cols <- paste0("div", 75:81)

    group_summary <- sim_data %>%
      filter(Group %in% c("success", "non-success"), !is.na(HapDiv)) %>%
      group_by(Group) %>%
      summarise(
        across(all_of(div_cols), ~ mean(.x, na.rm = TRUE)),
        HapDiv = mean(HapDiv, na.rm = TRUE),
        .groups = "drop"
      )

    total_summary <- sim_data %>%
      filter(Group %in% c("success", "non-success"), !is.na(HapDiv)) %>%
      summarise(
        across(all_of(div_cols), ~ mean(.x, na.rm = TRUE)),
        HapDiv = mean(HapDiv, na.rm = TRUE)
      ) %>%
      mutate(Group = "Total") %>%
      relocate(Group)

    final_summary <- bind_rows(group_summary, total_summary) %>%
      mutate(Group = factor(Group, levels = c("success", "non-success", "Total"))) %>%
      arrange(Group) %>%
      mutate(across(-Group, ~ round(.x, 3)))

    names(final_summary) <- c("1830", "1855", "1880", "1905", "1930", "1955", "1980", "2005", "Current")

    final_summary
  })
 
  
  
  ####################
  ####################
  #################### Random examples
  ####################
  ####################
  
 ####################
 ####### Data table 
  
  sampled_data <- reactive({
    
    sim_data <- selected_data() 
    
    set.seed(input$numeric)
    
    sampled_data <- sim_data %>%
      group_by(Group) %>%
      filter(Group %in% c("success", "non-success")) %>%
      sample_n(5) %>%
      ungroup()
  
    
  }) 


  output$sample_table <- renderTable({
    sampled_data <- sampled_data()
    
    sampled_data %>%
      select(Group,
             hap75, div75,
             hap76, div76,
             hap77, div77,
             hap78, div78,
             hap79, div79,
             hap80, div80,
             hap81, div81,
             HapDiv, HapDiv, 
             `75`, `76`,
             `77`, `78`,
             `79`, `80`,
             `81`, `82`,
             FG1, FG2, FG3,
             numf1, numf2, numf3, ratio) %>%
      mutate(
        across(
          .cols = where(is.numeric) & !matches("div", ignore.case = TRUE),
          .fns  = ~ format(round(.), nsmall = 0)
        )
      )
    
  })
  
  
  
  output$sample_haps <- renderPlot({
    sampled_data <- sampled_data() %>%
      mutate(run_id = row_number())  # Add a unique ID for each run
    
    hap_long <- sampled_data %>%
      select(run_id, Group, hap75, hap76, hap77, hap78, hap79, hap80, hap81, hapCur) %>%
      rename(hap82 = hapCur) %>%
      pivot_longer(
        cols = hap75:hap82,
        names_to = "generation",
        values_to = "num_haplotypes"
      )
    
    hap_long %>%
      ggplot(aes(x = generation, y = num_haplotypes, group = run_id, colour = as.factor(run_id))) +  # group lines by run
      geom_line(alpha = 0.7) +
      facet_wrap(~ Group, ncol = 1) +
      theme_minimal() +
      theme(panel.spacing = unit(.25, "lines")) +
      labs(
        title = "Number of Haplotypes Trajectories",
        x = "Generation",
        y = "Number of Haplotypes"
      ) +
      theme(legend.position = "none")
  })
  
  
  output$sample_FGs <- renderPlot({
    
    if (input$fg_prevalence == "without") {
      data <- FG_losswithoutPrev
    } else {
      data <- FG_losswithPrev
    }
    
    sampled_data <- sampled_data()
    
    group_map <- sampled_data %>%
      select(Group, simSEED)
    
    # filter, join, reshape
    df_subset <- data %>%
      filter(simSEED %in% sampled_data$simSEED) %>%
      left_join(group_map, by = "simSEED") %>%
      pivot_longer(cols = starts_with("FG"),
                   names_to = "FG",
                   values_to = "value")
    
    # plot
    ggplot(df_subset, aes(x = generation, y = value, color = factor(simSEED))) +
      geom_line(alpha = 0.7) +
      facet_grid(Group ~ FG, scales = "free_y") +  # rows = Group, cols = FG
      theme_minimal() +
      theme(legend.position = "none") +  # hide legend if too busy
      labs(title = "FG trajectories by Group",
           x = "Generation",
           y = "number in sample")
    
  })
  
   
   


}

shinyApp(ui, server)
