facet_no_y_theme <- theme( # for the first hit-calling plot, the most austere
  text = element_text(size = 8),
  axis.title = element_blank(), # don't label the axes
  axis.text.x = element_text(), # don't label the numbers on the axes
  axis.ticks = element_blank(), # dont have ticks
  legend.position = "right", # put a legent on the right of the plot
  plot.title = element_text(lineheight=.8, face="bold", size = 12), # have a title
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),
  panel.spacing.y = unit(0.1, "lines"),
  aspect.ratio = (1/1.618)
)

### the main plot
# facet_wrap_linetyped2(values$df_all, title = values$plot_title, facets_wide = 6 )
facet_wrap_linetyped2 <- function(df_melt, title, facets_wide) {
  p <- ggplot(df_melt, aes(x = Temperature, # temperature on X
                           y = value, # RFU on y
                           color = channel_f, # colored by the state
                           linetype = type,
                           group = dye_conc_type_channel # group means series, as in, this defines the unique data sets
  )) +
    geom_line(size = 0.3, alpha = 0.8) + # change the line type depending on the dye concentration # linetype = df_melt$conc #
    facet_wrap(~dye, scales = "free", ncol = facets_wide) +
    labs(title = title, color = "Channel") +
    theme_bw() +
    scale_color_manual( values = c("Cy5.5" = "#67000d", "Cy5" = "#a50f15", "ROX" = "#ef3b2c", "TAMRA" = "#f16913", "JOE" = "#74c476", "FAM" = "#2171b5")) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    facet_no_y_theme
  
  p # return the plot
}

facet_grid_validation <- function(df_melt, title, facets_wide) {
  p <- ggplot(df_melt, aes(x = Temperature, # temperature on X
                           y = value, # RFU on y
                           color = conc, # colored by the state
                           linetype = protein,
                           group = dye_conc_type_channel # group means series, as in, this defines the unique data sets
  )) +
    geom_line(size = 0.3, alpha = 0.8) + # change the line type depending on the dye concentration # linetype = df_melt$conc #
    facet_grid(rows = dye, cols = channel_f, scales = "free") +
   # facet_wrap(~dye, scales = "free", ncol = facets_wide) +
    labs(title = title, color = "[Dye] (µM)") +
    scale_color_viridis_c(begin = 0.2, end = 1) +
    theme_bw() +
    #scale_color_manual( values = c("Cy5.5" = "#67000d", "Cy5" = "#a50f15", "ROX" = "#ef3b2c", "TAMRA" = "#f16913", "JOE" = "#74c476", "FAM" = "#2171b5")) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    facet_no_y_theme
  
  p # return the plot
}

save_screen_meta_data <- function (path_in,
                                   protein,
                                   exp_num,
                                   daughter_num,
                                   script_version,
                                   instrument,
                                   plate_type,
                                   plate_lot,
                                   additional_notes,
                                   save_outputs,
                                   buffer_type,
                                   facets_wide
) {
  
  # define and/or create directory logic
  path_int <- paste0(path_in, "/intermediate/") # holds intermediate data files
  path_fin <- paste0(path_in, "/final/") # holds all final materials
  path_fin_facets <- paste0(path_fin, "/hits/") # holds all final materials
  
  #mkdir(path_raw) # if exists, fails silently ("FALSE")
  mkdir(path_int)
  mkdir(path_fin)
  mkdir(path_fin_facets)
  
  # write the read_me files
  # generate a read me  and session Info fils e for this experiment with relevant specs immediately
  fileConn <- file(paste0(path_fin, exp_num, "_", as.character(base::Sys.Date()),"_",  "readme.txt"))
  writeLines(c(paste("Experiment", paste0(exp_num, protein), sep = ": "),
               paste("Dye daughter plates:", daughter_num, sep = ": "),
               paste("Buffer used for screen:", buffer_type, sep = ":"),
               paste("Hit calling script version", script_version, sep = ": "),
               paste("Instrument used", instrument, sep = ": "),
               paste("Plate type", plate_type, sep = ": "),
               paste("Plate lot", plate_lot, sep = ": "),
               paste("Analysis completed on", as.character(base::Sys.Date()), sep = ": "),
               additional_notes), fileConn)
  close(fileConn)
  
  writeLines(capture.output(sessionInfo()), paste0(path_fin, exp_num, as.character(base::Sys.Date()),"_","sessionInfo.txt")) # save the session info into the final materials folder
}

make_channel_vec <- function( df ) { # make the vector which specifies channel for each reading
  channels <- df %>%
    group_by(`0`) %>%
    filter(n() == 1) %>%
    select(`0`) %>%
    as_vector()
  
  n_meas <- df %>%  # the number of wells measured per channel (always the same for all channels )
    group_by(`0`) %>%
    filter(n() > 1) %>%
    tally() %>%
    nrow()
  
  rep(channels , each = n_meas + 1) # add one, for the row which will still contain the channel itself
}

read_qTower <- function( file_path ) {
  
  df_raw <- read_csv(file_path,
                     col_names = c(0:500) %>% as.character()) %>% # read in 500 columns; this should exceed any actual run, and fill in columsn to right as NAs
    select_if(~sum(!is.na(.)) > 0) #%>% # remove all columns which are all NAs
  
  
  df <- df_raw %>%
    drop_na( tail(names(.), 1))  %>% # drop the header, which is empty in the tailing cols
    ## drop_na( tail(names(.), 1) %>% var() ) %>%  # this line worked in an earlier version of dplyr i think?
    mutate( channel = make_channel_vec(.)) %>% # add channel as a column
    filter(!`0` %in% .$channel) %>%
    rename(well = `0`) %>%
    mutate_at(vars(-well, -channel), as.numeric) %>%
    pivot_longer(-c("well", "channel"), names_to = "Temperature", values_to = "value") %>%
    mutate_at(vars(well, channel), as.character) %>%
    mutate_at(vars(Temperature, value), as.numeric) %>%
    mutate(channel_f = factor(.$channel, levels = c("FAM", "JOE", "TAMRA", "ROX", "Cy5", "Cy5.5", "SyproOrange")))
}


# new daughter layout function
df_to_layout <- function(df, layout_type) {
  df_m <-   set_names( df ,  c("type","row",as.numeric( df [1,-c(1,2)]))) %>%
    . [ -1 , -1] %>%
    reshape2::melt( . ,id.vars = "row") %>%
    mutate( . , well = as_vector(map2( . $row,  . $variable, paste0)) ) %>%
    set_names( . , c("row", "column", layout_type, "well"))
  df_m
}

make_layout <- function( filename ) { # from path to raw layout to a final fomatted layout file
  # read the layout file, and split each layout into an individual
  layout_list <- data.table::fread( filename, header = TRUE) %>%
    as_tibble() %>%
    split( . ,  . $Type)
  
  # put into a merge-able form
  layout <- df_to_layout(layout_list[[1]], names(layout_list)[[1]])[c(1,2,4)] # initialize the list
  for (i in c(1:length(layout_list))) {
    layout <- layout %>%
      mutate("var" =  as_vector(df_to_layout(layout_list[[i]], layout_type = names(layout_list)[[i]])[3] )) %>% # append the column of interest
      set_names(c(names(layout), names(layout_list)[[i]])) # rename based on the column of interest
  }
  layout <- layout %>%
    unite("condition", c(4:ncol(.)), remove = FALSE) %>% # create a unique column, used to define groups after averaging
    mutate_if(is.factor, as.character)
  
  layout
}

make_well_names <- function(row_style, num_style) {
  if (row_style == "rows") { rows <-  letters[1:16] } else {rows <- LETTERS[1:16] }
  if (num_style == "01") { cols <- c(c("01", "02", "03", "04", "05", "06", "07", "08", "09"), c(10:24)) } else {cols <- c(1:24) }
  
  int <-  list(row = rows,
               col =cols) %>%
    cross() %>% # general all possible combos of rows and columns
    map(lift(paste0)) %>% # pate them together
    
    as_vector() %>%
    as_tibble()  %>%
    mutate(nchar = lapply(.$value, nchar) %>% as_vector(),
           char =  str_extract_all(.$value, "[A-Z; a-z]", simplify = TRUE) 
           %>% str_to_upper(locale = "en") 
           %>% as_vector()) 
  
  
  if (num_style == "01") {
    ordered <- int %>%
      dplyr::arrange(value) %>%
      .[[1]]
  } else {
    ordered <- int %>%
      dplyr::arrange(char) %>%
      .[[1]]
  }
  ordered
}

alpha_to_num <- function(alpha) {
  lets <- LETTERS[1:26]
  match(alpha, lets)
} 

set_dye_num <- function(row_nums) {
  floor(row_nums/4 - 0.1) + 1
}

convert_numerics <- function( vec ) {
  
  if(all(varhandle::check.numeric(vec))){
    # convert the vector to numeric
    vec <- as.numeric(vec)
  }
  vec
}

df_to_layout_maker <- function(df_layout, col_name, set_custom_name = FALSE, custom_name = "var_name") {
  name_of_col <- col_name # set the name of the column
  if (set_custom_name == TRUE) {
    name_of_col <- custom_name
  }
  
  df_dye_raw <- df_layout %>%
    mutate(Type = rep(name_of_col, times = nrow(.))) %>%
    select(Type, row, col, !! sym(col_name)) %>%
    pivot_wider(names_from = col, values_from = !! sym(col_name))
  
  df_dye_b <- as_tibble(t(df_dye_raw), rownames = "row_names")
  
  df_dye_c <-  as_tibble(t(df_dye_b )) %>%
    set_names(c("Type", names(.)[-1])) %>%
    mutate(Type = rep(.[4,1], times = nrow(.)) %>% as_vector())
}

# test_maker <- df_to_layout_maker(test_func, "conc", FALSE, "whatever")
# test_maker

df_to_layout_file <- function(df_layout_list, var_names = c(dye, conc)) {
  df_list <- list()
  layout_list <- map(var_names, df_to_layout_maker, df_layout_list)
  
}

make_layout_tibble <- function(df_valid, 
                               layout, 
                               final_vol = 500, 
                               protein_name, 
                               validate_all = TRUE,
                               fold_dil = 2,
                               high_conc_fold = 2) {
  # create the well table which is used to make the layouts
  well_table <- make_well_names("ROWS", "1") %>%
    tibble(well = .) %>%
    mutate(well_f = factor(well, levels = well),
           row = str_extract_all(well, "[A-Z; a-z]", simplify = TRUE) %>% as.character,
           col = parse_number(well)) %>%
    mutate(row_num = alpha_to_num(row)) %>%
    mutate(block_num = set_dye_num(row_num)) %>%
    mutate(conc_num = high_conc_fold/fold_dil^(row_num - block_num)) %>%
    group_by(col, block_num)
  
  if (validate_all == TRUE) {
    dyes_test_df <- df_valid %>%
      filter(assignment != "none")
  } else {
    dyes_test_df <- df_valid %>%
      filter(assignment == "validate")
  }
  
  dyes_test <- dyes_test_df %>%
    select(dye) %>%
    add_row(dye = "SYPRO", .before = 1) %>%
    distinct(dye) %>%
    mutate(dye_num = c(1:nrow(.))) # add a validation #
  
  df <-  well_table %>%
    nest() %>%
    ungroup() %>%
    mutate(dye_num = c(1:nrow(.))) %>%
    right_join(dyes_test, by = "dye_num") %>%
    unnest(cols = c(data)) %>%
    
    mutate(volume = final_vol,
           protein = protein_name) %>%
    left_join(layout %>% select(dye, conc), by = "dye")   %>%
    mutate_all(convert_numerics) %>%
    mutate(conc = conc*conc_num)%>%
    right_join(well_table) %>%
    replace_na(list(dye = "Empty", conc = "Empty", volume = "Empty", protein = "Empty")) %>%
    select(row, col, dye, conc, volume, protein) 
  
}