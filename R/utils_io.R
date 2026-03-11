# ==========================================
# Utility functions for file input/output
# ==========================================

# Read contour anchor coordinates
read_contour_anchors <- function(csv_file){

  if(!file.exists(csv_file)){
    stop("Contour anchor file not found")
  }

  df <- read.csv(csv_file, stringsAsFactors = FALSE)

  required_cols <- c("x","y")

  if(!all(required_cols %in% colnames(df))){
    stop("Contour anchor CSV must contain columns: x,y")
  }

  return(df[,required_cols])
}



# Save spatial anchors selected interactively
save_spatial_anchors <- function(anchor_df, output_file){

  write.csv(anchor_df, output_file, row.names = FALSE)

}



# Load previously saved spatial anchors
load_spatial_anchors <- function(input_file){

  if(!file.exists(input_file)){
    stop("Spatial anchor file not found")
  }

  df <- read.csv(input_file, stringsAsFactors = FALSE)

  required_cols <- c("x","y")

  if(!all(required_cols %in% colnames(df))){
    stop("Spatial anchor CSV must contain columns: x,y")
  }

  return(df[,required_cols])
}



# Save R session info for reproducibility
write_session_info <- function(output_file){

  info <- capture.output(sessionInfo())

  writeLines(info, output_file)

}
