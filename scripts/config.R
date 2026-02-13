######################################################
############    Libraries      #######################
######################################################
library(ggplot2)
library(Seurat)
library(SeuratObject)
library(rlang)
library(forcats)
library(patchwork)
library(gridExtra)
library(ggcorrplot)
library(clustree)
library(gtools)
library(ArchR)
library(RColorBrewer)
library(scales)
library(tools)
library(rlang)
library(tidyverse)
library(ggthemes)
library(future)
library(ggsankey)
library(pheatmap)
library(ontologyIndex)
library(rrvgo)
library(ggrepel)
library(ggh4x)
library(corto)
library(ggpubr)
library(muscat)
library(SingleCellExperiment)
library(edgeR)
library(xlsx)
library(muscatWrapper)
library(ggh4x)
library(igraph)
library(ggraph)
library(purrr)
library(ggVennDiagram)
library(scater)
library(harmony)
library(lisi)

######################################################
############  Data directories #######################
######################################################


# Directory paths
DATA_DIR <- file.path("..","..", "data","seurat_objects")

# File names
ATH_FILE <- "merge_at_nema.rds"
OSA_FILE <- "merge_os_nema.rds"
OSA_MOCK_FILE <- "merge_os_mock.rds"
SPECIE_FILE <- "specie_combined_sct.rds"


# Full paths
ATH_DATA_PATH <- file.path(DATA_DIR, ATH_FILE)
OSA_DATA_PATH <- file.path(DATA_DIR, OSA_FILE)
OSA_DATA_PATH_MOCK <- file.path(DATA_DIR, OSA_MOCK_FILE)
SPECIE_DATA_PATH <- file.path(DATA_DIR, SPECIE_FILE)

######################################################
######### Plots themes and colors     ################
######################################################

theme_paper <- theme_bw() + 
  theme(
    legend.position = "none",
    text = element_text(size=5, family="sans"),
    axis.text = element_text(size=5, family="sans"),
    element_line(linewidth = 0.1),
    panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.1),
    strip.background = element_rect(colour = "black", linewidth = 0.1)
  )

theme_notebook <- theme_bw() + 
  theme(
    legend.position = "none",
    text = element_text(size=10, family="sans"),
    axis.text = element_text(size=10, family="sans"),
    element_line(linewidth = 0.1),
    panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.1),
    strip.background = element_rect(colour = "black", linewidth = 0.1)
  )

theme_strip <- theme_bw() + 
  theme(
    legend.position = "none",
    text = element_text(size=5, family="sans"),
    axis.text = element_text(size=5, family="sans"),
    element_line(linewidth = 0.1),
    strip.text.x = element_text(margin = margin(b = 1, t = 2), size=5),
    panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.1),
    strip.background = element_rect(colour = "black", linewidth = 0.1)
  )

#Colors paper
cols_treat <- c("Infected"="#DE77AE","Mock"="grey90")
cols_rep <- c("R1"="#F28E2B","R2" = "59A14F","R3"="#E15759")
cols_rep_r0 <- c("R0"="#4E79A7","R1"="#F28E2B","R2" = "59A14F","R3"="#E15759")
cols_nema_at <- c("At28"="#DFC27D","At8"="#F6E8C3","Other"="grey90")
cols_nema_os<-c("Os17"="#C7EAE5","Os16"="#80CDC1","Os10"="#35978F","Os9"="#01665E","Os3"="#003C30","Other"="grey90")
cols_nema_combined <- c("At28"="#DFC27D","At8"="#F6E8C3","Other"="grey90", 
                        "Os17"="#C7EAE5","Os16"="#80CDC1","Os10"="#35978F","Os9"="#01665E","Os3"="#003C30")
cols_specie <- c("Arabidopsis thaliana"="#FFDB99","Oryza sativa"="#5a7667")
cols_sankey <- c("Other"="grey90",
                   "ath_Other"="grey90",
                   "osa_Other"="grey90",
                   "Os17"="#C7EAE5","Os16"="#80CDC1","Os10"="#35978F","Os9"="#01665E","Os3"="#003C30",
                   "At28"="#DFC27D","At8"="#F6E8C3",
                   "Int6"="#B2ABD2","Int7"="#8073AC")
cols_phase <- c("G1"="grey90","G2M"="#C51B7D", "S"= "#E6B400FF")
##correlation colors
paletteLength <- 50
corr_color <- colorRampPalette(c("grey10","grey30","grey50","grey70","grey90",
                              "#F7F7F7",
                              "#BFD3E6","#9EBCDA","#8C96C6","#8856A7","#810F7C"))(paletteLength)
cols_umap <- c('#324B63',
               '#daddee',
               '#b6cde3',
               '#92bbcd',
               '#6ba9b4',
               '#51887c',
               '#66a96d',
               '#a6c87b',
               '#d9e68f',
               '#ecaa6d',
               '#d56247',
               '#a43c3c',
               '#59160e',
               '#773c35',
               '#9e6257',
               '#d3847a',
               '#db96ad',
               '#c277a3',
               '#9b5d99',
               '#62417d',
               '#0c0f32',
               '#252a6f',
               '#34539e',
               '#457bb4',
               '#68a5d5',
               '#92cade',
               '#97a6b7',
               '#7b8aa4',
               '#5a6685')

cols_ct_at <- c("Epidermis"="#FABFD2",
                "Atrichoblast"="#D37295",
                "Trichoblast"="#B07AA1",
                "Columella"="#4E79a7",
                "Cortex"="#F1CE63",
                "Endodermis"="#9D7660",
                "Pericycle"="#E15759",
                "Procambium"="#59A14f",
                "Phloem"="#9467bd",
                "Xylem"="#ff7f0e",
                "LRC"="#A0CBE8",
                "Initials"="#499894",
                "SC"="black")

# cols_ct_at <- c("Epidermis"="#B5CF6B",
#                 "Atrichoblast"="#BD9E39",
#                 "Trichoblast"="#CEDB9C",
#                 "Columella"="#4E79a7",
#                 "Cortex"="#E7CB94",
#                 "Endodermis"="#843C39",
#                 "Pericycle"="#D6616B",
#                 "Procambium"="#E7969C",
#                 "Phloem"="#A55194",
#                 "Xylem"="#DE9ED6",
#                 "LRC"="#9C9EDE",
#                 "Initials"="#637939",
#                 "SC"="black")

cols_ct_os <- c("Epidermis"="#FABFD2",
                "Atrichoblast"="#D37295",
                "Trichoblast"="#B07AA1",
                "Exodermis"= "#C8AA52",
                "Sclerenchyma"="#847036",
                "RC"="#4E79a7",
                "Cortex"="#F1CE63",
                "Endodermis"="#9D7660",
                "Pericycle"="#E15759",
                "Phloem"="#9467bd",
                "Xylem"="#ff7f0e",
                "Initials"="#499894")

# cols_ct_os <- c("Epidermis"="#B5CF6B",
#                 "Atrichoblast"="#BD9E39",
#                 "Trichoblast"="#CEDB9C",
#                 "RC"="#4E79a7",
#                 "Cortex"="#E7CB94",
#                 "Endodermis"="#843C39",
#                 "Pericycle"="#D6616B",
#                 "Exodermis"="#E7969C",
#                 "Phloem"="#A55194",
#                 "Xylem"="#DE9ED6",
#                 "Sclerenchyma"="#9C9EDE",
#                 "Initials"="#637939")


cols_treatment_regulons <- c("Mock"="#7FBC41",
                             "Shared"="#CFD2DA",
                             "Infected"="#DE77AE")
palette_pink <- rev(brewer.pal(11,"PiYG")[1:6])
palette_cyan <- c("#F7F7F7","#6DCFCFFF","#2D9C9CFF","#0A6969FF","#00393AFF")
palette_yellow <-  c("#F7F7F7","#F2E8C4FF","#EED682FF","#EAC541FF","#E6B400FF")

######################################################
###Functions for subsampling and clustering  #########
######################################################

## Function to find the cluster resolution for a determined number of communities in a Seurat object

FindClusterResolution <- function(object, nclusters = 20, initial_res = 0.1, res_increment = 0.1){
  # Check if the inputs are valid
  if (!is(object, "Seurat")) stop("object must be a Seurat object")
  if (!is.numeric(nclusters) || nclusters <= 0) stop("nclusters must be a positive number")
  if (!is.numeric(initial_res) || initial_res <= 0) stop("initial_res must be a positive number")
  if (!is.numeric(res_increment) || res_increment <= 0) stop("res_increment must be a positive number")
  
  res <- initial_res
  nclusters_y <- 0
  
  while(nclusters_y < nclusters){
    new_res_obj <- FindClusters(object, resolution = res)
    nclusters_y <- nlevels(new_res_obj@active.ident)
    res <- res + res_increment
  }
  
  gc() # Call garbage collector once at the end
  return(new_res_obj)
}

DefaultSeuratPipeline <- function(count_matrix, dims_umap=1:30, metadata=NULL,...){
  # Check if the inputs are valid
  if (!is.matrix(count_matrix) && !is(count_matrix, 'dgCMatrix')) stop("count_matrix must be a matrix or a sparseMatrix")
  if (!is.null(metadata) && !is.data.frame(metadata)) stop("metadata must be a data frame or NULL")
  
  object_seurat <- CreateSeuratObject(count_matrix) %>% 
    SCTransform() %>% 
    RunPCA() %>% 
    RunUMAP(dims=dims_umap) %>% 
    FindNeighbors() %>% 
    FindClusterResolution(...)
  
  if(!is.null(metadata)){
    colnames(metadata) <- paste("reference", colnames(metadata),sep="_")
    object_seurat <- AddMetaData(object_seurat, metadata)
    object_metadata <- object_seurat@meta.data %>% select(c(starts_with("SCT_snn"),colnames(metadata)))
  } else{
    object_metadata <- object_seurat@meta.data %>% select(starts_with("SCT_snn"))
  }
  
  return(object_metadata)
}

SubsampleSeurat <- function(sc_obj, variable, clusters, size_subset = 10000, boot_rep = 10,nclusters=NULL,...) {
  original_metadata <- sc_obj[[]][,c(variable,clusters)]
  original_matrix <- sc_obj@assays$RNA@counts
  
  if(is.null(nclusters)){
    nclusters <- nlevels(original_metadata[,clusters])
  } else {
    nclusters <- nclusters
  }
  
  
  # Divide matrix by variable of interest
  grouped_cells <- split(rownames(original_metadata), original_metadata[,variable])
  
  list_matrix <- lapply(grouped_cells, function(x){
    original_matrix[,colnames(original_matrix) %in% x]})
  
  # Create subsamplings of the cells
  list_subset_matrix <- lapply(1:boot_rep, function(i) {
    do.call(cbind, lapply(list_matrix, function(x) {x[,sample(ncol(x), size_subset, replace = F)]}))
  })
  
  # Generate new Seurat objects
  list_subset_objects <- lapply(list_subset_matrix, 
                                DefaultSeuratPipeline, 
                                dims_umap = 1:30, 
                                metadata = original_metadata,
                                nclusters = nclusters)
  
  colnames<-(c("new_cluster",variable,"ref_cluster"))
  list_subset_objects<-lapply(list_subset_objects, setNames, colnames)
  return(list_subset_objects)
}





#' Create a Sankey-style UMAP visualization
#' @param seurat_obj A Seurat object
#' @param group.by The metadata column to group cells by
#' @param clusters Optional vector of clusters to subset
#' @param umap_cols Names of UMAP dimensions (default: c("UMAP_1", "UMAP_2"))
#' @param colors Named vector of colors for groups
#' @param point_size Size of points in UMAP (default: 1)
#' @param point_alpha Alpha of points in UMAP (default: 0.7)
#' @param line_size Size of connecting lines (default: 0.1)
#' @param line_alpha Alpha of connecting lines (default: 0.1)
#' @param sankey_position Position of Sankey bar (default: 8)
#' @return A ggplot object
sankey_umap <- function(seurat_obj,
                        group.by,
                        clusters = NULL,
                        umap_cols = c("UMAP_1", "UMAP_2"),
                        colors = NULL,
                        point_size = 1,
                        point_alpha = 0.7,
                        line_size = 0.1,
                        line_alpha = 0.1,
                        sankey_position = 8) {
  
  # Extract UMAP coordinates and metadata
  umap_data <- as.data.frame(Embeddings(seurat_obj, "umap"))
  colnames(umap_data) <- umap_cols
  metadata <- seurat_obj[[group.by]]
  
  # Combine data
  plot_data <- cbind(umap_data, group = metadata[,group.by])
  plot_data$Cluster <- Idents(seurat_obj)
  
  # Filter clusters if specified
  if (!is.null(clusters)) {
    plot_data <- plot_data[plot_data$Cluster %in% clusters,]
  }
  
  # Order groups by Y position
  order_groups <- plot_data %>%
    group_by(group) %>%
    summarize(mean_y = mean(!!sym(umap_cols[2]))) %>%
    arrange(desc(mean_y))
  
  order_cells <- plot_data %>%
    rownames_to_column("cell_name") %>%
    arrange(group = factor(group, levels = unique(order_groups$group))) %>%
    group_by(group) %>%
    arrange(desc(!!sym(umap_cols[2])), .by_group = TRUE)
  
  # Calculate Sankey positions
  y_range <- range(plot_data[[umap_cols[2]]])
  diff <- abs(diff(y_range)) * 2
  y_positions <- rev(seq(y_range[1] * 2, y_range[2] * 2, 
                         length.out = nrow(order_cells)))
  
  # Create Sankey data
  data_sankey <- data.frame(
    cell_name = order_cells$cell_name,
    group = order_cells$group,
    Cluster = order_cells$Cluster
  ) %>%
    mutate(
      !!sym(umap_cols[2]) := y_positions,
      !!sym(umap_cols[1]) := sankey_position
    )
  
  # Combine data for plotting
  all_data <- rbind(
    order_cells[, c("cell_name", "group", umap_cols, "Cluster")],
    data_sankey
  )
  
  # Use provided colors or generate default
  if (is.null(colors)) {
    n_groups <- length(unique(plot_data$group))
    colors <- setNames(
      hcl.colors(n_groups, "Temps"),
      unique(plot_data$group)
    )
  }
  
  # Create plot
  ggplot() +
    geom_point(data = plot_data,
               aes(x = !!sym(umap_cols[1]), y = !!sym(umap_cols[2])),
               size = point_size, shape = 21, stroke = 0,
               alpha = point_alpha, fill = "black") +
    geom_line(data = all_data,
              aes(x = !!sym(umap_cols[1]), y = !!sym(umap_cols[2]),
                  color = group, group = cell_name),
              linewidth = line_size, alpha = line_alpha) +
    geom_tile(data = data_sankey,
              aes(x = !!sym(umap_cols[1]) + 0.5,
                  y = !!sym(umap_cols[2]),
                  fill = group),
              linewidth = 0) +
    theme_void() +
    theme(legend.position = "none") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors)
}




sankey_umap_fraction <- function(seurat_obj,
                                 group.by,
                                 clusters = NULL,
                                 umap_cols = c("UMAP_1", "UMAP_2"),
                                 colors = NULL,
                                 point_size = 1,
                                 point_alpha = 0.7,
                                 line_size = 0.1,
                                 line_alpha = 0.1,
                                 sankey_position = 8,
                                 label_size = 3,
                                 show_labels = TRUE,
                                 x_text_margin = 3) { 
  
  # Extract UMAP coordinates and metadata
  umap_data <- as.data.frame(Embeddings(seurat_obj, "umap"))
  colnames(umap_data) <- umap_cols
  metadata <- seurat_obj[[group.by]]
  
  # Combine data
  plot_data <- cbind(umap_data, group = metadata[,group.by])
  plot_data$Cluster <- Idents(seurat_obj)
  
  # Filter clusters if specified
  if (!is.null(clusters)) {
    plot_data <- plot_data[plot_data$Cluster %in% clusters,]
  }
  
  # Order groups by Y position
  order_groups <- plot_data %>%
    group_by(group) %>%
    summarize(mean_y = mean(!!sym(umap_cols[2])), .groups = "drop") %>%
    arrange(desc(mean_y))
  
  order_cells <- plot_data %>%
    rownames_to_column("cell_name") %>%
    arrange(group = factor(group, levels = unique(order_groups$group))) %>%
    group_by(group) %>%
    arrange(desc(!!sym(umap_cols[2])), .by_group = TRUE)
  
  # Calculate Sankey positions
  y_range <- range(plot_data[[umap_cols[2]]])
  diff <- abs(diff(y_range)) * 2
  y_positions <- rev(seq(y_range[1] * 2, y_range[2] * 2, 
                         length.out = nrow(order_cells)))
  
  # Create Sankey data
  data_sankey <- data.frame(
    cell_name = order_cells$cell_name,
    group = order_cells$group,
    Cluster = order_cells$Cluster
  ) %>%
    mutate(
      !!sym(umap_cols[2]) := y_positions,
      !!sym(umap_cols[1]) := sankey_position
    )
  
  # Combine data for plotting
  all_data <- rbind(
    order_cells[, c("cell_name", "group", umap_cols, "Cluster")],
    data_sankey
  )
  
  # Calculate percentages and mean positions for each group
  group_percentages <- data_sankey %>%
    group_by(group) %>%
    summarise(
      count = n(),
      percentage = n() / nrow(data_sankey) * 100,
      y_pos = mean(!!sym(umap_cols[2])),
      label = sprintf("%s (%.1f%%)", group, n() / nrow(data_sankey) * 100),
      .groups = "drop"
    )
  
  # Use provided colors or generate default
  if (is.null(colors)) {
    n_groups <- length(unique(plot_data$group))
    colors <- setNames(
      hcl.colors(n_groups, "Temps"),
      unique(plot_data$group)
    )
  }

  
  # Get the x-axis range from the data
  x_range <- range(plot_data[[umap_cols[1]]])
  x_limit_right <- sankey_position + x_text_margin # Extend right limit for labels
  # Calculate percentages ONCE per group (not per cell)
  group_stats <- plot_data %>%
    group_by(group) %>%
    summarise(
      count = n(),
      percentage = n() / nrow(plot_data) * 100,
      .groups = "drop"
    )
  
  # Calculate mean positions for labels
  group_positions <- data_sankey %>%
    group_by(group) %>%
    summarise(
      y_pos = mean(!!sym(umap_cols[2])),
      .groups = "drop"
    )
  
  # Combine stats and positions
  label_data <- group_stats %>%
    left_join(group_positions, by = "group") %>%
    mutate(label = sprintf("%s (%.1f%%)", group, percentage))
  
  # Create plot
  p <- ggplot() +
    geom_point(data = plot_data,
               aes(x = !!sym(umap_cols[1]), y = !!sym(umap_cols[2])),
               size = point_size, shape = 21, stroke = 0,
               alpha = point_alpha, fill = "black") +
    geom_line(data = all_data,
              aes(x = !!sym(umap_cols[1]), y = !!sym(umap_cols[2]),
                  color = group, group = cell_name),
              linewidth = line_size, alpha = line_alpha) +
    geom_tile(data = data_sankey,
              aes(x = !!sym(umap_cols[1]) + 0.5,
                  y = !!sym(umap_cols[2]),
                  fill = group),
              linewidth = 0) +
    coord_cartesian(xlim = c(min(plot_data[[umap_cols[1]]]), 
                             sankey_position + x_text_margin))
  
  # Add labels if requested
  if (show_labels) {
    p <- p + 
      geom_text(data = label_data,  # Use the new label_data
                aes(x = sankey_position + 1,
                    y = y_pos,
                    label = label),
                size = label_size,
                hjust = 0)
  }
  
  # Add theme and scales
  p + theme_void() +
    theme(legend.position = "none",
          plot.margin = margin(5.5, 40, 5.5, 5.5)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors)
}

#Jaccard index function

#Function to calculate jaccard indexes when comparing two lists

jaccard_pvalue <- function(list1, list2, num_samples = 1000) {
  # Initialize matrices
  jaccard_matrix <- array(NA, dim = c(length(list1), length(list2), num_samples),
                          dimnames = list(names(list1), names(list2), paste0("Rep_", 1:num_samples)))
  original_jaccard_indexes <- matrix(NA, nrow = length(list1), ncol = length(list2), 
                                     dimnames = list(names(list1), names(list2)))
  
  # Calculate original Jaccard indexes and random samples
  for (i in 1:length(list1)) {
    for (j in 1:length(list2)) {
      intersection_size <- length(intersect(list1[[i]], list2[[j]]))
      union_size <- length(union(list1[[i]], list2[[j]]))
      original_jaccard_indexes[i, j] <- intersection_size / union_size
      
      for (rep in 1:num_samples) {
        shuffled_list1 <- sample(unlist(list1), replace = FALSE)
        shuffled_list2 <- sample(unlist(list2), replace = FALSE)
        
        shuffled_list1 <- split(shuffled_list1, rep(1:length(list1), sapply(list1, length)))
        shuffled_list2 <- split(shuffled_list2, rep(1:length(list2), sapply(list2, length)))
        
        intersection_size <- length(intersect(shuffled_list1[[i]], shuffled_list2[[j]]))
        union_size <- length(union(shuffled_list1[[i]], shuffled_list2[[j]]))
        jaccard_matrix[i, j, rep] <- intersection_size / union_size
      }
    }
  }
  
  # Calculate p-values
  p_val_matrix <- matrix(NA, nrow = length(list1), ncol = length(list2), 
                         dimnames = list(names(list1), names(list2)))
  
  for (i in 1:length(list1)) {
    for (j in 1:length(list2)) {
      p_val_matrix[i, j] <- mean(original_jaccard_indexes[i, j] < jaccard_matrix[i, j, ])
    }
  }
  
  return(list(original_jaccard_indexes = original_jaccard_indexes, p_val_matrix = p_val_matrix))
}


######################################
##Hypergeometric test function ######
#######################################

hypergeometric_test <- function(genes, sets, background = NULL, adjust_method = "BH",lower.tail=FALSE) {
  # Define the full background if not provided
  if (is.null(background)) {
    background <- unique(unlist(sets))
  }

  # Filter relevant sets (those with overlap)
  relevant_sets <- sets[sapply(sets, function(x) any(x %in% genes))]

  # Initialize p-values and fold-change list
  results <- list()

  # Perform the hypergeometric test
  for (i in seq_along(relevant_sets)) {
    set_name <- names(relevant_sets)[i]
    current_set <- relevant_sets[[i]]
    overlap <- intersect(genes, current_set)
    observed_overlap <- length(overlap)

    # Hypergeometric test
    pval <- phyper(
      q = observed_overlap - 1,
      m = length(current_set),
      n = length(background) - length(current_set),
      k = length(genes),
      lower.tail = lower.tail
    )

    # Calculate expected overlap
    expected_overlap <- (length(current_set) * length(genes)) / length(background)

    # Calculate fold-change enrichment
    fold_change <- observed_overlap / expected_overlap

    # Store results
    results[[set_name]] <- c(
      pval = pval,
      fold_change = fold_change,
      observed = observed_overlap,
      expected = expected_overlap
    )
  }

  # Convert to data frame
  res <- do.call(rbind, results) %>%
    as.data.frame() %>%
    rownames_to_column("Set") %>%
    mutate(p_adj = p.adjust(pval, method = adjust_method)) # Adjust p-values

  return(res)
}
