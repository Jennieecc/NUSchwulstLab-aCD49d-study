# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----   Single cell transcriptomics reveals cerebrospinal fluid immune   -----
# -----  dysregulation during healthy brain aging and cognitive impairment -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 04-20-2022
# Written by: Emma Tapp, Natalie Piehl
# Summary: Define helper functions used across the project
#
#-------------------------------------------------------------------------------
# Define functions

# Load in libraries
suppressMessages({
  library("ggrepel")
  library("tidyverse")
  library("ggpubr")
  library("ggthemes")
  library("grid")
})

# Negation of %in%
'%!in%' <- Negate('%in%')

# Define publication theme plotting function
theme_Publication_blank <- function(base_size=12, base_family="") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(fill = "transparent",colour = NA),
           plot.background = element_rect(fill = "transparent",colour = NA),
           panel.border = element_rect(colour = NA, fill = "transparent"),
           axis.title = element_text(size = rel(1)),
           axis.title.y = element_text(angle=90,margin=margin(0,10,0,0)),
           axis.title.x = element_text(margin=margin(10,0,0,0)),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(size = 0.3),
           axis.line.x = element_line(size = 0.3, linetype = "solid", colour = "black"),
           axis.line.y = element_line(size = 0.3, linetype = "solid", colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA, fill="transparent"),
           legend.position = "bottom",
           legend.margin = margin(t = 10, unit='pt'),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#d8d8d8",fill="#d8d8d8")
   ))
} 

# Define plot exporting helper function
set_panel_size <- function(p=NULL, g=ggplotGrob(p), file=NULL, 
                           margin = unit(1,"mm"),
                           width=unit(7, "inch"), 
                           height=unit(5, "inch")){
  panels <- grep("panel", g$layout$name)
  panel_index_w<- unique(g$layout$l[panels])
  panel_index_h<- unique(g$layout$t[panels])
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  
  g$widths[panel_index_w] <-  rep(width,  nw)
  g$heights[panel_index_h] <- rep(height, nh)
  
  if(!is.null(file)) {
    ggplot2::ggsave(file, g,
                    width = convertWidth(sum(g$widths) + margin,
                                         unitTo = "in", valueOnly = TRUE),
                    height = convertHeight(sum(g$heights) + margin,
                                           unitTo = "in", valueOnly = TRUE), useDingbats=F,
                    dpi=300)
    invisible(g)
  }
}

# Define volcano plotting function
volcano_plot <- function(data, file = NULL, title = NULL,
                         padj.lim = NULL, lfc.lim = NULL, lfc.asymmetric = NULL,
                         padj.thresh = 0.01, lfc.thresh = 0.25, x_title = "avg_log2FC",
                         width=unit(4, "inch"), height=unit(4, "inch")) {
  # Create gene name column
  data$gene <- rownames(data)
  
  # Generate PFC scores
  data$PFC <- -log10(data$BH) * abs(data$avg_log2FC)
  PFC <- unique(data$PFC)
  data$PFC[data$PFC == Inf] <- sort(PFC, partial=length(PFC)-1)[length(PFC)-1]
  
  #Generate log padj column
  data$log_padj <- -log10(data$BH)
  
  # Define limits if not provided
  if (is.null(padj.lim)) {
    log.padj.lim <- unique(data$log_padj)[order(-unique(data$log_padj))][2]
  } else {
    log.padj.lim <- -log10(padj.lim)
  }
  if (is.null(lfc.lim)) {
    lfc.lim <- abs(data[order(-abs(data$avg_log2FC)),"avg_log2FC"][1])
  }
  
  # Generate color column
  data$color <- rep("black", nrow(data))
  data[which(data$BH <= padj.thresh &
               data$avg_log2FC > lfc.thresh), 'color'] <- "red"
    data[which(data$BH <= padj.thresh &
                 data$avg_log2FC < -lfc.thresh), 'color'] <- "blue"
      
    # Scale down genes outside of bounds
    data[which(data$log_padj > log.padj.lim), 'log_padj'] <- log.padj.lim
    data[which(data$avg_log2FC > lfc.lim), 'avg_log2FC'] <- lfc.lim
    data[which(data$avg_log2FC < -lfc.lim), 'avg_log2FC'] <- -lfc.lim
    
    # Generate asymmetric lfc limits if necessary
    if (is.null(lfc.asymmetric)) {
      lfc_lims <- c(-lfc.lim, lfc.lim)
    } else {
      lfc_lims <- lfc.asymmetric
    }
    
    # Plot data
    p <-
      ggplot(data,
             aes(
               x = avg_log2FC,
               y = log_padj,
               color = color,
               label = gene,
               size = PFC
             )) +
      theme_Publication_blank() +
      geom_hline(
        yintercept = -log10(padj.thresh),
        linetype = 2,
        color = "gray"
      ) +
      geom_vline(xintercept = lfc.thresh,
                 linetype = 2,
                 color = "gray") +
      geom_vline(
        xintercept = -lfc.thresh,
        linetype = 2,
        color = "gray"
      ) +
      geom_point(aes(size = PFC), alpha = 0.5) +
      scale_color_manual(values = c("red" = "red",
                                          "black" = "black",
                                          "blue" = "blue")) +
                                            geom_text_repel(
                                              data = data[which(data$color != 'black'), ],
                                              inherit.aes = T,
                                              color = 'black',
                                              size = 6,
                                              force = 3
                                            ) +
      theme(legend.position = "none") +
      labs(title = title,
           x = x_title) +
      scale_x_continuous(limits = lfc_lims) +
      scale_y_continuous(limits = c(0, log.padj.lim)) 
    theme(axis.text.x = element_text(size = 21)) +
      theme(axis.text.y = element_text(size = 21))
    
    # Export plot
    if (is.null(file)) {
      return(p)
    } else {
      set_panel_size(
        p,
        file = file,
        width = width,
        height = height
      )
      return(NULL)
    }
}





suppressPackageStartupMessages({
  library(rlang)
})

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))  
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}