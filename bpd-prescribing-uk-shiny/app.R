#===============================================================================
# Create Shiny app that displays Sankey plots and heatmaps of treatment 
# patterns, including an interactive Sankey plot
# ==============================================================================

# ======== Packages ========
suppressPackageStartupMessages({
  library(shiny)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
  library(scales)
  library(stringr)
  library(networkD3)
  library(htmlwidgets)
  library(tibble)
  library(jsonlite)
})

# ======== Load precomputed data ========
read_rds_or_stop <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing required file: %s", path))
  readRDS(path)
}

coerce_to_df <- function(x) {
  if (inherits(x, "data.frame")) return(x)
  if (is.list(x)) {
    if (length(x) == 1 && is.list(x[[1]]) && !inherits(x[[1]], "data.frame"))
      return(coerce_to_df(x[[1]]))
    for (el in x) if (inherits(el, "data.frame")) return(el)
    stop("Loaded object is a list with no data.frame elements.")
  }
  stop(sprintf("Unsupported object class: %s", paste(class(x), collapse = ", ")))
}

coerce_atomic_chr <- function(x) {
  if (is.list(x)) {
    vapply(
      x,
      function(v) {
        if (is.null(v) || length(v) == 0) return(NA_character_)
        if (is.list(v)) {
          if (length(v) == 0) return(NA_character_)
          v <- v[[1]]
        }
        as.character(v)
      },
      character(1)
    )
  } else {
    as.character(x)
  }
}


# Load aggregated count data
trajectories_limit <- read_rds_or_stop("data/trajectories_limit_counts.rds")
trajectories_limit_li <- read_rds_or_stop("data/trajectories_limit_li_counts.rds")

# Check required columns
req_cols <- c("state_1", "state_2", "state_3", "n")
miss1 <- setdiff(req_cols, names(trajectories_limit))
miss2 <- setdiff(req_cols, names(trajectories_limit_li))
if (length(miss1)) stop(paste("Missing cols in trajectories_limit:", paste(miss1, collapse = ", ")))
if (length(miss2)) stop(paste("Missing cols in trajectories_limit_li:", paste(miss2, collapse = ", ")))

# ======== Common settings ========
default_levels <- c(
  "AD", "AP", "MS",
  "AD; AP", "AD; MS", "AP; MS",
  "AD; AP; MS", "Int/Disc",
  "No Switch"
)

state_colors_default <- c(
  "AD" = "#E41A1C",
  "AP" = "#377EB8",
  "MS" = "#FFD92F",
  "AD; AP" = "#984EA3",
  "AD; MS" = "#FF7F00",
  "AP; MS" = "#4DAF4A",
  "AD; AP; MS" = "#994F00",
  "Int/Disc" = "#999999",
  "No Switch" = "#F5F5F5"
)

state_full_labels_default <- c(
  "AD" = "Antidepressants",
  "AP" = "Antipsychotics",
  "MS" = "Mood Stabilisers",
  "AD; AP" = "Antidepressants & Antipsychotics",
  "AD; MS" = "Antidepressants & Mood Stabilisers",
  "AP; MS" = "Antipsychotics & Mood Stabilisers",
  "AD; AP; MS" = "Antidepressants & Antipsychotics &\n     Mood Stabilisers",
  "Int/Disc" = "Interruption/Discontinuation",
  "No Switch" = "No Switch"
)

# Lithium settings
# Order of states in each step
state_levels_li <- c(
  "Lithium", "Lithium; AD", "Lithium; AP", "Lithium; MS", "Lithium; 2+", 
  "Non-Lithium", "Switch Non-Lithium", "Int/Disc", 
  "No Switch"
)

# Define colors
state_colors_li <- c(
  "Lithium" = "#FFC33B",          # bright red
  "Lithium; AD" = "#E20134",      # purple
  "Lithium; AP" = "#009f81",      # orange
  "Lithium; MS" = "#FF6E3A",      # green
  "Lithium; 2+" = "#994F00",  # brown
  "Non-Lithium" = "#008DF9",          # yellow
  "Switch Non-Lithium" = "#00C2F9",
  "Int/Disc" = "#999999",    # gray
  "No Switch" = "#F5F5F5"
)


# Map abbreviation states
state_full_labels_li <- c(
  "Lithium" = "Lithium", 
  "Lithium; AD" = "Lithium & Antidepressants", 
  "Lithium; AP" = "Lithium & Antipsychotics", 
  "Lithium; MS" = "Lithium & Mood Stabilisers",
  "Lithium; 2+" = "Lithium & 2 or More Classes",
  "Non-Lithium" = "Non-Lithium",
  "Switch Non-Lithium" = "Switch Within Non-Lithium",
  "Int/Disc" = "Interruption/Discontinuation",    
  "No Switch" = "No Switch")


# ======== Static Sankey helpers ========
prep_alluvial_table <- function(tr, state_levels = default_levels) {
  tr <- tr %>%
    mutate(
      state_2 = ifelse(is.na(.data$state_2), "No Switch", .data$state_2),
      state_3 = ifelse(is.na(.data$state_3), "No Switch", .data$state_3)
    ) %>%
    mutate(across(c(.data$state_1, .data$state_2, .data$state_3), as.character))
  # Use existing counts if found
  if (!"n" %in% names(tr)) {
    tr <- count(tr, state_1, state_2, state_3, name = "n")
  } else {
    tr <- tr %>% select(state_1, state_2, state_3, n)
  }
  tr %>%
    mutate(traj_id = row_number()) %>%
    pivot_longer(
      cols = starts_with("state_"),
      names_to = "step",
      values_to = "state"
    ) %>%
    mutate(state = factor(.data$state, levels = state_levels)) %>%
    group_by(.data$step) %>%
    mutate(step_total = sum(.data$n), pct = .data$n / .data$step_total) %>%
    group_by(.data$step, .data$state) %>%
    mutate(
      stratum_pct = sum(.data$n) / first(.data$step_total),
      label = sprintf("%s (%s)", as.character(.data$state),
                      scales::percent(.data$stratum_pct, accuracy = 1))
    ) %>%
    ungroup()
}

alpha_for_hidden <- function(flow, hide_no_switch) {
  if (!hide_no_switch) return(rep(0.8, nrow(flow)))
  is_hidden <- (flow$step == "state_2" & flow$state == "No Switch") |
    (flow$step == "state_3" & flow$state == "No Switch")
  ifelse(is_hidden, 0, 0.8)
}

create_sankey <- function(
    tr,
    state_levels = default_levels,
    state_full_labels = state_full_labels_default,
    state_colors = state_colors_default,
    highlight_state = NULL,
    fade_others = TRUE,
    hide_no_switch = FALSE
) {
  flow <- prep_alluvial_table(tr, state_levels)
  
  init_map <- flow %>%
    dplyr::filter(.data$step == "state_1") %>%
    dplyr::select(.data$traj_id, init = .data$state)
  flow <- dplyr::left_join(flow, init_map, by = "traj_id")
  
  base_alpha <- alpha_for_hidden(flow, hide_no_switch)
  if (!is.null(highlight_state) && highlight_state != "(none)" && isTRUE(fade_others)) {
    is_hi <- as.character(flow$init) == as.character(highlight_state)
    base_alpha <- ifelse(is_hi, base_alpha, pmin(base_alpha, 0.15))
  }
  
  ggplot(
    flow,
    aes(x = .data$step, stratum = .data$state,
        alluvium = .data$traj_id, y = .data$pct,
        fill = .data$state)
  ) +
    ggalluvial::geom_flow(
      stat = "alluvium", lode.guidance = "frontback", knot.pos = 0.5,
      aes(alpha = base_alpha)
    ) +
    scale_alpha_identity(guide = "none") +
    ggalluvial::geom_stratum(width = 0.40, color = "gray40") +
    geom_text(
      stat = "stratum",
      aes(label = .data$label),
      size = 3, fontface = "bold",
      check_overlap = TRUE, color = "black"
    ) +
    scale_x_discrete(labels = c("First", "Second", "Third"),
                     expand = expansion(mult = c(0.01, 0.01))) +
    scale_fill_manual(values = state_colors,
                      breaks = names(state_full_labels),
                      labels = state_full_labels,
                      name = "Drug Combination") +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(y = "Percent of Patients", x = "Treatment Step",
         fill = "Drug Combination") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 13),
          axis.text.x = element_text(size = 11),
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
}

create_sankey_first_state <- function(
    tr,
    state_levels = default_levels,
    state_full_labels = state_full_labels_default,
    state_colors = state_colors_default,
    highlight_state = NULL,
    fade_others = TRUE,
    hide_no_switch = FALSE
) {
  tr <- tr %>%
    mutate(
      state_2 = ifelse(is.na(.data$state_2), "No Switch", .data$state_2),
      state_3 = ifelse(is.na(.data$state_3), "No Switch", .data$state_3)
    ) %>%
    mutate(across(c(.data$state_1, .data$state_2, .data$state_3), as.character))
  if (!"n" %in% names(tr)) {
    tr <- count(tr, state_1, state_2, state_3, name = "n")
  } else {
    tr <- tr %>% select(state_1, state_2, state_3, n)
  }
  flow <- tr %>%
    mutate(traj_id = dplyr::row_number(),
           first_state = .data$state_1) %>%
    tidyr::pivot_longer(cols = tidyselect::starts_with("state_"),
                        names_to = "step", values_to = "state") %>%
    dplyr::mutate(
      state = factor(.data$state, levels = state_levels),
      first_state = factor(.data$first_state, levels = state_levels)
    ) %>%
    dplyr::group_by(.data$step) %>%
    dplyr::mutate(step_total = sum(.data$n), pct = .data$n / .data$step_total) %>%
    dplyr::group_by(.data$step, .data$state) %>%
    dplyr::mutate(
      stratum_pct = sum(.data$n) / dplyr::first(.data$step_total),
      label = sprintf("%s (%s)", as.character(.data$state),
                      scales::percent(.data$stratum_pct, accuracy = 1))
    ) %>%
    dplyr::ungroup()
  
  base_alpha <- alpha_for_hidden(flow, hide_no_switch)
  if (!is.null(highlight_state) && highlight_state != "(none)" && isTRUE(fade_others)) {
    is_hi <- as.character(flow$first_state) == as.character(highlight_state)
    base_alpha <- ifelse(is_hi, base_alpha, pmin(base_alpha, 0.15))
  }
  
  ggplot(
    flow,
    aes(x = .data$step, stratum = .data$state,
        alluvium = .data$traj_id, y = .data$pct,
        fill = .data$first_state)
  ) +
    ggalluvial::geom_flow(
      stat = "alluvium", lode.guidance = "frontback", knot.pos = 0.5,
      aes(alpha = base_alpha)
    ) +
    scale_alpha_identity(guide = "none") +
    ggalluvial::geom_stratum(fill = "gray90", width = 0.40, color = "gray40") +
    geom_text(
      stat = "stratum",
      aes(label = .data$label),
      size = 3, fontface = "bold",
      check_overlap = TRUE, color = "black"
    ) +
    scale_x_discrete(labels = c("First", "Second", "Third"),
                     expand = expansion(mult = c(0.01, 0.01))) +
    scale_fill_manual(values = state_colors,
                      breaks = names(state_full_labels),
                      labels = state_full_labels,
                      name = "Drug Combination") +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(y = "Percent of Patients", x = "Treatment Step",
         fill = "Drug Combination") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 13),
          axis.text.x = element_text(size = 11),
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
}


# ======== Interactive Sankey (networkD3) helpers with colors ========

create_sankey_interactive <- function(
    tr,
    state_levels  = default_levels,
    state_colors  = state_colors_default,
    highlight_state = NULL,
    fade_others     = TRUE,
    hide_no_switch  = FALSE
) {
  # Ensure explicit "No Switch" for missing next states
  tr <- tr %>%
    mutate(
      state_2 = ifelse(is.na(state_2), "No Switch", state_2),
      state_3 = ifelse(is.na(state_3), "No Switch", state_3)
    )
  
  # Optionally hide "No Switch"
  if (hide_no_switch) {
    tr <- tr %>%
      mutate(
        state_2 = ifelse(state_2 == "No Switch", NA_character_, state_2),
        state_3 = ifelse(state_3 == "No Switch", NA_character_, state_3)
      )
  }
  
  # Aggregate if necessary
  if (!"n" %in% names(tr)) {
    tr <- tr %>% count(state_1, state_2, state_3, name = "n")
  } else {
    tr <- tr %>% select(state_1, state_2, state_3, n)
  }
  
  tr <- tr %>% mutate(init_state = state_1)
  
  # ======== Build links ========
  # First -> Second
  links12 <- tr %>%
    filter(!is.na(state_1), !is.na(state_2)) %>%
    group_by(init_state, state_1, state_2) %>%
    summarise(value = sum(n), .groups = "drop") %>%
    mutate(
      source       = gsub("; ", "&", state_1),
      target       = gsub(" ", "_", gsub("; ", "&", paste0(state_2, "_2"))),
      group_init   = gsub("; ", "&", init_state),
      group_current = source,
      # 1 -> 2 tooltip: "First → Second"
      tooltip_base = paste(state_1, "\u2192", state_2)
    )
  
  # Second -> Third
  links23 <- tr %>%
    filter(!is.na(state_2), !is.na(state_3)) %>%
    group_by(init_state, state_2, state_3) %>%
    summarise(value = sum(n), .groups = "drop") %>%
    mutate(
      source       = gsub(" ", "_", gsub("; ", "&", paste0(state_2, "_2"))),
      target       = gsub(" ", "_", gsub("; ", "&", paste0(state_3, "_3"))),
      group_init   = gsub("; ", "&", init_state),
      group_current = source,
      # 2 -> 3 tooltip: "First → Second → Third"
      tooltip_base = paste(init_state, "\u2192", state_2, "\u2192", state_3)
    )
  
  links <- as.data.frame(bind_rows(links12, links23))
  
  # Compute percentage of total patients for each flow
  total_n <- sum(tr$n, na.rm = TRUE)
  links <- links %>%
    mutate(
      pct = if (total_n > 0) round(100 * value / total_n, 1) else NA_real_,
      # Human-readable first state (for tooltips)
      init_label = gsub("&", "; ", group_init)
    )
  
  # Apply highlighting
  links$group <- links$group_current
  if (!is.null(highlight_state) && highlight_state != "(none)" && fade_others) {
    hi_enc <- gsub("; ", "&", highlight_state)
    links$group <- ifelse(links$group_init == hi_enc, links$group, "other")
  }
  
  # ======== Build nodes ========
  nodes <- data.frame(name = unique(c(links$source, links$target)))
  
  nodes <- nodes %>%
    mutate(
      step = case_when(
        grepl("_2$", name) ~ 2L,
        grepl("_3$", name) ~ 3L,
        TRUE               ~ 1L
      ),
      base        = gsub("_", " ", sub("_[0-9]+$", "", name)),
      state       = gsub("&", "; ", base),
      order_state = match(state, state_levels),
      order_state = ifelse(
        is.na(order_state),
        max(order_state, na.rm = TRUE) + 1L,
        order_state
      ),
      order = step * 100 + order_state
    ) %>%
    arrange(order)
  
  nodes$group    <- "my_unique_group"
  links$IDsource <- match(links$source, nodes$name) - 1L
  links$IDtarget <- match(links$target, nodes$name) - 1L
  
  # ======== Colour scale ========
  link_groups <- unique(links$group)
  all_groups  <- unique(c(link_groups, "my_unique_group"))
  
  map_group_to_hex <- function(g) {
    if (g == "my_unique_group") return("#FFFFFF")
    if (g == "other")          return("#DDDDDD")
    
    base_state <- gsub("&", "; ", gsub("_", " ", sub("_[0-9]+$", "", g)))
    col <- unname(state_colors[base_state])
    ifelse(is.na(col), "#999999", col)
  }
  
  palette  <- vapply(all_groups, map_group_to_hex, character(1))
  my_color <- htmlwidgets::JS(sprintf(
    'd3.scaleOrdinal().domain(%s).range(%s)',
    jsonlite::toJSON(all_groups, auto_unbox = TRUE),
    jsonlite::toJSON(unname(palette), auto_unbox = TRUE)
  ))
  
  # ======== Sankey ========
  p <- networkD3::sankeyNetwork(
    Links       = links,
    Nodes       = nodes,
    Source      = "IDsource",
    Target      = "IDtarget",
    Value       = "value",
    NodeID      = "name",
    LinkGroup   = "group",
    NodeGroup   = "group",
    fontSize    = 12,
    nodeWidth   = 30,
    sinksRight  = FALSE,
    colourScale = my_color,
    iterations  = 0
  )
  
  # Zoom/pan + custom tooltips + disable hover on grey flows
  p <- htmlwidgets::onRender(
    p,
    '
  function(el, x) {
    var svg = d3.select(el).select("svg");
    var g   = svg.select("g");
    
    // Helper: clean labels, add space after "_2"/"_3", replace "&" with "; "
    function prettyName(name) {
      if (!name) return "";
      name = name.replace("_2", "_2 ");
      name = name.replace("_3", "_3 ");
      name = name.replace(/&/g, "; ");
      return name;
    }
    
    // Original R link data (post-sankeyNetwork), including tooltip_base
    var linksData = x.links || [];
    
    // Enable zoom & pan
    svg.call(
      d3.zoom().on("zoom", function() {
        g.attr("transform", d3.event.transform);
      })
    );
    
    // Compute total patients as the sum of flows from step 1 -> step 2
    var total = 0;
    d3.select(el).selectAll(".link").each(function(d) {
      var src = d.source.name || "";
      var tgt = d.target.name || "";
      var isStep1to2 = !/_\\d$/.test(src) && /_2$/.test(tgt);
      if (isStep1to2) {
        total += d.value;
      }
    });
    
    // Set tooltips and disable hover for greyed-out flows (links)
    d3.select(el).selectAll(".link").each(function(d, i) {
      var link  = d3.select(this);
      var title = link.select("title");
      if (title.empty()) {
        title = link.append("title");
      }
      
      if (d.group === "other") {
        // Grey flows: no tooltip, no hover
        title.text("");
        link.style("pointer-events", "none");
      } else {
        // Look up the original R row for this link
        var row = (linksData && linksData[i]) ? linksData[i] : {};
        
        // Use tooltip_base from R:
        //  - 1 -> 2 flows: "First → Second"
        //  - 2 -> 3 flows: "First → Second → Third"
        var baseTxt = row.tooltip_base ||
          (prettyName(d.source.name) + " \u2192 " + prettyName(d.target.name));
        
        var txt = baseTxt + "\\n" + d.value;
        
        if (total > 0) {
          var pct = d.value / total * 100;
          txt += " (" + pct.toFixed(1) + "%)";
        }
        
        title.text(txt);
      }
    });
    
    // Add tooltips to nodes (boxes): state name + count + percentage
    d3.select(el).selectAll(".node").each(function(d) {
      var node  = d3.select(this);
      var title = node.select("title");
      if (title.empty()) {
        title = node.append("title");
      }
      
      // Total patients represented by this node:
      // prefer outgoing flows; if none, use incoming flows
      var nodeTotal = 0;
      if (d.sourceLinks && d.sourceLinks.length) {
        d.sourceLinks.forEach(function(l) { nodeTotal += l.value; });
      } else if (d.targetLinks && d.targetLinks.length) {
        d.targetLinks.forEach(function(l) { nodeTotal += l.value; });
      }
      
      var label = prettyName(d.name || "");
      if (total > 0 && nodeTotal > 0) {
        var pctNode = nodeTotal / total * 100;
        label += "\\n" + nodeTotal + " (" + pctNode.toFixed(1) + "%)";
      }
      
      title.text(label);
    });
  }
  '
  )
  
  p
}


# ======== Interactive Sankey (networkD3) colored by first state ========

create_sankey_interactive_first_state <- function(
    tr,
    state_levels = default_levels,
    state_colors = state_colors_default,
    highlight_state = NULL,   # e.g. "AD", "Lithium"
    fade_others = TRUE,
    hide_no_switch = FALSE
) {
  # Ensure explicit "No Switch" for missing next states
  tr <- tr %>%
    mutate(
      state_2 = ifelse(is.na(.data$state_2), "No Switch", .data$state_2),
      state_3 = ifelse(is.na(.data$state_3), "No Switch", .data$state_3)
    )
  
  # Optionally hide "No Switch" in later steps
  if (isTRUE(hide_no_switch)) {
    tr <- tr %>%
      mutate(
        state_2 = ifelse(state_2 == "No Switch", NA_character_, state_2),
        state_3 = ifelse(state_3 == "No Switch", NA_character_, state_3)
      )
  }
  
  # Ensure we have aggregated counts
  if (!"n" %in% names(tr)) {
    tr <- tr %>%
      count(state_1, state_2, state_3, name = "n")
  } else {
    tr <- tr %>% select(state_1, state_2, state_3, n)
  }
  
  # Keep initial state explicitly
  tr <- tr %>%
    mutate(init_state = state_1)
  
  # ================-- LINKS (same structure, but colour group = initial state) ================--
  # First -> Second
  links12 <- tr %>%
    filter(!is.na(state_1), !is.na(state_2)) %>%
    group_by(init_state, state_1, state_2) %>%
    summarise(value = sum(n), .groups = "drop") %>%
    mutate(
      source = state_1,
      source = gsub("; ", "&", source),
      target = paste0(state_2, "_2"),
      target = gsub("; ", "&", target),
      target = gsub(" ", "_", target),
      group_initial = gsub("; ", "&", init_state)
    )
  
  # Second -> Third
  links23 <- tr %>%
    filter(!is.na(state_2), !is.na(state_3)) %>%
    group_by(init_state, state_2, state_3) %>%
    summarise(value = sum(n), .groups = "drop") %>%
    mutate(
      source = paste0(state_2, "_2"),
      source = gsub("; ", "&", source),
      source = gsub(" ", "_", source),
      target = paste0(state_3, "_3"),
      target = gsub("; ", "&", target),
      target = gsub(" ", "_", target),
      group_initial = gsub("; ", "&", init_state)
    )
  
  links <- as.data.frame(bind_rows(links12, links23))
  # Compute percentage of total patients for each flow, rounded to nearest tenth
  total_n <- sum(tr$n, na.rm = TRUE)
  links <- links %>%
    mutate(
      pct = if (total_n > 0) round(100 * value / total_n, 1) else NA_real_,
      # Human-readable first state (for tooltips)
      init_label = gsub("&", "; ", group_initial)
    )
  
  # By default colour by initial state
  links$group <- links$group_initial
  
  # Apply highlighting: grey out non-highlight flows
  if (!is.null(highlight_state) && highlight_state != "(none)" && isTRUE(fade_others)) {
    hi_enc <- gsub("; ", "&", highlight_state)
    links$group <- ifelse(links$group_initial == hi_enc, links$group, "other")
  }
  
  links$group <- factor(links$group)
  
  # ================-- NODES (ordered by step + state_levels) ================--
  nodes <- data.frame(
    name = c(as.character(links$source), as.character(links$target)) %>%
      unique(),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      # which step does this node belong to?
      step = dplyr::case_when(
        grepl("_2$", name) ~ 2L,
        grepl("_3$", name) ~ 3L,
        TRUE               ~ 1L
      ),
      # canonical state label (matching state_levels)
      base  = sub("_[0-9]+$", "", name),     # drop "_2"/"_3"
      base  = gsub("_", " ", base),          # "No_Switch" -> "No Switch"
      state = gsub("&", "; ", base),       # "AD+AP" -> "AD; AP"
      order_state = match(state, state_levels),
      order_state = ifelse(
        is.na(order_state),
        max(order_state, na.rm = TRUE) + 1L,
        order_state
      ),
      order = step * 100 + order_state
    ) %>%
    arrange(order, state, name)
  
  # nodes in a single group; flows carry the colour meaning
  nodes$group <- as.factor("my_unique_group")
  
  # ================-- ID mapping AFTER ordering ================--
  links$IDsource <- match(links$source, nodes$name) - 1L
  links$IDtarget <- match(links$target, nodes$name) - 1L
  
  # ================-- COLOUR SCALE from state_colors (+ highlight grey) ================--
  link_groups <- levels(links$group)
  all_groups  <- unique(c(link_groups, "my_unique_group"))
  
  map_group_to_hex <- function(g) {
    if (g == "my_unique_group") {
      return("#FFFFFF")  # node fill
    }
    if (g == "other") {
      return("#DDDDDD")  # faded / non-highlight flows
    }
    
    # Decode group back to initial state label
    base <- sub("_[0-9]+$", "", g)        # (mostly irrelevant here; groups are initial)
    base <- gsub("_", " ", base)
    base_state <- gsub("&", "; ", base) # Lithium+AD -> Lithium; AD (if ever used)
    
    # For initial-state colouring, we only care about the *initial* state;
    # if state_colors is keyed by initial categories (Lithium, Non-Lithium, etc),
    # base_state will match those.
    col <- unname(state_colors[base_state])
    if (is.na(col) || length(col) == 0) {
      return("#999999")
    }
    col
  }
  
  palette <- vapply(all_groups, map_group_to_hex, character(1))
  
  domain_js <- jsonlite::toJSON(all_groups, auto_unbox = TRUE)
  range_js  <- jsonlite::toJSON(unname(palette), auto_unbox = TRUE)
  
  my_color <- htmlwidgets::JS(
    sprintf(
      'd3.scaleOrdinal()
         .domain(%s)
         .range(%s)',
      domain_js,
      range_js
    )
  )
  
  # ================-- Build sankey widget ================--
  p <- networkD3::sankeyNetwork(
    Links       = links,
    Nodes       = nodes,
    Source      = "IDsource",
    Target      = "IDtarget",
    Value       = "value",
    NodeID      = "name",
    LinkGroup   = "group",
    NodeGroup   = "group",
    fontSize    = 12,
    nodeWidth   = 30,
    sinksRight  = FALSE,
    colourScale = my_color,
    iterations  = 0 # disable automatic reordering
  )
  
  p <- htmlwidgets::onRender(
    p,
    '
  function(el, x) {
    var svg = d3.select(el).select("svg");
    var g   = svg.select("g");
    
    // Helper: clean labels, add space after "_2"/"_3"
    function prettyName(name) {
      if (!name) return "";
      // keep the step suffixes readable
      name = name.replace("_2", "_2 ");
      name = name.replace("_3", "_3 ");
      return name;
    }
    
    // Enable zoom & pan
    svg.call(
      d3.zoom().on("zoom", function() {
        g.attr("transform", d3.event.transform);
      })
    );
    
    // First, compute total patients as the sum of flows
    // from step 1 (no "_2"/"_3" suffix) to step 2 ("_2" suffix)
    var total = 0;
    d3.select(el).selectAll(".link").each(function(d) {
      var src = d.source.name || "";
      var tgt = d.target.name || "";
      var isStep1to2 = !/_\\d$/.test(src) && /_2$/.test(tgt);
      if (isStep1to2) {
        total += d.value;
      }
    });
    
    // Now set tooltips and disable hover for greyed-out flows
    d3.select(el).selectAll(".link").each(function(d) {
      var link  = d3.select(this);
      var title = link.select("title");
      
      if (d.group === "other") {
        // Grey flows: no tooltip, no hover
        title.text("");
        link.style("pointer-events", "none");
      } else {
        var txt = d.source.name + " \u2192 " + d.target.name + "\\n" + d.value;
        
        if (total > 0) {
          var pct = d.value / total * 100;
          txt += " (" + pct.toFixed(1) + "%)";
        }
        
        title.text(txt);
      }
    });
    
    // Add tooltips to nodes (boxes): state name + count + percentage
    d3.select(el).selectAll(".node").each(function(d) {
      var node  = d3.select(this);
      var title = node.select("title");
      if (title.empty()) {
        title = node.append("title");
      }
      
      // Total patients represented by this node:
      // prefer outgoing flows; if none, use incoming flows
      var nodeTotal = 0;
      if (d.sourceLinks && d.sourceLinks.length) {
        d.sourceLinks.forEach(function(l) { nodeTotal += l.value; });
      } else if (d.targetLinks && d.targetLinks.length) {
        d.targetLinks.forEach(function(l) { nodeTotal += l.value; });
      }
      
      var label = prettyName(d.name || "");
      if (total > 0 && nodeTotal > 0) {
        var pctNode = nodeTotal / total * 100;
        label += "\\n" + nodeTotal + " (" + pctNode.toFixed(1) + "%)";
      }
      
      title.text(label);
    });
  }
  '
  )
  
  p
  
}



# ======== Heatmap first-to-second transition ========
make_transition_heatmap <- function(tr, axis_order, title = NULL, pct_digits = 1) {
  tr <- tr %>%
    mutate(state_2 = ifelse(is.na(.data$state_2), "No Switch", .data$state_2))
  # Get counts, summing over state_3
  if (!"n" %in% names(tr)) {
    tr <- count(tr, state_1, state_2, name = "n")
  } else {
    tr <- tr %>% group_by(., state_1, state_2) %>%
      summarise(n = sum(n, na.rm = TRUE), .groups = "drop")
  }
  
  # State 1 percentages + total per state 1
  df <- tr %>%
    group_by(.data$state_1) %>%
    mutate(percentage = (n / sum(n)) * 100, total_n = sum(n)) %>%
    ungroup()
  # Grand total across all state 1 groups
  grand_total <- df %>%
    distinct(state_1, total_n) %>%
    summarise(gt = sum(total_n, na.rm = TRUE), .groups = "drop") %>%
    pull(gt)
  # Row labels that include n and % of grand total
  sample_sizes <- df %>%
    distinct(state_1, total_n) %>%
    mutate(
      state_1 = factor(state_1, levels = axis_order),
      pct_of_total = if (grand_total > 0) (total_n / grand_total) * 100 else NA_real_,
      label_simple = paste0(
        state_1, "\n",
        "n=", format(total_n, big.mark = ","), " (",
        sprintf(paste0("%.", pct_digits, "f%%"), pct_of_total), ")"
      )
    ) %>%
    arrange(state_1)
  
  df_final <- df %>%
    left_join(sample_sizes %>% select(state_1, label_simple), by = "state_1") %>%
    mutate(
      state_2 = factor(.data$state_2, levels = axis_order),
      label_simple = factor(.data$label_simple, levels = sample_sizes$label_simple)
    )
  # Dynamic max for color scale
  maxp <- max(df_final$percentage, na.rm = TRUE)
  color_max <- ceiling(maxp / 10) * 10 + 10
  
  ggplot(df_final, aes(x = .data$state_2, y = .data$label_simple,
                       fill = .data$percentage)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(
      aes(label = sprintf(paste0("%.", pct_digits, "f%%"), percentage)),
      color = "black", size = 3.5, fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = "#f7fbff", mid = "#4292c6", high = "#08519c",
      midpoint = color_max / 2, limits = c(0, color_max),
      name = "Percentage (%)"
    ) +
    labs(x = "Second-Line Treatment", y = "First-Line Treatment", title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                 face = "bold", size = 10),
      axis.text.y = element_text(face = "bold", size = 8, lineheight = 0.9),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    ) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                 barwidth = 15, barheight = 0.8))
}

#================== Heatmap second and third line (overall) ====================

# Ensure explicit labels for missing states
trajectories_limit_23 <- trajectories_limit %>%
  mutate(
    state_2 = ifelse(is.na(state_2), "No Switch", state_2),
    state_3 = ifelse(is.na(state_3), "No Switch", state_3)
  )

# Aggregate second → third, then compute within–state_2 percentages
transition_matrix_23 <- trajectories_limit_23 %>%
  group_by(state_2, state_3) %>%
  summarise(n = sum(n, na.rm = TRUE), .groups = "drop") %>%
  group_by(state_2) %>%
  mutate(percentage = (n / sum(n)) * 100, total_n = sum(n)) %>%
  ungroup()

# Order for both axes
axis_order_23 <- c("AD", "AP", "MS", "AD; AP", "AD; MS", "AP; MS",
                   "AD; AP; MS", "No Switch", "Int/Disc")

# Row labels: include n and % of grand total (based on second-line totals)
grand_total_23 <- sum(unique(transition_matrix_23$total_n), na.rm = TRUE)

sample_sizes_23 <- transition_matrix_23 %>%
  distinct(state_2, total_n) %>%
  mutate(
    state_2      = factor(state_2, levels = axis_order_23),
    pct_of_total = if (grand_total_23 > 0) (total_n / grand_total_23) * 100 else NA_real_,
    label_simple = paste0(
      state_2, "\n",
      "n=", format(total_n, big.mark = ","), " (",
      sprintf("%.1f%%", pct_of_total), ")"
    )
  ) %>%
  arrange(state_2)

transition_matrix_final_23 <- transition_matrix_23 %>%
  left_join(sample_sizes_23 %>% select(state_2, label_simple), by = "state_2") %>%
  mutate(
    state_3      = factor(state_3, levels = axis_order_23),
    label_simple = factor(label_simple, levels = sample_sizes_23$label_simple)
  )

# Hide "No Switch" row only at plot time (keep percentages unchanged) 
transition_matrix_plot <- transition_matrix_final_23 %>%
  dplyr::filter(state_2 != "No Switch") %>%
  droplevels()  # drop the now-empty y-factor level

# Dynamic max for color scale
max_percentage_23 <- max(transition_matrix_plot$percentage, na.rm = TRUE)
color_max_23 <- ceiling(max_percentage_23 / 5) * 5 + 5

# Plot: Fixed row heights (Second → Third)
second_to_third <- ggplot(
  transition_matrix_plot,
  aes(x = state_3, y = label_simple, fill = percentage)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = sprintf("%.1f%%", percentage)),
    color = "black", size = 3.5, fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "#f7fbff", mid = "#4292c6", high = "#08519c",
    midpoint = color_max_23 / 2, limits = c(0, color_max_23),
    name = "Percentage (%)"
  ) +
  labs(
    x = "Third-Line Treatment",
    y = "Second-Line Treatment"
    # title = "Second → Third (Fixed Row Heights)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(
      angle = 45, hjust = 1, vjust = 1,
      face = "bold", size = 10
    ),
    axis.text.y  = element_text(face = "bold", size = 8, lineheight = 0.9),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
    panel.grid   = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    plot.title   = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 10))
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = 15, barheight = 0.8
    )
  )

# ================ Second → Third transitions (Lithium-focused) ================

trajectories_limit_li_23 <- trajectories_limit_li %>%
  mutate(
    state_2 = ifelse(is.na(state_2), "No Switch", state_2),
    state_3 = ifelse(is.na(state_3), "No Switch", state_3)
  )

transition_matrix_23_li <- trajectories_limit_li_23 %>%
  group_by(state_2, state_3) %>%
  summarise(n = sum(n, na.rm = TRUE), .groups = "drop") %>%
  group_by(state_2) %>%
  mutate(percentage = (n / sum(n)) * 100, total_n = sum(n)) %>%
  ungroup()

axis_order_li_23 <- c(
  "Lithium", "Lithium; AD", "Lithium; AP", "Lithium; MS", "Lithium; 2+",
  "Non-Lithium", "Switch Non-Lithium", "No Switch", "Int/Disc"
)

grand_total_23_li <- sum(unique(transition_matrix_23_li$total_n), na.rm = TRUE)

sample_sizes_23_li <- transition_matrix_23_li %>%
  distinct(state_2, total_n) %>%
  mutate(
    state_2      = factor(state_2, levels = axis_order_li_23),
    pct_of_total = if (grand_total_23_li > 0) (total_n / grand_total_23_li) * 100 else NA_real_,
    label_simple = paste0(
      state_2, "\n",
      "n=", format(total_n, big.mark = ","), " (",
      sprintf("%.1f%%", pct_of_total), ")"
    )
  ) %>%
  arrange(state_2)

transition_matrix_final_23_li <- transition_matrix_23_li %>%
  left_join(sample_sizes_23_li %>% select(state_2, label_simple), by = "state_2") %>%
  mutate(
    state_3      = factor(state_3, levels = axis_order_li_23),
    label_simple = factor(label_simple, levels = sample_sizes_23_li$label_simple)
  )

# Hide "No Switch" row only at plot time (keep percentages unchanged) 
transition_matrix_plot_li <- transition_matrix_final_23_li %>%
  dplyr::filter(state_2 != "No Switch") %>%
  droplevels()  # drop the now-empty y-factor level

max_percentage_23_li <- max(transition_matrix_plot_li$percentage, na.rm = TRUE)
color_max_23_li <- ceiling(max_percentage_23_li / 5) * 5 + 5

second_to_third_li <- ggplot(
  transition_matrix_plot_li,
  aes(x = state_3, y = label_simple, fill = percentage)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = sprintf("%.1f%%", percentage)),
    color = "black", size = 3.5, fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "#f7fbff", mid = "#4292c6", high = "#08519c",
    midpoint = color_max_23_li / 2, limits = c(0, color_max_23_li),
    name = "Percentage (%)"
  ) +
  labs(
    x = "Third-Line Treatment",
    y = "Second-Line Treatment"
    # title = "Second → Third (Lithium, Fixed Row Heights)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(
      angle = 45, hjust = 1, vjust = 1,
      face = "bold", size = 10
    ),
    axis.text.y  = element_text(face = "bold", size = 8, lineheight = 0.9),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
    panel.grid   = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    plot.title   = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 10))
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = 15, barheight = 0.8
    )
  )




# ======== UI ========

ui <- fluidPage(
  titlePanel("Psychiatric prescribing patterns in patients with newly diagnosed bipolar disorder in the UK: 2000 to 2022"),
  wellPanel(
    id = "overview-box",
    style = "background:#f9fbff; border-color:#e3eefc;",
    tags$p(
      "This interactive dashboard summarizes psychiatric prescribing 
      patterns for patients newly diagnosed with bipolar disorder in the UK 
      (2000-2022). Prescriptions are followed from for one year from diagnosis."
    ),
    tags$p(
      HTML("The <strong>Trajectory Controls</strong> sidebar allows you to 
           highlight specific treatments, switch colour modes, or focus on 
           lithium-specific regimens.")
    ),
    tags$p(
      HTML("The <strong>Prescribing patterns visualisations</strong> tab provides static plots that can be downloaded:")
    ),
    tags$ul(
      style = "margin:10px 10px 10px 22px;",
      tags$li(HTML("<strong>Overall prescribing patterns:</strong> Sankey diagram showing up to three treatment changes, with flow height representing the proportion of patients following each path.")),
      tags$li(HTML("<strong>Filtered Sankey diagram:</strong> Focuses on patients with a specific initial treatment state.")),
      tags$li(HTML("<strong>Heat map:</strong> Summarizes the proportion of patients moving from first- to second-line combinations.")),
    ),
    tags$p(
      HTML("The <strong>Interactive Sankey diagram</strong> tab provides an interactive plot where you can hover over each flow to see the percentage of patients in that flow.")
    ),
    tags$p(
      HTML(paste0("Abbreviations used in treatment states: ",
                  "<strong>AD</strong> = Antidepressants; ",
                  "<strong>AP</strong> = Antipsychotics; ",
                  "<strong>MS</strong> = Mood Stabilisers; ",
                  "<strong>AD; AP</strong> = Antidepressants &amp; Antipsychotics; ",
                  "<strong>AD; MS</strong> = Antidepressants &amp; Mood Stabilisers; ",
                  "<strong>AP; MS</strong> = Antipsychotics &amp; Mood Stabilisers; ",
                  "<strong>AD; AP; MS</strong> = Antidepressants &amp; Antipsychotics &amp; Mood Stabilisers; ",
                  "<strong>Int/Disc</strong> = Interruption/Discontinuation; ",
                  "<strong>No Switch</strong> = No change in combination between steps."))
    )
    # tags$ul(
    #   style = "margin:6px 0 0 22px;",
    #   tags$li(HTML("<strong>AD</strong> = Antidepressants")),
    #   tags$li(HTML("<strong>AP</strong>: Antipsychotics")),
    #   tags$li(HTML("<strong>MS</strong>: Mood Stabilisers")),
    #   tags$li(HTML("<strong>AD; AP</strong>: Antidepressants &amp; Antipsychotics")),
    #   tags$li(HTML("<strong>AD; MS</strong>: Antidepressants &amp; Mood Stabilisers")),
    #   tags$li(HTML("<strong>AP; MS</strong>: Antipsychotics &amp; Mood Stabilisers")),
    #   tags$li(HTML("<strong>AD; AP; MS</strong>: Antidepressants &amp; Antipsychotics &amp; Mood Stabilisers")),
    #   tags$li(HTML("<strong>Int/Disc</strong>: Interruption/Discontinuation")),
    #   tags$li(HTML("<strong>No Switch</strong>: No change in combination between steps"))
    # )
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Trajectory Controls"),
      tags$div(
        class = "help-block",
        tags$p("To explore different prescribing pathways, use the options below to:"),
        tags$ol(
          style = "margin: 6px 0 10px 10px; padding-left: 14px;",
          tags$li("Pick an initial state to highlight specific trajectories."),
          tags$li("Toggle how flows are coloured."),
          tags$li("Optionally focus on lithium.")
        )
      ),
      
      tags$hr(),
      h5("Highlight and Filter by Initial State"),
      helpText(paste0("Choose a first-line treatment state to highlight. The overall ", 
                      "Sankey will dim other starting states, and a filtered Sankey ", 
                      "below will show only trajectories that begin with your selection.")),
      selectInput(
        "highlight_state",
        NULL,
        choices = c("(none)", default_levels),
        selected = "(none)"
      ),
      
      tags$hr(),
      h5("Colour by Initial State"),
      helpText(paste0("Check the box below to colour the flows by the initial ", 
                      "prescribing state and compare how different starting ", 
                      "regimens evolve. By default, the flows are coloured ",
                      "by the current treatment at each step")),
      checkboxInput("color_by_first", "Color flows by first state", FALSE),
      
      tags$hr(),
      h5("Focus on Lithium"),
      helpText(paste0("Click here to switch the Sankey diagrams to focus on ", 
                      "Lithium-specific prescribing categories. By default, ",
                      "overall prescribing categories are shown.")),
      checkboxInput("lithium_mode", "Lithium-focused", FALSE),
      # checkboxInput("fade_unselected", "Fade other pathways (static)", TRUE),
      # checkboxInput("hide_no_switch", "Hide 'No Switch' flows (all views)", FALSE),
      
      # tags$hr(),
      # checkboxInput("mosaic_mode", "Heatmap → Mosaic layout", FALSE),
      # sliderInput("min_prop", "Mosaic minimum segment proportion",
      #             min = 0, max = 0.05, value = 0.01, step = 0.002, ticks = FALSE),
      tags$hr(),
      
      ### Downloading alluvial and heatmap plots
      
      # Choose which static alluvial plot to download
      selectInput(
        "which_sankey_download", "Sankey to download:",
        choices = c(
          "Main Sankey (top plot)" = "main",
          "Filtered Sankey (by selected initial state)" = "filtered"
        ),
        selected = "main"
      ),
      downloadButton("dl_alluvial", label = HTML("Download<br>Sankey (PNG)")),
      
      # Choose which heatmap to download
      selectInput(
        "which_heatmap_download", "Heatmap to download:",
        choices = c(
          "First → Second (overall)"          = "fs_overall",
          "First → Second (lithium-focused)"  = "fs_lithium",
          "Second → Third (overall)"          = "st_overall",
          "Second → Third (lithium-focused)"  = "st_lithium"
        ),
        selected = "fs_overall"
      ),
      downloadButton("dl_heatmap", label = HTML("Download<br>Heatmap (PNG)")),
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "Prescribing patterns visualisations",
          uiOutput("sankey_title"),
          uiOutput("sankey_description"),
          plotOutput("sankey_plot", height = 520),
          tags$hr(),
          h4("Filtered Sankey highlighting selected initial state"),
          uiOutput("selected_caption"),
          uiOutput("filtered_sankey_description"),
          conditionalPanel(
            "input.highlight_state != '(none)'",
            plotOutput("sankey_filtered", height = 520)
          ),
          tags$hr(),
          uiOutput("heatmap_title"),
          uiOutput("heatmap_description"),
          
          # First → Second heatmaps
          conditionalPanel(
            "input.lithium_mode == false",
            plotOutput("heat_overall", height = 520)
          ),
          conditionalPanel(
            "input.lithium_mode == true",
            plotOutput("heat_lithium", height = 520)
          ),
          
          tags$br(),
          h5("Second- to third-line prescribing transitions"),
          # Second → Third heatmaps
          conditionalPanel(
            "input.lithium_mode == false",
            plotOutput("heat_overall_23", height = 520)
          ),
          conditionalPanel(
            "input.lithium_mode == true",
            plotOutput("heat_lithium_23", height = 520)
          )
          
        ),
        
        # tabPanel(
        #   "Sankey – Color by First State",
        #   plotOutput("sankey_first_plot", height = 600)
        # )
        tabPanel(
          "Interactive Sankey diagram",
          uiOutput("sankey_interactive_caption"),
          uiOutput("sankey_interactive_description"),
          sankeyNetworkOutput("sankey_interactive", height = "640px")
        )
      )
    )
  ),
  tags$hr(),
  tags$p(
    style = "font-size:13px; line-height:1.3; margin-top:-10px; color:#555;",
    "Note: This demonstration app uses aggregated, de-identified data only. 
   No individual patient records are included."
  ),
  tags$div(
    style = "font-size:12px; color:#555; line-height:1.3; margin-bottom:30px;",
    HTML(paste0(
      "<strong>Cite this app:</strong> <br>",
      "Wu, S.M., Barnett J.F., Launders, N., Bramon, E., Osborn, D.P.J., Hayes, J.F., & Richards-Belle, A. (2025). ",
      "<em>Psychiatric prescribing patterns in patients with newly diagnosed bipolar disorder in the UK: 2000–2022</em> ",
      "(Version ", "1.0.0", ") [Shiny web application]. ",
      paste0("URL: https://dop-mhds.shinyapps.io/bpd-prescribing-uk-shiny/.")
    ))
  )
)

# ======== Server ========
server <- function(input, output, session) {
  
  current_traj <- reactive({
    tr <- if (isTRUE(input$lithium_mode)) trajectories_limit_li else trajectories_limit
    if (isTRUE(input$filter_AD) && "state_1" %in% names(tr)) {
      tr <- dplyr::filter(tr, .data$state_1 == "AD")
    }
    tr
  })
  
  # Keep the highlight dropdown synced to available first states in current data/mode
  observe({
    tr <- current_traj()
    choices <- c("(none)", sort(unique(as.character(na.omit(tr$state_1)))))
    sel <- if (as.character(input$highlight_state) %in% choices) as.character(input$highlight_state) else "(none)"
    updateSelectInput(session, "highlight_state", choices = choices, selected = sel)
  })
  
  # Sankey diagram header
  output$sankey_title <- renderUI({
    if (isTRUE(input$color_by_first)) {
      if (isTRUE(input$lithium_mode)) {
        h4("Lithium prescribing patterns colored by initial state")
      } else {
        h4("Overall prescribing patterns colored by initial state")
      }
    } else {
      if (isTRUE(input$lithium_mode)) {
        h4("Lithium prescribing patterns")
      } else {
        h4("Overall prescribing patterns")
      }
    }
  })
  
  output$sankey_description <- renderUI({
    color_by_first <- isTRUE(input$color_by_first)
    lithium_mode   <- isTRUE(input$lithium_mode)
    highlight      <- as.character(input$highlight_state)
    
    # Construct description dynamically
    desc <- ""
    
    if (lithium_mode) {
      desc <- paste0(
        "This Sankey diagram summarizes treatment transitions in ",
        strong("lithium-focused trajectories"), 
        " among patients newly diagnosed with bipolar disorder."
      )
    } else {
      desc <- paste0(
        "This Sankey diagram illustrates overall ",
        strong("psychiatric medication switching patterns"),
        " among patients newly diagnosed with bipolar disorder."
      )
    }
    
    if (color_by_first) {
      desc <- paste0(
        desc,
        " Flows are colored by the patient's ",
        strong("initial treatment state"), 
        ", highlighting how different starting regimens evolve over time."
      )
    } else {
      desc <- paste0(
        desc,
        " Flows are colored by the ",
        strong("current treatment state"),
        " at each step."
      )
    }
    
    if (!is.null(highlight) && highlight != "(none)") {
      desc <- paste0(
        desc,
        " The plot currently ",
        strong("highlights patients starting in "), strong(highlight), "."
      )
    }
    
    # Return formatted HTML paragraph
    HTML(paste0("<p style='font-size:14px; line-height:1.4;'>", desc, "</p>"))
  })
  
  
  # ======== Static Sankey – Overall (fades others if a highlight is chosen)
  output$sankey_plot <- renderPlot({
    tr <- current_traj()
    highlight <- as.character(input$highlight_state)
    fade <- TRUE  # fade unselected
    hide_ns <- FALSE # do not hide no switch
    color_by_first <- isTRUE(input$color_by_first)  # color by initial state
    
    if (isTRUE(input$lithium_mode)) { # Lithium focus
      if (color_by_first) { # Color by initial
        create_sankey_first_state(
          tr,
          state_levels = state_levels_li,
          state_full_labels = state_full_labels_li,
          state_colors = state_colors_li,
          highlight_state = as.character(input$highlight_state),
          fade_others = fade,
          hide_no_switch = hide_ns
        )
      } else { # Do not color by initial
        create_sankey(
          tr,
          state_levels = state_levels_li,
          state_full_labels = state_full_labels_li,
          state_colors = state_colors_li,
          highlight_state = as.character(input$highlight_state),
          fade_others = fade,
          hide_no_switch = hide_ns
        )
      }
      
    } else { # No lithium focus
      if (color_by_first) { # Color by initial state
        create_sankey_first_state(
          tr,
          highlight_state = as.character(input$highlight_state),
          fade_others = fade,
          hide_no_switch = hide_ns
        )
      } else { # Do not color by initial state
        create_sankey(
          tr,
          highlight_state = as.character(input$highlight_state),
          fade_others = fade,
          hide_no_switch = hide_ns
        )
      }
    }
  })
  
  # ======== Filtered Sankey (subset where state_1 == selected highlight)
  output$selected_caption <- renderUI({
    if (is.null(input$highlight_state) || input$highlight_state == "(none)") {
      HTML("<em>Select an initial state to highlight using the panel on the left.</em>")
    } else {
      HTML(paste0("<strong>Selected initial state:</strong> ", input$highlight_state))
    }
  })
  
  output$filtered_sankey_description <- renderUI({
    hs <- as.character(input$highlight_state)
    lithium_mode   <- isTRUE(input$lithium_mode)
    color_by_first <- isTRUE(input$color_by_first)
    
    if (is.null(hs) || hs == "(none)") return(NULL)
    
    # Base description
    desc <- if (lithium_mode) {
      paste0(
        "This filtered Sankey shows the treatment transitions among patients whose ",
        "initial combination included ", strong(hs),
        ", within the lithium-focused analysis."
      )
    } else {
      paste0(
        "This filtered Sankey illustrates the treatment pathways of patients who ",
        "began with ", strong(hs), " as their first-line therapy."
      )
    }
    
    # Color scheme reference
    desc <- paste0(
      desc, " Flows are colored by the ",
      if (color_by_first) "initial treatment state, " else "current treatment state, ",
      "and their thickness reflects the proportion of patients following each trajectory. ",
      "Percentages are updated to only consider patients with the highlighted initial state."
    )
    
    HTML(sprintf("<p style='font-size:14px; line-height:1.4;'>%s</p>", desc))
  })
  
  
  output$sankey_filtered <- renderPlot({
    
    hs <- as.character(input$highlight_state)
    if (is.null(hs) || hs == "(none)") return(NULL)
    
    tr_all <- current_traj()
    tr <- dplyr::filter(tr_all, as.character(.data$state_1) == hs)
    if (nrow(tr) == 0) return(NULL)
    
    if (isTRUE(input$lithium_mode)) {
      create_sankey(
        tr,
        state_levels = state_levels_li,
        state_full_labels = state_full_labels_li,
        state_colors = state_colors_li,
        highlight_state = NULL,                 # no fading inside the subset
        fade_others = FALSE,
        hide_no_switch = FALSE
      )
    } else {
      create_sankey(
        tr,
        highlight_state = NULL,                 # no fading inside the subset
        fade_others = FALSE,
        hide_no_switch = FALSE
      )
    }
  })
  
  
  # ======== Heatmap 
  output$heatmap_title <- renderUI({
    if (isTRUE(input$lithium_mode)) {
      h4("Heat map of prescribing transitions with lithium focus")
    } else {
      h4("Heat map of prescribing transitions")
    }
  })
  
  output$heatmap_description <- renderUI({
    lithium_mode <- isTRUE(input$lithium_mode)
    
    if (lithium_mode) {
      desc <- paste0(
        "This heatmap summarizes the proportion of patients transitioning from each ",
        "first-line lithium-specific treatment to each second-line combination. ",
        "Darker tiles indicate more frequent transitions, while lighter tiles indicate less common changes. ",
        "Below, a second heatmap shows transitions from the second- to the third-line lithium-specific combinations, ",
        "using the same color scale and labelling."
      )
    } else {
      desc <- paste0(
        "This heatmap depicts transitions from the initial to the second-line treatment combination ",
        "among patients newly diagnosed with bipolar disorder. ",
        "The color intensity represents the proportion of patients switching between specific treatment categories. ",
        "Below, a second heatmap summarizes transitions from the second- to the third-line treatment combinations ",
        "among the same overall prescribing categories."
      )
    }
    
    HTML(paste0("<p style='font-size:14px; line-height:1.4;'>", desc, "</p>"))
  })
  
  
  output$heat_overall <- renderPlot({
    tr <- current_traj()
    axis_order <- c("AD", "AP", "MS", "AD; AP", "AD; MS", "AP; MS",
                    "AD; AP; MS", "No Switch", "Int/Disc")
    make_transition_heatmap(
      tr, axis_order = axis_order,
      title = "Initial \u2192 Second-line (Heat map)"
    )
  })
  
  # Lithium heatmap (only shown in UI when lithium_mode is TRUE)
  output$heat_lithium <- renderPlot({
    axis_order <- c(
      "Lithium", "Lithium; AD", "Lithium; AP",
      "Lithium; MS", "Lithium; 2+",
      "Non-Lithium", "Switch Non-Lithium", "No Switch", "Int/Disc"
    )
    make_transition_heatmap(
      trajectories_limit_li, axis_order,
      title = "Initial \u2192 Second-line (Lithium)"
    )
  })
  
  # Second → Third heatmap (overall)
  output$heat_overall_23 <- renderPlot({
    second_to_third
  })
  
  # Second → Third heatmap (Lithium-focused)
  output$heat_lithium_23 <- renderPlot({
    second_to_third_li
  })
  
  # ======== Interactive Sankey (percent tooltips; hide flows visually)
  output$sankey_interactive_caption <- renderUI({
    hs <- as.character(input$highlight_state)
    mode_label <- if (isTRUE(input$lithium_mode)) " (Lithium-focused)" else ""
    
    if (is.null(hs) || hs == "(none)") {
      txt <- paste0("No initial state selected for highlighting", mode_label,
                    ". Showing all pathways.")
    } else {
      txt <- paste0("Highlighting trajectories with initial state: <strong>",
                    hs, "</strong>", mode_label, ".")
    }
    
    HTML(paste0("<p><em>", txt, "</em></p>"))
  })
  
  output$sankey_interactive_description <- renderUI({
    color_by_first <- isTRUE(input$color_by_first)
    lithium_mode   <- isTRUE(input$lithium_mode)
    highlight      <- as.character(input$highlight_state)
    
    # Construct description dynamically
    desc <- ""
    
    if (lithium_mode) {
      desc <- paste0(
        "This interactive Sankey diagram summarizes treatment transitions in ",
        strong("lithium-focused trajectories"), 
        " among patients newly diagnosed with bipolar disorder."
      )
    } else {
      desc <- paste0(
        "This Sankey diagram illustrates overall ",
        strong("psychiatric medication switching patterns"),
        " among patients newly diagnosed with bipolar disorder. " 
      )
    }
    
    if (color_by_first) {
      desc <- paste0(
        desc,
        " Flows are colored by the patient's ",
        strong("initial treatment state"), 
        ", highlighting how different starting regimens evolve over time."
      )
    } else {
      desc <- paste0(
        desc,
        " Flows are colored by the ",
        strong("current treatment state"),
        " at each step."
      )
    }
    
    if (!is.null(highlight) && highlight != "(none)") {
      desc <- paste0(
        desc,
        " The plot currently ",
        strong("highlights patients starting in "), strong(highlight), ". "
      )
    }
    
    desc <- paste0(desc, "<br><br>Hover over the flows to see the number and ", 
                   "percentage of all patients following each path. ",
                   "You can also zoom in and out of the plot and move the boxes around.")
    
    # Return formatted HTML paragraph
    HTML(paste0("<p style='font-size:14px; line-height:1.4;'>", desc, "</p>"))
  })
  
  
  output$sankey_interactive <- networkD3::renderSankeyNetwork({
    tr <- current_traj()
    highlight <- as.character(input$highlight_state)
    fade <- TRUE  # fade unselected
    hide_ns <- FALSE # do not hide no switch
    color_by_first <- isTRUE(input$color_by_first)  # color by initial state
    
    # Choose color palette depending on lithium mode
    if (isTRUE(input$lithium_mode)) {
      colors <- state_colors_li
      levels <- state_levels_li
    } else {
      colors <- state_colors_default
      levels <- default_levels
    }
    
    if (color_by_first) { # Color by initial state
      create_sankey_interactive_first_state(
        tr = tr,
        state_levels = levels,
        state_colors = colors,
        highlight_state = highlight,
        fade_others = fade,
        hide_no_switch = hide_ns
      )
    } else { # Color by current state
      create_sankey_interactive(
        tr = tr,
        state_levels = levels,
        state_colors = colors,
        highlight_state = highlight,
        fade_others = fade,
        hide_no_switch = hide_ns
      )
    }
  })
  
  # ======== Downloads
  output$dl_alluvial <- downloadHandler(
    filename = function() sprintf("sankey_%s.png", Sys.Date()),
    content  = function(file) {
      tr <- current_traj()
      hs <- as.character(input$highlight_state)
      fade <- TRUE
      hide_ns <- FALSE
      color_by_first <- isTRUE(input$color_by_first)
      li <- isTRUE(input$lithium_mode)
      
      # If user chose "filtered" and a valid highlight is selected
      if (identical(input$which_sankey_download, "filtered") &&
          !is.null(hs) && hs != "(none)") {
        
        tr_filt <- dplyr::filter(tr, as.character(.data$state_1) == hs)
        
        # Fallback: if no data for that first state, use full data
        if (nrow(tr_filt) == 0) tr_filt <- tr
        
        if (li) {
          p <- create_sankey(
            tr_filt,
            state_levels = state_levels_li,
            state_full_labels = state_full_labels_li,
            state_colors = state_colors_li,
            highlight_state = NULL,   # no fading inside subset
            fade_others = FALSE,
            hide_no_switch = hide_ns
          )
        } else {
          p <- create_sankey(
            tr_filt,
            highlight_state = NULL,   # no fading inside subset
            fade_others = FALSE,
            hide_no_switch = hide_ns
          )
        }
        
      } else {
        # Main Sankey: mirror the logic used in output$sankey_plot
        if (li) { # Lithium focus
          if (color_by_first) {
            p <- create_sankey_first_state(
              tr,
              state_levels = state_levels_li,
              state_full_labels = state_full_labels_li,
              state_colors = state_colors_li,
              highlight_state = hs,
              fade_others = fade,
              hide_no_switch = hide_ns
            )
          } else {
            p <- create_sankey(
              tr,
              state_levels = state_levels_li,
              state_full_labels = state_full_labels_li,
              state_colors = state_colors_li,
              highlight_state = hs,
              fade_others = fade,
              hide_no_switch = hide_ns
            )
          }
        } else {  # Overall (non-lithium) focus
          if (color_by_first) {
            p <- create_sankey_first_state(
              tr,
              highlight_state = hs,
              fade_others = fade,
              hide_no_switch = hide_ns
            )
          } else {
            p <- create_sankey(
              tr,
              highlight_state = hs,
              fade_others = fade,
              hide_no_switch = hide_ns
            )
          }
        }
      }
      
      ggsave(file, p, dpi = 300, width = 13, height = 9, bg = "white")
    }
  )
  
  output$dl_heatmap <- downloadHandler(
    filename = function() sprintf("heatmap_%s.png", Sys.Date()),
    content  = function(file) {
      
      choice <- input$which_heatmap_download
      
      if (choice == "fs_overall") {
        # First → Second (overall)
        axis_order <- c("AD", "AP", "MS", "AD; AP", "AD; MS", "AP; MS",
                        "AD; AP; MS", "No Switch", "Int/Disc")
        p <- make_transition_heatmap(
          trajectories_limit, axis_order = axis_order,
          title = "Initial \u2192 Second-line (Heat map)"
        )
        
      } else if (choice == "fs_lithium") {
        # First → Second (lithium-focused)
        axis_order <- c(
          "Lithium", "Lithium; AD", "Lithium; AP",
          "Lithium; MS", "Lithium; 2+",
          "Non-Lithium", "Switch Non-Lithium", "No Switch", "Int/Disc"
        )
        p <- make_transition_heatmap(
          trajectories_limit_li, axis_order = axis_order,
          title = "Initial \u2192 Second-line (Lithium)"
        )
        
      } else if (choice == "st_overall") {
        # Second → Third (overall)
        p <- second_to_third
        
      } else {  # "st_lithium"
        # Second → Third (lithium-focused)
        p <- second_to_third_li
      }
      
      ggsave(file, p, dpi = 300, width = 10, height = 6, bg = "white")
    }
  )
  
}

shinyApp(ui, server)
