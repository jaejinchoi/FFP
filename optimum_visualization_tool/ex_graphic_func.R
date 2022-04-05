# Script Author: JaeJin Choi, OCT-2020

# required for actual work
require(ggplot2)
require(dplyr)
require(plotrix) #two y-axis plot; twoord.plot
require(ggtree)
require(phangorn) #contains upgma
require(ape)
require(tibble)

ggplot_draw_main <- function(plot_obj, input_df, x, x_title, y, y_title, aes_color, title, subtitle, add_smooth_line_flag=FALSE
                             , smooth_method, smooth_span, smooth_formula)
{
    # print(paste0(x, "\t", y))
    # print(input_df[x]) #correct
    # print(input_df$y) #NULL, is problem

    # https://www.rdocumentation.org/packages/ggplot2/versions/1.0.0/topics/aes_string
    
    plot_obj <- ggplot(data=input_df, aes_string(x=x, y=y, colour=aes_color)) + #aes not accepting string type variants, instead use aes_string
        theme_bw() + #white background with black edges
        
        ## how to place multiple columns of data?
        geom_point() +
        geom_line() +

        # theme: labels
        labs(x=x_title
             , y=y_title
             , title=title
             , subtitle=subtitle
        ) +

        # theme: orientation, style, size
        theme(
            axis.text.x=element_text(size=13, family='Arial')
            , axis.text.y=element_text(size=13, family='Arial')
            , axis.title.x=element_text(size=17, family='Arial', vjust=-2, hjust=0.5)
            , axis.title.y=element_text(size=17, family='Arial', vjust=2, hjust=0.5)
            , plot.title=element_text(size=20, vjust=4, hjust=0.5)
            , plot.subtitle=element_text(size=17, vjust=4, hjust=0.5)
            , plot.margin=unit(c(1,1,1,1), "cm")
            , plot.caption=element_text(size=13, hjust=0, vjust=-5)
        )

    ## https://stats.oarc.ucla.edu/r/faq/how-can-i-explore-different-smooths-in-ggplot2/
    if (add_smooth_line_flag==TRUE)
    {
        plot_obj <- plot_obj + stat_smooth(
          data=input_df
          , aes_string(x=x, y=y)
          , method=smooth_method #method="loess"
          , formula=smooth_formula
          , span=smooth_span #span=0.75 #default
          , linetype="dashed"
          # , se=F
          , na.rm=T
          )
    }

    return (plot_obj)

}

ggplot_add_indicators <- function(plot_obj, input_df, trend_indicators)
{
    # geom_vline(xintercept, linetype, color, size) #alone with embbeding does not work
    ## abline() does not work with ggplot
    if (1 %in% trend_indicators && sum(input_df$matrix_upperlimit)!=0)
    {
        plot_obj <- plot_obj + geom_vline(
            xintercept = input_df$feature_length[input_df$matrix_upperlimit!=0][1]
            # , linetype = "fine" #this causes error?!
            , color = "blue"
            , size = 0.7
        )

    }

    #show a point of maximum JSD matrix standard deviation (or deviation), in red vertical lines
    if (2 %in% trend_indicators)
    {
        plot_obj <- plot_obj + geom_vline(
            xintercept = input_df$feature_length[input_df$matrix_wout_upperlimit!=0][1]
            , linetype = "dashed"
            , color = "blue"
            , size = 0.7
        )
    }

    # #show a point of minimum Robinson-Foulds topological distance, in blue vertical lines
    # if (3 %in% trend_indicators)
    # {
    #     plot_obj <- plot_obj + geom_vline(
    #         xintercept = input_df$feature_length[input_df$rf_tree_dist==min(input_df$rf_tree_dist, na.rm=T)]
    #         , linetype = "dashed"
    #         , color = "black"
    #         , size = 0.7
    #     )
    # }
    # 
    # #show a point of minimum Subtree Pruning-and-Regrafting tree distance, in blue vertical lines
    # if (4 %in% trend_indicators)
    # {
    #     plot_obj <- plot_obj + geom_vline(
    #         xintercept = input_df$feature_length[input_df$spr_tree_dist==min(input_df$spr_tree_dist, na.rm=T)]
    #         , linetype = "dashed"
    #         , color = "orange"
    #         , size = 0.7
    #     )
    # }
    
    return (plot_obj)
}

## on plots, add line and text indication of interest
add_geom_indicate <- function(plot_obj, at_x, at_y, line_type, label, color)
{
    ## draw vertical line
    plot_obj <- plot_obj + geom_vline(
        xintercept = at_x
        , colour = color
        , linetype = line_type
        , size = 0.7
    )
    
    #a vertical line with controled height
    # http://www.sthda.com/english/wiki/ggplot2-add-straight-lines-to-a-plot-horizontal-vertical-and-regression-lines
    #https://ggplot2.tidyverse.org/reference/geom_segment.html
    
    # print(at_x[1])
    # print(at_y[1])
    # 
    # plot_obj <- plot_obj + geom_segment(aes(x = 3, y = 0, xend = 3, yend=1
    #                  # aes(x = at_x[1], y=0, xend=at_x[1], yend=at_y[1])
    #                  , colour = color
    #                  , linetype = line_type
    #                  , size = 0.7)
    # )
    
    
    ## place text
    plot_obj <- plot_obj + annotate(
        geom = "text"
        , label = label
        , x = at_x
        , y = at_y
        , colour = color
        , hjust=1.5
        , vjust=0.0
    )
    
    # plot_obj <- plot_obj + annotate(
    #     geom = "text"
    #     , label = "median max"
    #     , x = at_x
    #     , y = at_y
    #     , colour = "blue"
    #     , hjust=1
    #     , vjust=0
    # )
    
    return (plot_obj)
}

twoord_draw_main <- function(plot_obj, input_df, x, x_title, left_y, left_y_title, right_y, right_y_title, title, subtitle)
{
    # Error in stripchart.default: invalid plotting method; occurs when a function look for vectors and not a single dataframe
    plot_obj <- twoord.plot(
        data = input_df
        , lx = x
        , ly = left_y #ly tick scale not shown!!
        , rx = x
        , ry = right_y

        # , rylim=range(input_df$jsd_sd, na.rm=T)
        # , lylim=range(input_df$tree_dist, na.rm=T)

        , xlab = x_title
        , ylab = left_y_title
        , rylab = right_y_title
        , main = title
        , sub = subtitle #sub(subscription)
        , axislab.cex=1.0 #axis and tick size multiplier
        , lcol = "black" #black
        , rcol = "red" #red
    ) #+

    return (plot_obj)
}


tw_add_abline <- function(plot_obj, at_x, line_type, color, line_width) #using abline
{
    abline(
        v = at_x
        , col = color
        , lty = line_type #e.g., dashed, fine, dotted
        , lwd = line_width #2 #long fine line
    )
}

tw_add_indicators <- function(plot_obj, input_df, trend_indicators) #using abline
{
    ### display other_indicators (check box). Can add additional indicators of interests but require to update ui part
    #show a point of JSD-upperlimit, in gray vertical lines
    if (1 %in% trend_indicators & sum(input_df$matrix_upperlimit)!=0)
    {
        abline(
            v = input_df$feature_length[input_df$matrix_upperlimit!=0][1] #the first feature length showing JSD_upperlimit
            , col = "blue"
            , lty = "fine"
            , lwd = 2 #long fine line
        )
    }

    #show a point of maximum JSD matrix standard deviation (or deviation), in red vertical lines
    if (2 %in% trend_indicators)
    {
        abline(
            v = input_df$feature_length[input_df$matrix_sd==max(input_df$matrix_sd)] #the first feature length showing JSD_upperlimit
            , col = "red" #color blue==3
            , lty = "dashed" #line type
            , lwd = 1 #line thickness; line width
        ) #long dashed line

    }

    #show a point of minimum Robinson-Foulds topological distance, in blue vertical lines
    if (3 %in% trend_indicators)
    {
        abline(
            v = input_df$feature_length[input_df$rf_tree_dist==min(input_df$rf_tree_dist, na.rm=T)] #point(s) of minimal RF (e.g., RF==0)
            , col = "black" #color blue==3
            , lty = "dashed"
            , lwd = 1
        ) #long dashed line
    }

    return (plot_obj)

}


ggtree_draw_main <- function(input_dist, input_tree, tree_algorithm, tree_shape, branch_switch)
{
    #how to use ggtree to plot a phylo-tree in various shapes
    #https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html
    #https://stackoverflow.com/questions/61143008/add-another-layer-to-ggplot2-ggtree-based-on-user-input-rshiny

    ##2021-6-9. In case of caller_env missing (do the update instructed as below); still not working?
    #https://community.rstudio.com/t/creating-phylogenetic-trees-using-ggtree-caller-env-error-unable-to-run/105076/6
    ##latest solution report
    #https://stackoverflow.com/questions/67542581/ggtree-2-4-2-error-error-in-datamasknew-data-caller-env-argument-caller
    #temporary solution while R version below 4.1.0. Retain dplyr version to 1.0.5
    
    output_tree <- NULL
    tree_plot <- NULL
    
    ## list some other tree building algorithms
    if (!is.null(input_tree))
    {
        output_tree <- input_tree #just copy given input_tree, skip full calculation
        
    } else if (tree_algorithm=='bionj')
    {
        output_tree <- bionj(as.dist(input_dist))
        
    }else if (tree_algorithm=='nj')
    {
        output_tree <- nj(as.dist(input_dist))
        
    }else if (tree_algorithm=='upgma')
    {
        output_tree <- upgma(as.dist(input_dist))
        
    }

    ## parse and prepare tip labels
    labels <- c()

    if (isS4(output_tree)) {
        labels <- output_tree@phylo$tip.label
    } else {
        labels <- output_tree$tip.label
    }
    # print(labels) #tip labels

    #refer: https://github.com/tbradley1013/tree-subset-shiny/blob/master/server.R
    labels_df <- tibble( #parsing tips (regular expression for internal- and end-node parsing, in presence)
        label = labels,
        genus = str_extract(label, "[^;]+;[^;]+$") %>% str_replace(";[^;]+$", ""),
        species = str_extract(label, "[^;]+$")
    )  %>%
        mutate(
            species = if_else(is.na(genus), "", str_replace(species, "s__", "")),
            genus = if_else(is.na(genus), label, str_replace(genus, "g__", ""))
        )

    tree_dim.tip_size <- 4
    tree_dim.height <- length(labels) * 150 * tree_dim.tip_size

    # creating the plot and label
    # tree_plot <- output_tree %>%
    #     ggtree(layout=tree_shape
    #            , branch.length = branch_switch
    #            # , height = tree_dim.height #not used parameter
    #            # , yscale = tree_dim.height
    #     ) %<+%
    # 
    #     labels_df +
    #     geom_tiplab(aes(label = paste(genus, species))
    #                 , size = tree_dim.tip_size
    #                 ) + #least 9p is readable
    #     
    #     theme_tree2() + #theme for x-axis scaling
    # 
    #     # coord_cartesian(xlim = c(0, 40), ylim = c(0, 19)) +
    #     
    #     geom_text2(aes(label = branch.length), hjust = -.3)
    # 
    # #scale_color_manual(values = c(`1` = "red", `0` = "black"))
    # tree_plot <- lims(x = c(0, max(tree_plot$data$x) * length(labels))) #set x-axis limit, y-axis as no meaning in trees

    
    ##information on phylo(tree) object
    # http://phytools.org/eqg/Exercise_3.2/
    
    output_tree$edge.length_negative <- output_tree$edge.length>0
    
    # print(output_tree$edge.length_negative)
    # print(output_tree$edge.length)
    # print(output_tree$edge) #node pairs
    # print(output_tree$branch.length) #null
    
    ## plot examples: https://nbisweden.github.io/workshop-plotting-in-r/2109/lab_phylo.html
    tree_plot <- ggtree(tr=output_tree
                       , layout=tree_shape
                       , branch.length = branch_switch
                       , aes(color=branch.length>0)
                       ) +
        geom_tiplab(inherit.aes=F
                    , align=T
                    ) + #place tip(taxon) labels; all in black color
        
        theme_tree2() +

        # geom_text2(aes(label = sprintf(branch.length, fmt="%.2e") #edge.length not compatible?
        #                # , color=branch.length>0
        #                , x=branch #position to center
        #                )
        #                # , vjust = -1.0
        #                # , hjust = 1.0
        #                , show.legend=F
        #            ) +
        
        geom_label(aes(x=branch #at middle
                       , label = sprintf(branch.length, fmt="%.2e"))
                   ) +
    
        # scale_color_continuous(low="red", high="black") +
        scale_color_manual(values=c("red", "black")) +
        theme(legend.position="none") + #remove all legend
        
        scale_x_continuous(position = "top")  #place x-axis scale bar on top
    
    # print(output_tree$edge.length) #edge length
    # print(output_tree$edge.length<0)
    # print(output_tree$edge[output_tree$edge.length<0,]) #pickup edge pairs the lengths are negative
    
    # edge_df <- data.frame(output_tree$edge, output_tree$edge.length)
    # print(edge_df)
    
    # d = data.frame(node=1:59, color=sample(c('red', 'blue', 'green'), 59, replace=T))
    # ggtree(tr) %<+% d + aes(color=I(color))
    
    return(tree_plot)
    
}