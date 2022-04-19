# Script Author: JaeJin Choi, OCT-2020

# required for actual work
require(ggplot2)
require(dplyr)
# require(plotrix) #two y-axis plot; twoord.plot
require(ggtree)
# require(phangorn) #contains upgma
# require(ape)
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

    #show a point of minimum Robinson-Foulds topological distance, in blue vertical lines
    if (3 %in% trend_indicators)
    {
        plot_obj <- plot_obj + geom_vline(
            xintercept = input_df$feature_length[input_df$rf_tree_dist==min(input_df$rf_tree_dist, na.rm=T)]
            , linetype = "dashed"
            , color = "black"
            , size = 0.7
        )
    }
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
    
    ## place vertical line and text
    plot_obj <- plot_obj + annotate(
        geom = "text"
        , label = label
        , x = at_x
        , y = at_y
        , colour = color
        , hjust=1.5
        , vjust=0.0
    )

    return (plot_obj)
}
