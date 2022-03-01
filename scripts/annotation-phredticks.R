#' Annotation: Phred Scaled Metrics tick marks
#'
#' This annotation adds log tick marks for phred scaled metrics, which are on a decible scale
#'
#' @param base the base of the log (default 10)
#' @param sides a string that controls which sides of the plot the log ticks appear on.
#'   It can be set to a string containing any of `"trbl"`, for top, right,
#'   bottom, and left.
#' @param outside logical that controls whether to move the log ticks outside
#' of the plot area. Default is off (`FALSE`). You will also need to use
#' `coord_cartesian(clip = "off")`. See examples.
#' @param short a [grid::unit()] object specifying the length of the
#'   short tick marks
#' @param mid a [grid::unit()] object specifying the length of the
#'   middle tick marks. In base 10, these are the "5" ticks.
#' @param long a [grid::unit()] object specifying the length of the
#'   long tick marks. In base 10, these are the "1" (or "10") ticks.
#' @param scaled is the data already log-scaled? This should be `TRUE`
#'   (default) when the data is already transformed with `log10()` or when
#'   using `scale_y_log10()`. It should be `FALSE` when using
#'   `coord_trans(y = "log10")`.
#' @param colour Colour of the tick marks.
#' @param size Thickness of tick marks, in mm.
#' @param linetype Linetype of tick marks (`solid`, `dashed`, etc.)
#' @param alpha The transparency of the tick marks.
#' @param color An alias for `colour`.
#' @param ... Other parameters passed on to the layer
#'
#' @export
#' @seealso [coord_trans()] for coordinate transformations.
#'
#' @examples
#' # TODO
#'
annotation_phredticks <- function(base = 10, sides = "bl", outside = FALSE, scaled = TRUE,
                                short = unit(0.1, "cm"), mid = unit(0.2, "cm"), long = unit(0.3, "cm"),
                                colour = "black", size = 0.5, linetype = 1, alpha = 1, color = NULL, ...)
{
    if (!is.null(color))
        colour <- color

    layer(
        data = ggplot2:::dummy_data(),
        mapping = NULL,
        stat = StatIdentity,
        geom = GeomPhredticks,
        position = PositionIdentity,
        show.legend = FALSE,
        inherit.aes = FALSE,
        params = list(
            base = base,
            sides = sides,
            outside = outside,
            scaled = scaled,
            short = short,
            mid = mid,
            long = long,
            colour = colour,
            size = size,
            linetype = linetype,
            alpha = alpha,
            ...
        )
    )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomPhredticks <- ggproto("GeomPhredticks", Geom,
                        extra_params = "",
                        handle_na = function(data, params) {
                            data
                        },
            
                        draw_panel = function(data, panel_params, coord, base = 10, sides = "bl",
                                              outside = FALSE, scaled = TRUE, short = unit(0.1, "cm"),
                                              mid = unit(0.2, "cm"), long = unit(0.3, "cm"))
                        {
                            require(grid)
                            ticks <- list()
                            flipped <- inherits(coord, "CoordFlip")
                            x_name <- if (flipped) "y" else "x"
                            y_name <- if (flipped) "x" else "y"

                            # Convert these units to numbers so that they can be put in data frames
                            short <- grid::convertUnit(short, "cm", valueOnly = TRUE)
                            mid   <- grid::convertUnit(mid,   "cm", valueOnly = TRUE)
                            long  <- grid::convertUnit(long,  "cm", valueOnly = TRUE)

                            if (grepl("[b|t]", sides)) {

                                # Get positions of x tick marks
                                xticks <- calc_phredticks(
                                    minphred = floor(panel_params$x.range[1]),
                                    maxphred = ceiling(panel_params$x.range[2]),
                                    start = 0,
                                    shortend = short,
                                    midend = mid,
                                    longend = long
                                )
                                # TODO implement for phred
                                #if (scaled)
                                #    xticks$value <- log(xticks$value, base)
                        
                                # Rename to 'x' for coordinates$transform
                                names(xticks)[names(xticks) == "value"] <- x_name
                                xticks <- coord$transform(xticks, panel_params)
                                xticks = xticks[xticks$x <= 1 & xticks$x >= 0,]

                                if (outside)
                                    xticks$end = -xticks$end

                                # Make the grobs
                                if (grepl("b", sides)) {
                                    ticks$x_b <- with(data, grid::segmentsGrob(
                                        x0 = unit(xticks$x, "native"), x1 = unit(xticks$x, "native"),
                                        y0 = unit(xticks$start, "cm"), y1 = unit(xticks$end, "cm"),
                                        gp = grid::gpar(col = alpha(colour, alpha), lty = linetype, lwd = size * .pt)
                                    ))
                                }
                                if (grepl("t", sides)) {
                                    ticks$x_t <- with(data, grid::segmentsGrob(
                                        x0 = unit(xticks$x, "native"), x1 = unit(xticks$x, "native"),
                                        y0 = unit(1, "npc") - unit(xticks$start, "cm"), y1 = unit(1, "npc") - unit(xticks$end, "cm"),
                                        gp = grid::gpar(col = alpha(colour, alpha), lty = linetype, lwd = size * .pt)
                                    ))
                                }
                            }

                        if (grepl("[l|r]", sides)) {
                            yticks <- calc_phredticks(
                                minphred = 0,#floor(panel_params$y.range[1]),
                                maxphred = 30,#ceiling(panel_params$y.range[2]),
                                start = 0,
                                shortend = short,
                                midend = mid,
                                longend = long
                            )

                            ## TODO  - To implement for phred
                            if (scaled)
                               yticks$value <- log(yticks$value, base)

                            names(yticks)[names(yticks) == "value"] <- y_name   # Rename to 'y' for coordinates$transform
                            yticks <- coord$transform(yticks, panel_params)
                            yticks = yticks[yticks$y <= 1 & yticks$y >= 0,]

                            if (outside)
                                yticks$end = -yticks$end

                            # Make the grobs
                                if (grepl("l", sides)) {
                                    ticks$y_l <- with(data, grid::segmentsGrob(
                                        y0 = unit(yticks$y, "native"), y1 = unit(yticks$y, "native"),
                                        x0 = unit(yticks$start, "cm"), x1 = unit(yticks$end, "cm"),
                                        gp = grid::gpar(col = alpha(colour, alpha), lty = linetype, lwd = size * .pt)
                                    ))
                                }
                                if (grepl("r", sides)) {
                                    ticks$y_r <- with(data, grid::segmentsGrob(
                                        y0 = unit(yticks$y, "native"), y1 = unit(yticks$y, "native"),
                                        x0 = unit(1, "npc") - unit(yticks$start, "cm"), x1 = unit(1, "npc") - unit(yticks$end, "cm"),
                                        gp = grid::gpar(col = alpha(colour, alpha), lty = linetype, lwd = size * .pt)
                                    ))
                                }
                            }

                            grid::gTree(children = do.call("gList", ticks))
                        },

                        default_aes = aes(colour = "black", size = 0.5, linetype = 1, alpha = 1)
)


# Calculate the position of phred scale tick marks
# Returns data frame with:
# - value: the position of the log tick on the data axis, for example 1, 2, ..., 9, 10, 20, ...
# - start: on the other axis, start position of the line (usually 0)
# - end: on the other axis, end position of the line (for example, .1, .2, or .3)
calc_phredticks <- function(
                          minphred = 0, maxphred = minphred + 1,
                          start = 0,
                          shortend = .1, midend = .2, longend = .3) {

    ## Phred tick df
    tick_df <-
        tibble(x_pct = c(0, seq(10, 90, 10),
                         seq(91, 99, 1),
                         seq(99.1, 99.9, 0.1),
                         seq(99.91, 99.99, 0.01),
                         seq(99.991, 99.999, 0.001)
                         ),
               ticknums = c(0, rep(1:9, times = 5
                              ))) %>%
        mutate(x = x_pct / 100,
               phred_x = transform_phred(x))

    ## Adding tick length info
    tick_df <- tick_df %>%
        mutate(
            start = start,
            tickend = case_when(
                x %in% c(0) ~ longend,
                ticknums == 9 ~ longend,
                ticknums == 5  ~ midend,
                TRUE ~ shortend
            )
        )

    ## Limiting to data range
    tickdf <- ggplot2:::new_data_frame(
        list(
            value = tick_df$phred_x,
            start = tick_df$start,
            end = tick_df$tickend
        ),
        n = nrow(tick_df)
    )
    return(tickdf)
}

transform_phred <- function(x) -10 * log(1 - x) / log(10)
