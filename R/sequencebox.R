#' Align sequences
#'
#' Align a list of sequences to a reference
#'
#' @param reference Reference sequence or an array of reference sequences (array is only allowed if mode is 'individual')
#' @param sequences An array of sequences to be aligned to the reference(s)
#' @param mode 'individual' to align strings to reference(s) individually or 'common' to align to one common reference.
#' @param matrices Return alignment matrices
#'
#' @examples
#' align("ABCDEFGH", c("AAACEF", "ABCDEA"), "common")
#'
#' @import grid
#'
#' @export
align <- function(reference, sequences, mode = "individual", matrices = FALSE) {
    .Call("align_c", reference, sequences, mode, matrices)
}

#' Plot sequences
#'
#' Plot a list of sequences coloring each letter accoring to the reference.
#'
#' @param reference Reference sequence
#' @param sequences An array of sequences to be plotted
#'
#' @examples
#' plot_sequences("ABCDEFGH", c("AAACEF", "ABCDEA"))
#'
#' @import grid
#'
#' @export
plot_sequences <- function(reference, sequences) {
    if (length(reference) != 1) {
        stop(paste0("reference should be a single string, not a vector of length ", length(reference)))
    }

    reference_split <- strsplit(reference, "")[[1]]
    sequences_split <- strsplit(sequences, "")
    max_strlen <- length(reference_split)
    for (seqsplit in sequences_split) {
        max_strlen <- max(max_strlen, length(seqsplit))
    }

    hstep <- 1 / max_strlen
    vstep <- 1 / (length(sequences) + 1)

    letter_grob <- function(letter, x, y, rectFill = "white") {
        grobTree(
            rectGrob(
                x = x, y = y, hjust = 0, vjust = 0,
                width = hstep, height = vstep, gp = gpar(col = "white", fill = rectFill)
            ),
            textGrob(
                letter,
                x = x + hstep / 2, y = y + vstep / 2, hjust = 0.5, vjust = 0.5,
                gp = gpar(col = "black", fontsize = 20, fontfamily = "mono")
            )
        )
    }

    grobs <- grobTree()

    draw_seqline <- function(grobs, y, splt) {
        index <- 0
        for (char in splt) {
            fill <- "white"
            if (index + 1 < length(reference_split)) {
                refchar <- reference_split[[index + 1]]
                if (char != refchar) {
                    fill <- "gray90"
                }
            }
            grobs <- gTree(children = gList(grobs, letter_grob(char, index * hstep, y, fill)))
            index <- index + 1
        }
        grobs
    }

    grobs <- draw_seqline(grobs, 1 - vstep, reference_split)
    {
        index <- 2
        for (seqsplit in sequences_split) {
            grobs <- draw_seqline(grobs, 1 - vstep * index, seqsplit)
            index <- index + 1
        }
    }

    # TODO(sen) Probably hide this behind an option
    grid.newpage()
    grid.draw(grobs)

    grobs
}

#' Plot alignment matrix
#'
#' Plot a Needleman–Wunsch alignment matrix. Consists of 2 parts -
#' a score matrix and a direction matrix.
#'
#' @param score Matrix with scores
#' @param direction Matrix with directions
#' @param ref Reference sequence
#' @param seq Sequence that was aligned
#'
#' @examples
#' ref <- "ABCD"
#' seq <- "ACC"
#' align_result <- align(ref, seq, matrices = TRUE)
#' plot_alignment_matrix(align_result$scores[[1]], align_result$directions[[1]], ref, seq)
#'
#' @import grid
#'
#' @export
plot_alignment_matrix <- function(score, direction, ref, seq) {
    stopifnot(nrow(score) == nrow(direction))
    stopifnot(ncol(score) == ncol(direction))
    stopifnot(nrow(score) == nchar(seq) + 1)
    stopifnot(ncol(score) == nchar(ref) + 1)

    grobs <- grobTree()

    hstep <- 1 / (ncol(score) + 1)
    vstep <- 1 / (nrow(score) + 1)

    createEntryText <- function(text, x, y, col = "black") {
        textGrob(
            text,
            x = x, y = 1 - y, hjust = 0.5, vjust = 0.5,
            gp = gpar(col = col, fontsize = 20, fontfamily = "mono")
        )
    }

    {
        onpath_row <- nrow(score)
        onpath_col <- ncol(score)
        for (row in seq(nrow(score), 1)) {
            for (col in seq(ncol(score), 1)) {
                topleftx = col * hstep
                toplefty = row * vstep

                color <- "gray20"
                dir <- direction[row, col]
                if (row == onpath_row & col == onpath_col) {
                    color <- "blue"
                    if (dir == 0) {
                        onpath_row = onpath_row - 1
                        onpath_col = onpath_col - 1                            
                    } else if (dir == 1) {
                        onpath_row = onpath_row - 1
                    } else if (dir == 2) {
                        onpath_col = onpath_col - 1
                    }
                }

                scoregrob <- createEntryText(score[row, col], topleftx + hstep / 2, toplefty + vstep / 2, color)

                top <- 1 - toplefty
                left <- topleftx
                centerY = 1 - (toplefty + vstep / 2)
                centerX <- topleftx + hstep / 2
                dx <- hstep / 6
                dy <- - vstep / 6
                arr <- arrow(length = unit(0.02, "npc"))
                
                pars <- gpar(col = color)
                dirgrob <- switch(
                    dir + 1,
                    segmentsGrob(left + dx, top + dy, left - dx, top - dy, arrow = arr, gp = pars),
                    segmentsGrob(centerX, top + dy, centerX, top - dy, arrow = arr, gp = pars),
                    segmentsGrob(left + dx, centerY, left - dx, centerY, arrow = arr, gp = pars),
                )

                if (row > 1 | col > 1) {
                    grobs <- gTree(children = gList(grobs, scoregrob, dirgrob))
                } else {
                    grobs <- gTree(children = gList(grobs, scoregrob))
                }
            }
        }
    }

    {
        refsplit <- strsplit(ref, "")[[1]]
        index <- 2
        for (char in refsplit) {
            lettergrob <- createEntryText(char, index * hstep + hstep / 2, vstep / 2)
            grobs <- gTree(children = gList(grobs, lettergrob))
            index = index + 1
        }
    }

    {
        seqsplit <- strsplit(seq, "")[[1]]
        index <- 2
        for (char in seqsplit) {
            lettergrob <- createEntryText(char, hstep / 2, index * vstep + vstep / 2)
            grobs <- gTree(children = gList(grobs, lettergrob))
            index = index + 1
        }
    }

    # TODO(sen) Probably hide this behind an option
    grid.newpage()
    grid.draw(grobs)

    grobs
}

#' Generate random sequence
#'
#' @param src Source string to get characters out of
#' @param len Length of the random sequences
#' 
#' @examples
#' generate_random_sequence("ATGC", 10)
#'
#' @export
generate_random_sequence <- function(src, len) {
    .Call("generate_random_sequence_c", src, len)
}
