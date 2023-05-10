#' Align sequences 
#'
#' Align a list of sequences to a reference
#'
#' @param reference Reference sequence or an array of reference sequences (array is only allowed if mode is 'individual')
#' @param sequences An array of sequences to be aligned to the reference(s)
#' @param mode 'individual' to align strings to reference(s) individually or 'common' to align to one common reference.
#'
#' @examples
#' align("ABCDEFGH", c("AAACEF", "ABCDEA"), "common")
#' 
#' @import grid
#' 
#' @export
align <- function(reference, sequences, mode = "individual") {
    .Call("align_c", reference, sequences, mode)
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
                letter, x = x + hstep / 2, y = y + vstep / 2, hjust = 0.5, vjust = 0.5,
                gp = gpar(col = "black", fontsize = 20, fontfamily = "mono")
            )
        )
    }

    grobs <- grobTree()

    draw_seqline <- function(grobs, y, splt) {
        index = 0
        for (char in splt) {
            fill <- "white"
            if (index + 1 < length(reference_split)) {
                refchar <- reference_split[[index + 1]]
                if (char != refchar) {
                    fill <- "gray90"
                }
            }
            grobs <- gTree(children = gList(grobs, letter_grob(char, index * hstep, y, fill)))
            index = index + 1
        }
        grobs
    }

    grobs <- draw_seqline(grobs, 1 - vstep, reference_split)
    {
        index = 2
        for (seqsplit in sequences_split) {
            grobs <- draw_seqline(grobs, 1 - vstep * index, seqsplit)
            index = index + 1
        }
    }

    grid.newpage()
    grid.draw(grobs)

    grobs
}
