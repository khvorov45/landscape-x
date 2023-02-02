align <- function(reference, sequences) {
    .Call("align_c", reference, sequences)
}