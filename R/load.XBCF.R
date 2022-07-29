load.longBet <- function(fileName) {
    json_str = readChar(fileName, file.info(fileName)$size)
    obj = .Call(`_longBet_json_to_r`, json_str)  # model$tree_pnt
    class(obj) = "longBet"
    return(obj)
}