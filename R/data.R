#' Round-robbin game statistics from the Eastern NBA conference during the 2014-2015 season.
#'
#' A dataset of tensors containing pair-wise statistics between teams in the first half
#' and second half of the 2014-2015 season. This dataset only contains statistics for
#' the Eastern Conference of the National Basketball League (NBA).
#'
#' @format A list of tensors with the following elements:
#' \describe{
#'   \item{east_array}{A four-dimensional tensor. Element (i, j, k, l) contains
#'   the arc-sin transformed statistic k (k = free-throw proportion,
#'   two-point field goal proportion, or three-point field goal proportion)
#'   between team i and team j in either their first home (l = 1) or away (l = 2)
#'   game of the 2014-2015 season.}
#'   \item{east_array_n}{A four-dimensional tensor. Element (i, j, k, l) contains
#'   the number of attempts of k (k = free-throws,
#'   two-point field goals, or three-point field goals)
#'   between team i and team j in either their first home (l = 1) or away (l = 2)
#'   game of the 2014-2015 season.}
#'   \item{east_array_sub}{A four-dimensional tensor. Element (i, j, k, l) contains
#'   the number of made attempts of k (k = free-throws,
#'   two-point field goals, or three-point field goals)
#'   between team i and team j in either their first home (l = 1) or away (l = 2)
#'   game of the 2014-2015 season.}
#'   \item{east_true_as}{A four-dimensional tensor. Element (i, j, k, l) contains
#'   the arc-sin transformed statistic k (k = free-throw proportion,
#'   two-point field goal proportion, or three-point field goal proportion)
#'   between team i and team j in either their last home (l = 1) or away (l = 2)
#'   game of the 2014-2015 season.}
#'   \item{east_true_p}{A four-dimensional tensor. Element (i, j, k, l) contains
#'   statistic k (k = free-throw proportion,
#'   two-point field goal proportion, or three-point field goal proportion)
#'   between team i and team j in either their first home (l = 1) or away (l = 2)
#'   game of the 2014-2015 season.}
#' }
#' @source \url{http://www.basketball-reference.com/}
#'
#' @seealso \code{\link{west}}
#'
"east"


#' Round-robbin game statistics from the Western NBA conference during the 2014-2015 season.
#'
#' A dataset of tensors containing pair-wise statistics between teams in the first half
#' and second half of the 2014-2015 season. This dataset only contains statistics for
#' the Western Conference of the National Basketball League (NBA).
#'
#' @format A list of tensors with the following elements:
#' \describe{
#'   \item{west_array}{A four-dimensional tensor. Element (i, j, k, l) contains
#'   the arc-sin transformed statistic k (k = free-throw proportion,
#'   two-point field goal proportion, or three-point field goal proportion)
#'   between team i and team j in either their first home (l = 1) or away (l = 2)
#'   game of the 2014-2015 season.}
#'   \item{west_array_n}{A four-dimensional tensor. Element (i, j, k, l) contains
#'   the number of attempts of k (k = free-throws,
#'   two-point field goals, or three-point field goals)
#'   between team i and team j in either their first home (l = 1) or away (l = 2)
#'   game of the 2014-2015 season.}
#'   \item{west_array_sub}{A four-dimensional tensor. Element (i, j, k, l) contains
#'   the number of made attempts of k (k = free-throws,
#'   two-point field goals, or three-point field goals)
#'   between team i and team j in either their first home (l = 1) or away (l = 2)
#'   game of the 2014-2015 season.}
#'   \item{west_true_as}{A four-dimensional tensor. Element (i, j, k, l) contains
#'   the arc-sin transformed statistic k (k = free-throw proportion,
#'   two-point field goal proportion, or three-point field goal proportion)
#'   between team i and team j in either their last home (l = 1) or away (l = 2)
#'   game of the 2014-2015 season.}
#'   \item{west_true_p}{A four-dimensional tensor. Element (i, j, k, l) contains
#'   statistic k (k = free-throw proportion,
#'   two-point field goal proportion, or three-point field goal proportion)
#'   between team i and team j in either their first home (l = 1) or away (l = 2)
#'   game of the 2014-2015 season.}
#' }
#' @source \url{http://www.basketball-reference.com/}
#'
#' @seealso \code{\link{east}}
#'
"west"
