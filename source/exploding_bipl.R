function(data,
         use_percentile_for_plotting=TRUE,
         percentile_included=0.99,
         use_convex_hull=FALSE,
         fact_=1.5, 
         draw_more_ticks=1,
         round_ticks_to=2,
         ticks_in_original_units=TRUE,
         plot_the_dens=TRUE,
         move_lines_dist=1,
         name=NULL,
         plot_all_the_data=TRUE,
         tick_size=8,
         cluster_var=NULL, 
         counts_or_densities=c('densities','counts')[1],
         smoothing_method=c('smooth.spline',0)[1],
         smooth_par=NULL,
         draw_ellipse=FALSE,
         variable_name_size=14,
         var_names=TRUE,
         gray_scale=TRUE,
         EVchoice = c(1,2),
         plotqual = TRUE
)
{
  if(use_percentile_for_plotting==TRUE & draw_ellipse==TRUE){
    plot_percentile_ellipse=TRUE
    plot_ellipse=FALSE
  }
  if(use_percentile_for_plotting==FALSE & draw_ellipse==TRUE) { 
    plot_percentile_ellipse=FALSE
    plot_ellipse=TRUE}
  if(draw_ellipse==FALSE){
    plot_percentile_ellipse=FALSE
    plot_ellipse=FALSE
  }
  
  compute.bagplot=function (x, y, factor = 3, na.rm = FALSE, approx.limit = 300, 
                            dkmethod = 2, precision = 1, verbose = FALSE, debug.plots = "no",alph=0.5) {
    "bagplot, version 2012/12/05, peter wolf"
    win <- function(dx, dy) {
      atan2(y = dy, x = dx)
    }
    
    find.hdepths.tp <- function(tp, data, number.of.directions=181){ # 121130
      ## standardize dimensions ##
      xy <- as.matrix(data); tp <- as.matrix(rbind(tp)); n.tp <- dim(tp)[1]
      for( j in 1:2) {
        xy[,j] <- xy[,j] - (h <- min(xy[,j], na.rm=TRUE))
        tp[,j] <- tp[,j] -  h
        if( 0 < (h <- max(xy[,j], na.rm=TRUE))){
          xy[,j] <- xy[,j]/h; tp[,j] <- tp[,j]/h
        }
      }
      ##loop over directions##
      phi    <- c(seq(0,180,length=number.of.directions)[-1]*(2*pi/360))
      sinphi <- c(sin(phi),1); cosphi <- c(cos(phi),0)
      RM1 <- round(digits=6,rbind(cosphi,sinphi))
      hdtp <- rep(length(xy[,1]),length(tp[,1]))
      for( j in seq(along=sinphi)){ #print(j)  
        xyt <- xy %*% RM1[,j]; tpt <- (tp %*% RM1[,j])[]
        xyt <- xyt[!is.na(xyt)] #; tpt <- sort(tpt)
        hdtp <- pmin(hdtp,(rank( c(tpt,xyt), ties.method="min"))[1:n.tp]
                     -rank( tpt,ties.method="min")
                     ,rank(-c(tpt,xyt), ties.method="min")[1:n.tp]
                     -rank(-tpt,ties.method="min")                
        )
      }
      hdtp
    }
    
    out.of.polygon <- function(xy, pg) {
      xy <- matrix(xy, ncol = 2)
      if (nrow(pg) == 1) 
        return(xy[, 1] == pg[1] & xy[, 2] == pg[2])
      m <- nrow(xy)
      n <- nrow(pg)
      limit <- -abs(1e-10 * diff(range(pg)))
      pgn <- cbind(diff(c(pg[, 2], pg[1, 2])), -diff(c(pg[, 
                                                          1], pg[1, 1])))
      S <- colMeans(xy)
      dxy <- cbind(S[1] - pg[, 1], S[2] - pg[, 2])
      if (!all(limit < apply(dxy * pgn, 1, sum))) {
        pg <- pg[n:1, ]
        pgn <- -pgn[n:1, ]
      }
      in.pg <- rep(TRUE, m)
      for (j in 1:n) {
        dxy <- xy - matrix(pg[j, ], m, 2, byrow = TRUE)
        in.pg <- in.pg & limit < (dxy %*% pgn[j, ])
      }
      return(!in.pg)
    }
    cut.z.pg <- function(zx, zy, p1x, p1y, p2x, p2y) {
      a2 <- (p2y - p1y)/(p2x - p1x)
      a1 <- zy/zx
      sx <- (p1y - a2 * p1x)/(a1 - a2)
      sy <- a1 * sx
      sxy <- cbind(sx, sy)
      h <- any(is.nan(sxy)) || any(is.na(sxy)) || any(Inf == 
                                                        abs(sxy))
      if (h) {
        if (!exists("verbose")) 
          verbose <- FALSE
        if (verbose) 
          cat("special")
        h <- 0 == zx
        sx <- ifelse(h, zx, sx)
        sy <- ifelse(h, p1y - a2 * p1x, sy)
        a1 <- ifelse(abs(a1) == Inf, sign(a1) * 123456789 * 
                       1e+10, a1)
        a2 <- ifelse(abs(a2) == Inf, sign(a2) * 123456789 * 
                       1e+10, a2)
        h <- 0 == (a1 - a2) & sign(zx) == sign(p1x)
        sx <- ifelse(h, p1x, sx)
        sy <- ifelse(h, p1y, sy)
        h <- 0 == (a1 - a2) & sign(zx) != sign(p1x)
        sx <- ifelse(h, p2x, sx)
        sy <- ifelse(h, p2y, sy)
        h <- p1x == p2x & zx != p1x & p1x != 0
        sx <- ifelse(h, p1x, sx)
        sy <- ifelse(h, zy * p1x/zx, sy)
        h <- p1x == p2x & zx != p1x & p1x == 0
        sx <- ifelse(h, p1x, sx)
        sy <- ifelse(h, 0, sy)
        h <- p1x == p2x & zx == p1x & p1x != 0
        sx <- ifelse(h, zx, sx)
        sy <- ifelse(h, zy, sy)
        h <- p1x == p2x & zx == p1x & p1x == 0 & sign(zy) == 
          sign(p1y)
        sx <- ifelse(h, p1x, sx)
        sy <- ifelse(h, p1y, sy)
        h <- p1x == p2x & zx == p1x & p1x == 0 & sign(zy) != 
          sign(p1y)
        sx <- ifelse(h, p1x, sx)
        sy <- ifelse(h, p2y, sy)
        h <- zx == p1x & zy == p1y
        sx <- ifelse(h, p1x, sx)
        sy <- ifelse(h, p1y, sy)
        h <- zx == p2x & zy == p2y
        sx <- ifelse(h, p2x, sx)
        sy <- ifelse(h, p2y, sy)
        h <- zx == 0 & zy == 0
        sx <- ifelse(h, 0, sx)
        sy <- ifelse(h, 0, sy)
        sxy <- cbind(sx, sy)
      }
      if (!exists("debug.plots")) 
        debug.plots <- "no"
      if (debug.plots == "all") {
        segments(sxy[, 1], sxy[, 2], zx, zy, col = "red")
        segments(0, 0, sxy[, 1], sxy[, 2], col = "green", 
                 lty = 2)
        points(sxy, col = "red")
      }
      return(sxy)
    }
    find.cut.z.pg <- function(z, pg, center = c(0, 0), debug.plots = "no") {
      if (!is.matrix(z)) 
        z <- rbind(z)
      if (1 == nrow(pg)) 
        return(matrix(center, nrow(z), 2, TRUE))
      n.pg <- nrow(pg)
      n.z <- nrow(z)
      z <- cbind(z[, 1] - center[1], z[, 2] - center[2])
      pgo <- pg
      pg <- cbind(pg[, 1] - center[1], pg[, 2] - center[2])
      if (!exists("debug.plots")) 
        debug.plots <- "no"
      if (debug.plots == "all") {
        plot(rbind(z, pg, 0), bty = "n")
        points(z, pch = "p")
        lines(c(pg[, 1], pg[1, 1]), c(pg[, 2], pg[1, 2]))
      }
      apg <- win(pg[, 1], pg[, 2])
      apg[is.nan(apg)] <- 0
      a <- order(apg)
      apg <- apg[a]
      pg <- pg[a, ]
      az <- win(z[, 1], z[, 2])
      segm.no <- apply((outer(apg, az, "<")), 2, sum)
      segm.no <- ifelse(segm.no == 0, n.pg, segm.no)
      next.no <- 1 + (segm.no%%length(apg))
      cuts <- cut.z.pg(z[, 1], z[, 2], pg[segm.no, 1], pg[segm.no, 
                                                          2], pg[next.no, 1], pg[next.no, 2])
      cuts <- cbind(cuts[, 1] + center[1], cuts[, 2] + center[2])
      return(cuts)
    }
    hdepth.of.points <- function(tp) {
      n.tp <- nrow(tp)
      tphdepth <- rep(0, n.tp)
      dpi <- 2 * pi - 1e-06
      for (j in 1:n.tp) {
        dx <- tp[j, 1] - xy[, 1]
        dy <- tp[j, 2] - xy[, 2]
        a <- win(dx, dy) + pi
        h <- a < 10
        a <- a[h]
        ident <- sum(!h)
        init <- sum(a < pi)
        a.shift <- (a + pi)%%dpi
        minusplus <- c(rep(-1, length(a)), rep(1, length(a)))
        h <- cumsum(minusplus[order(c(a, a.shift))])
        tphdepth[j] <- init + min(h) + 1
      }
      tphdepth
    }
    expand.hull <- function(pg, k) {
      if (1 >= nrow(pg)) 
        return(pg)
      resolution <- floor(20 * precision)
      pg0 <- xy[hdepth == 1, ]
      pg0 <- pg0[chull(pg0[, 1], pg0[, 2]), ]
      end.points <- find.cut.z.pg(pg, pg0, center = center, 
                                  debug.plots = debug.plots)
      lam <- ((0:resolution)^1)/resolution^1
      pg.new <- pg
      for (i in 1:nrow(pg)) {
        tp <- cbind(pg[i, 1] + lam * (end.points[i, 1] - 
                                        pg[i, 1]), pg[i, 2] + lam * (end.points[i, 2] - 
                                                                       pg[i, 2]))
        hd.tp <- find.hdepths.tp(tp, xy)
        ind <- max(sum(hd.tp >= k), 1)
        if (ind < length(hd.tp)) {
          tp <- cbind(tp[ind, 1] + lam * (tp[ind + 1, 1] - 
                                            tp[ind, 1]), tp[ind, 2] + lam * (tp[ind + 1, 
                                                                                2] - tp[ind, 2]))
          hp.tp <- find.hdepths.tp(tp, xy)
          ind <- max(sum(hd.tp >= k), 1)
        }
        pg.new[i, ] <- tp[ind, ]
      }
      pg.new <- pg.new[chull(pg.new[, 1], pg.new[, 2]), ]
      pg.add <- 0.5 * (pg.new + rbind(pg.new[-1, ], pg.new[1, 
                                                           ]))
      end.points <- find.cut.z.pg(pg.add, pg0, center = center)
      for (i in 1:nrow(pg.add)) {
        tp <- cbind(pg.add[i, 1] + lam * (end.points[i, 1] - 
                                            pg.add[i, 1]), pg.add[i, 2] + lam * (end.points[i, 
                                                                                            2] - pg.add[i, 2]))
        hd.tp <- find.hdepths.tp(tp, xy)
        ind <- max(sum(hd.tp >= k), 1)
        if (ind < length(hd.tp)) {
          tp <- cbind(tp[ind, 1] + lam * (tp[ind + 1, 1] - 
                                            tp[ind, 1]), tp[ind, 2] + lam * (tp[ind + 1, 
                                                                                2] - tp[ind, 2]))
          hd.tp <- find.hdepths.tp(tp, xy)
          ind <- max(sum(hd.tp >= k), 1)
        }
        pg.add[i, ] <- tp[ind, ]
      }
      pg.new <- rbind(pg.new, pg.add)
      pg.new <- pg.new[chull(pg.new[, 1], pg.new[, 2]), ]
    }
    cut.p.sl.p.sl <- function(xy1, m1, xy2, m2) {
      sx <- (xy2[2] - m2 * xy2[1] - xy1[2] + m1 * xy1[1])/(m1 - 
                                                             m2)
      sy <- xy1[2] - m1 * xy1[1] + m1 * sx
      if (!is.nan(sy)) 
        return(c(sx, sy))
      if (abs(m1) == Inf) 
        return(c(xy1[1], xy2[2] + m2 * (xy1[1] - xy2[1])))
      if (abs(m2) == Inf) 
        return(c(xy2[1], xy1[2] + m1 * (xy2[1] - xy1[1])))
    }
    pos.to.pg <- function(z, pg, reverse = FALSE) {
      if (reverse) {
        int.no <- apply(outer(pg[, 1], z[, 1], ">="), 2, 
                        sum)
        zy.on.pg <- pg[int.no, 2] + pg[int.no, 3] * (z[, 
                                                       1] - pg[int.no, 1])
      }
      else {
        int.no <- apply(outer(pg[, 1], z[, 1], "<="), 2, 
                        sum)
        zy.on.pg <- pg[int.no, 2] + pg[int.no, 3] * (z[, 
                                                       1] - pg[int.no, 1])
      }
      result <- ifelse(z[, 2] < zy.on.pg, "lower", "higher")
      return(result)
      if (all(result == "lower")) {
        result <- ifelse(((z[, 2] - zy.on.pg)/max(z[, 2] - 
                                                    zy.on.pg) + 1e-10) < 0, "lower", "higher")
      }
      if (all(result == "higher")) {
        result <- ifelse(((z[, 2] - zy.on.pg)/max(z[, 2] - 
                                                    zy.on.pg) - 1e-10) < 0, "lower", "higher")
      }
      print(result)
      return(result)
    }
    find.polygon.center <- function(xy) {
      if (length(xy) == 2) 
        return(xy[1:2])
      if (nrow(xy) == 2) 
        return(colMeans(xy))
      n <- length(xy[, 1])
      mxy <- colMeans(xy)
      xy2 <- rbind(xy[-1, ], xy[1, ])
      xy3 <- cbind(rep(mxy[1], n), mxy[2])
      S <- (xy + xy2 + xy3)/3
      F2 <- abs((xy[, 1] - xy3[, 1]) * (xy2[, 2] - xy3[, 2]) - 
                  (xy[, 2] - xy3[, 2]) * (xy2[, 1] - xy3[, 1]))
      lambda <- F2/sum(F2)
      SP <- colSums(cbind(S[, 1] * lambda, S[, 2] * lambda))
      return(SP)
    }
    xydata <- if (missing(y)) 
      x
    else cbind(x, y)
    if (is.data.frame(xydata)) 
      xydata <- as.matrix(xydata)
    if (any(is.na(xydata))) {
      if (na.rm) {
        xydata <- xydata[!apply(is.na(xydata), 1, any), , 
                         drop = FALSE]
        print("Warning: NA elements have been removed!!")
      }
      else {
        xy.medians <- apply(xydata, 2, function(x) median(x, 
                                                          na.rm = TRUE))
        for (j in 1:ncol(xydata)) xydata[is.na(xydata[, j]), 
                                         j] <- xy.medians[j]
        print("Warning: NA elements have been exchanged by median values!!")
      }
    }
    if (length(xydata) < 4) {
      print("not enough data points")
      return()
    }
    if ((length(xydata)%%2) == 1) {
      print("number of values isn't even")
      return()
    }
    if (!is.matrix(xydata)) 
      xydata <- matrix(xydata, ncol = 2, byrow = TRUE)
    very.large.data.set <- nrow(xydata) > approx.limit
    if (very.large.data.set) {
      step <- (n <- nrow(xydata))/approx.limit
      ind <- round(seq(1, n, by = step))
      xy <- xydata[ind, ]
    }
    else xy <- xydata
    n <- nrow(xy)
    points.in.bag <- floor(n*alph)
    if (verbose) 
      cat("end of initialization")
    prdata <- prcomp(xydata)
    is.one.dim <- (0 == max(prdata[[1]])) || (min(prdata[[1]])/max(prdata[[1]])) < 
      1e-05
    if (is.one.dim) {
      if (verbose) 
        cat("data set one dimensional")
      center <- colMeans(xydata)
      res <- list(xy = xy, xydata = xydata, prdata = prdata, 
                  is.one.dim = is.one.dim, center = center)
      class(res) <- "bagplot"
      return(res)
    }
    if (verbose) 
      cat("data not linear")
    if (nrow(xydata) <= 4) {
      if (verbose) 
        cat("only three or four data points")
      center <- colMeans(xydata)
      res <- list(xy = xy, xydata = xydata, prdata = prdata, 
                  hdepths = rep(1, n), hdepth = rep(1, n), is.one.dim = is.one.dim, 
                  center = center, hull.center = NULL, hull.bag = NULL, 
                  hull.loop = NULL, pxy.bag = NULL, pxy.outer = xydata, 
                  pxy.outlier = NULL, exp.dk = xydata)
      class(res) <- "bagplot"
      return(res)
    }
    xym <- apply(xy, 2, mean)
    xysd <- apply(xy, 2, sd)
    xyxy <- cbind((xy[, 1] - xym[1])/xysd[1], (xy[, 2] - xym[2])/xysd[2])
    dx <- (outer(xy[, 1], xy[, 1], "-"))
    dy <- (outer(xy[, 2], xy[, 2], "-"))
    alpha <- atan2(y = dy, x = dx)
    diag(alpha) <- 1000
    for (j in 1:n) alpha[, j] <- sort(alpha[, j])
    alpha <- alpha[-n, ]
    m <- n - 1
    if (debug.plots == "all") {
      plot(xy, bty = "n")
      xdelta <- abs(diff(range(xy[, 1])))
      dx <- xdelta * 0.3
      for (j in 1:n) {
        p <- xy[j, ]
        dy <- dx * tan(alpha[, j])
        segments(p[1] - dx, p[2] - dy, p[1] + dx, p[2] + 
                   dy, col = j)
        text(p[1] - xdelta * 0.02, p[2], j, col = j)
      }
    }
    if (verbose) 
      print("end of computation of angles")
    hdepth <- rep(0, n)
    dpi <- 2 * pi - 1e-06
    mypi <- pi - 1e-06
    minusplus <- c(rep(-1, m), rep(1, m))
    if (FALSE) {
      for (j in 1:n) {
        a <- alpha[, j] + pi
        h <- a < 10
        a <- a[h]
        init <- sum(a < mypi)
        a.shift <- (a + pi)%%dpi
        minusplus <- c(rep(-1, length(a)), rep(1, length(a)))
        h <- cumsum(minusplus[order(c(a, a.shift))])
        hdepth[j] <- init + min(h) + 1
      }
    }
    find.hdepths <- function(xy, number.of.directions = 181) {
      xy <- as.matrix(xy)
      for (j in 1:2) {
        xy[, j] <- xy[, j] - min(xy[, j])
        if (0 < (h <- max(xy[, j]))) 
          xy[, j] <- xy[, j]/max(xy[, j])
      }
      phi <- c(seq(0, 180, length = number.of.directions)[-1] * 
                 (2 * pi/360))
      sinphi <- c(sin(phi), 1)
      cosphi <- c(cos(phi), 0)
      RM1 <- round(digits = 6, rbind(cosphi, sinphi))
      hd <- rep(h <- length(xy[, 1]), h)
      for (j in seq(along = sinphi)) {
        xyt <- xy %*% RM1[, j]
        hd <- pmin(hd, rank(xyt, ties.method = "min"), rank(-xyt, 
                                                            ties.method = "min"))
      }
      hd
    }
    hdepth <- find.hdepths(xy, 181 * precision)
    if (verbose) {
      print("end of computation of hdepth:")
      print(hdepth)
    }
    if (debug.plots == "all") {
      plot(xy, bty = "n")
      xdelta <- abs(diff(range(xy[, 1])))
      dx <- xdelta * 0.1
      for (j in 1:n) {
        a <- alpha[, j] + pi
        a <- a[a < 10]
        init <- sum(a < pi)
        a.shift <- (a + pi)%%dpi
        minusplus <- c(rep(-1, length(a)), rep(1, length(a)))
        h <- cumsum(minusplus[ao <- (order(c(a, a.shift)))])
        no <- which((init + min(h)) == (init + h))[1]
        p <- xy[j, ]
        dy <- dx * tan(alpha[, j])
        segments(p[1] - dx, p[2] - dy, p[1] + dx, p[2] + 
                   dy, col = j, lty = 3)
        dy <- dx * tan(c(sort(a), sort(a))[no])
        segments(p[1] - 5 * dx, p[2] - 5 * dy, p[1] + 5 * 
                   dx, p[2] + 5 * dy, col = "black")
        text(p[1] - xdelta * 0.02, p[2], hdepth[j], col = 1)
      }
    }
    hd.table <- table(sort(hdepth))
    d.k <- cbind(dk = rev(cumsum(rev(hd.table))), k = as.numeric(names(hd.table)))
    k.1 <- sum(points.in.bag < d.k[, 1])
    k <- d.k[k.1, 2] + 1
    if (verbose) {
      cat("numbers of members of dk:")
      print(hd.table)
      print(d.k)
    }
    if (verbose) {
      cat("end of computation of k, k=", k, "k.1:", k.1)
    }
    center <- apply(xy[which(hdepth == max(hdepth)), , drop = FALSE], 
                    2, mean)
    hull.center <- NULL
    if (3 < nrow(xy) && length(hd.table) > 0) {
      n.p <- floor(1.5 * c(32, 16, 8)[1 + (n > 50) + (n > 200)] * 
                     precision)
      h <- unique(sort(hdepth, decreasing = TRUE))
      limit.hdepth.to.check <- sort(h)[min(length(h), 3)]
      h <- cands <- xy[limit.hdepth.to.check <= hdepth, , drop = FALSE]
      cands <- cands[chull(cands[, 1], cands[, 2]), ]
      n.c <- nrow(cands)
      if (is.null(n.c)) 
        cands <- h
      xyextr <- rbind(apply(cands, 2, min), apply(cands, 2, 
                                                  max))
      if ((xyextr[2, 1] - xyextr[1, 1]) < 0.2 * (h <- diff(range(xy[, 
                                                                    1])))) {
        xyextr[1:2, 1] <- mean(xyextr[, 1]) + c(-0.1, 0.1) * 
          h
      }
      if ((xyextr[2, 2] - xyextr[1, 2]) < 0.2 * (h <- diff(range(xy[, 
                                                                    2])))) {
        xyextr[1:2, 2] <- mean(xyextr[, 2]) + c(-0.1, 0.1) * 
          h
      }
      if (verbose) {
        cat("xyextr: looking for maximal depth")
        print(xyextr)
      }
      h1 <- seq(xyextr[1, 1], xyextr[2, 1], length = n.p)
      h2 <- seq(xyextr[1, 2], xyextr[2, 2], length = n.p)
      tp <- cbind(as.vector(matrix(h1, n.p, n.p)), as.vector(matrix(h2, 
                                                                    n.p, n.p, TRUE)))
      tphdepth <- max(find.hdepths.tp(tp, xy))
      if (verbose) 
        cat("points(TP,pch=c(letters,LETTERS)[TPD+1])")
      if (verbose) {
        cat("depth of testpoints")
        print(summary(tphdepth))
      }
      tphdepth <- max(tphdepth, d.k[, 2])
      num <- floor(2 * c(417, 351, 171, 85, 67, 43)[sum(n > 
                                                          c(1, 50, 100, 150, 200, 250))] * precision)
      num.h <- floor(num/2)
      angles <- seq(0, pi, length = num.h)
      ang <- tan(pi/2 - angles)
      kkk <- tphdepth
      if (verbose) {
        cat("max-hdepth found:")
        print(kkk)
      }
      if (verbose) 
        cat("find polygon with max depth")
      ia <- 1
      a <- angles[ia]
      xyt <- xyxy %*% c(cos(a), -sin(a))
      xyto <- order(xyt)
      ind.k <- xyto[kkk]
      cutp <- c(xyxy[ind.k, 1], -10)
      dxy <- diff(range(xyxy))
      pg <- rbind(c(cutp[1], -dxy, Inf), c(cutp[1], dxy, NA))
      ind.kk <- xyto[n + 1 - kkk]
      cutpl <- c(xyxy[ind.kk, 1], 10)
      pgl <- rbind(c(cutpl[1], dxy, -Inf), c(cutpl[1], -dxy, 
                                             NA))
      if (debug.plots == "all") {
        plot(xyxy, type = "p", bty = "n")
        text(xy, , 1:n, col = "blue")
        hx <- xy[ind.k, c(1, 1)]
        hy <- xy[ind.k, c(2, 2)]
        segments(hx, hy, c(10, -10), hy + ang[ia] * (c(10, 
                                                       -10) - hx), lty = 2)
        text(hx + rnorm(1, , 0.1), hy + rnorm(1, , 0.1), 
             ia)
      }
      if (verbose) 
        cat("start of computation of the directions: ", "kkk=", 
            kkk)
      for (ia in seq(angles)[-1]) {
        a <- angles[ia]
        angtan <- ang[ia]
        xyt <- xyxy %*% c(cos(a), -sin(a))
        xyto <- order(xyt)
        ind.k <- xyto[kkk]
        ind.kk <- xyto[n + 1 - kkk]
        pnew <- xyxy[ind.k, ]
        pnewl <- xyxy[ind.kk, ]
        if (debug.plots == "all") 
          points(pnew[1], pnew[2], col = "red")
        if (abs(angtan) > 1e+10) {
          if (verbose) 
            cat("kkk", kkk, "x=c case")
          pg.no <- sum(pg[, 1] < pnew[1])
          if (0 < pg.no) {
            cutp <- c(pnew[1], pg[pg.no, 2] + pg[pg.no, 
                                                 3] * (pnew[1] - pg[pg.no, 1]))
            pg <- rbind(pg[1:pg.no, ], c(cutp, angtan), 
                        c(cutp[1] + dxy, cutp[2] + angtan * dxy, 
                          NA))
          }
          else {
            if (verbose) 
              cat("!!! case degenerated UPPER polygon: pg.no==0")
            pg <- rbind(pg[1, ], c(pg[2, 1:2], NA))
          }
          pg.nol <- sum(pgl[, 1] >= pnewl[1])
          if (0 < pg.nol) {
            cutpl <- c(pnewl[1], pgl[pg.nol, 2] + pgl[pg.nol, 
                                                      3] * (pnewl[1] - pgl[pg.nol, 1]))
            pgl <- rbind(pgl[1:pg.nol, ], c(cutpl, angtan), 
                         c(cutpl[1] - dxy, cutpl[2] - angtan * dxy, 
                           NA))
          }
          else {
            if (verbose) 
              cat("!!! case degenerated LOWER polygon: pgl.no==0")
            pgl <- rbind(pgl[1, ], c(pgl[2, 1:2], NA))
          }
        }
        else {
          pg.inter <- pg[, 2] - angtan * pg[, 1]
          pnew.inter <- pnew[2] - angtan * pnew[1]
          pg.no <- sum(pg.inter < pnew.inter)
          if (is.na(pg[pg.no, 3])) 
            pg[pg.no, 3] <- -Inf
          cutp <- cut.p.sl.p.sl(pnew, ang[ia], pg[pg.no, 
                                                  1:2], pg[pg.no, 3])
          pg <- rbind(pg[1:pg.no, ], c(cutp, angtan), c(cutp[1] + 
                                                          dxy, cutp[2] + angtan * dxy, NA))
          pg.interl <- pgl[, 2] - angtan * pgl[, 1]
          pnew.interl <- pnewl[2] - angtan * pnewl[1]
          pg.nol <- sum(pg.interl > pnew.interl)
          if (is.na(pgl[pg.nol, 3])) 
            pgl[pg.nol, 3] <- Inf
          cutpl <- cut.p.sl.p.sl(pnewl, angtan, pgl[pg.nol, 
                                                    1:2], pgl[pg.nol, 3])
          pgl <- rbind(pgl[1:pg.nol, ], c(cutpl, angtan), 
                       c(cutpl[1] - dxy, cutpl[2] - angtan * dxy, 
                         NA))
        }
        if (debug.plots == "all") {
          points(pnew[1], pnew[2], col = "red")
          hx <- xyxy[ind.k, c(1, 1)]
          hy <- xyxy[ind.k, c(2, 2)]
          segments(hx, hy, c(10, -10), hy + ang[ia] * (c(10, 
                                                         -10) - hx), lty = 2)
          points(cutpl[1], cutpl[2], col = "red")
          hx <- xyxy[ind.kk, c(1, 1)]
          hy <- xyxy[ind.kk, c(2, 2)]
          segments(hx, hy, c(10, -10), hy + ang[ia] * (c(10, 
                                                         -10) - hx), lty = 2)
        }
      }
      if (2 < nrow(pg) && 2 < nrow(pgl)) {
        limit <- 1e-10
        idx <- c(TRUE, (abs(diff(pg[, 1])) > limit) | (abs(diff(pg[, 
                                                                   2])) > limit))
        if (any(idx == FALSE)) {
          pg <- pg[idx, ]
          pg[, 3] <- c(diff(pg[, 2])/diff(pg[, 1]), NA)
        }
        idx <- c((abs(diff(pgl[, 1])) > limit) | (abs(diff(pgl[, 
                                                               2])) > limit), TRUE)
        if (any(idx == FALSE)) {
          pgl <- pgl[idx, ]
          pgl[, 3] <- c(diff(pgl[, 2])/diff(pgl[, 1]), 
                        NA)
        }
        pgl[, 2] <- pgl[, 2] - 1e-05
        pg <- pg[-nrow(pg), ][-1, , drop = FALSE]
        pgl <- pgl[-nrow(pgl), ][-1, , drop = FALSE]
        indl <- pos.to.pg(round(pgl, digits = 10), round(pg, 
                                                         digits = 10))
        indu <- pos.to.pg(round(pg, digits = 10), round(pgl, 
                                                        digits = 10), TRUE)
        sr <- sl <- NULL
        if (indu[(npg <- nrow(pg))] == "lower" & indl[1] == 
            "higher") {
          rnuml <- which(indl == "lower")[1] - 1
          rnumu <- npg + 1 - which(rev(indu == "higher"))[1]
          if (is.na(rnuml)) 
            rnuml <- sum(pg[rnumu, 1] < pgl[, 1])
          if (is.na(rnumu)) 
            rnumu <- sum(pg[, 1] < pgl[rnuml, 1])
          xyl <- pgl[rnuml, ]
          xyu <- pg[rnumu, ]
          sr <- cut.p.sl.p.sl(xyl[1:2], xyl[3], xyu[1:2], 
                              xyu[3])
        }
        if (indl[(npgl <- nrow(pgl))] == "higher" & indu[1] == 
            "lower") {
          lnuml <- npgl + 1 - which(rev(indl == "lower"))[1]
          lnumu <- which(indu == "higher")[1] - 1
          if (is.na(lnuml)) 
            lnuml <- sum(pg[lnumu, 1] < pgl[, 1])
          if (is.na(lnumu)) 
            lnumu <- sum(pg[, 1] < pgl[lnuml, 1])
          xyl <- pgl[lnuml, ]
          xyu <- pg[lnumu, ]
          sl <- cut.p.sl.p.sl(xyl[1:2], xyl[3], xyu[1:2], 
                              xyu[3])
        }
        pg <- rbind(pg[indu == "higher", 1:2, drop = FALSE], 
                    sr, pgl[indl == "lower", 1:2, drop = FALSE], 
                    sl)
        if (debug.plots == "all") 
          lines(rbind(pg, pg[1, ]), col = "red")
        if (!any(is.na(pg))) 
          pg <- pg[chull(pg[, 1], pg[, 2]), ]
      }
      else {
        if (2 < nrow(pgl)) {
          pg <- rbind(pg[2, 1:2], pgl[-c(1, length(pgl[, 
                                                       1])), 1:2])
        }
        else {
          pg <- rbind(pg[-c(1, length(pg[, 1])), 1:2], 
                      pgl[2, 1:2])
        }
      }
      if (verbose) 
        cat("END of computation of the directions")
      hull.center <- cbind(pg[, 1] * xysd[1] + xym[1], pg[, 
                                                          2] * xysd[2] + xym[2])
      if (!any(is.na(hull.center))) 
        center <- find.polygon.center(hull.center)
      else hull.center <- rbind(center)
      if (verbose) {
        cat("CENTER")
        print(center)
      }
      if (verbose) {
        cat("hull.center", hull.center)
        print(table(tphdepth))
      }
    }
    if (verbose) 
      cat("center depth:", find.hdepths.tp(rbind(center), xy) - 
            1)
    if (verbose) {
      print("end of computation of center")
      print(center)
    }
    if (dkmethod == 1) {
      xyi <- xy[hdepth >= k, , drop = FALSE]
      if (0 < length(xyi)) 
        pdk <- xyi[chull(xyi[, 1], xyi[, 2]), , drop = FALSE]
      if (k > 1) {
        xyo <- xy[hdepth >= (k - 1), , drop = FALSE]
        pdk.1 <- xyo[chull(xyo[, 1], xyo[, 2]), , drop = FALSE]
      }
      else pdk.1 <- pdk
      if (0 == length(xyi)) 
        pdk <- pdk.1
      if (verbose) 
        cat("hull computed: pdk, pdk.1:")
      if (verbose) {
        print(pdk)
        print(pdk.1)
      }
      if (debug.plots == "all") {
        plot(xy, bty = "n")
        h <- rbind(pdk, pdk[1, ])
        lines(h, col = "red", lty = 2)
        h <- rbind(pdk.1, pdk.1[1, ])
        lines(h, col = "blue", lty = 3)
        points(center[1], center[2], pch = 8, col = "red")
      }
      exp.dk <- expand.hull(pdk, k)
      exp.dk.1 <- expand.hull(exp.dk, k - 1)
    }
    else {
      num <- floor(2 * c(417, 351, 171, 85, 67, 43)[sum(n > 
                                                          c(1, 50, 100, 150, 200, 250))] * precision)
      num.h <- floor(num/2)
      angles <- seq(0, pi, length = num.h)
      ang <- tan(pi/2 - angles)
      kkk <- k
      if (verbose) 
        print("find polygon with depth something higher than that of the bag")
      if (kkk <= max(d.k[, 2])) {
        ia <- 1
        a <- angles[ia]
        xyt <- xyxy %*% c(cos(a), -sin(a))
        xyto <- order(xyt)
        ind.k <- xyto[kkk]
        cutp <- c(xyxy[ind.k, 1], -10)
        dxy <- diff(range(xyxy))
        pg <- rbind(c(cutp[1], -dxy, Inf), c(cutp[1], dxy, 
                                             NA))
        ind.kk <- xyto[n + 1 - kkk]
        cutpl <- c(xyxy[ind.kk, 1], 10)
        pgl <- rbind(c(cutpl[1], dxy, -Inf), c(cutpl[1], 
                                               -dxy, NA))
        if (debug.plots == "all") {
          plot(xyxy, type = "p", bty = "n")
          text(xy, , 1:n, col = "blue")
          hx <- xy[ind.k, c(1, 1)]
          hy <- xy[ind.k, c(2, 2)]
          segments(hx, hy, c(10, -10), hy + ang[ia] * (c(10, 
                                                         -10) - hx), lty = 2)
          text(hx + rnorm(1, , 0.1), hy + rnorm(1, , 0.1), 
               ia)
        }
        if (verbose) 
          cat("start of computation of the directions: ", 
              "kkk=", kkk)
        for (ia in seq(angles)[-1]) {
          a <- angles[ia]
          angtan <- ang[ia]
          xyt <- xyxy %*% c(cos(a), -sin(a))
          xyto <- order(xyt)
          ind.k <- xyto[kkk]
          ind.kk <- xyto[n + 1 - kkk]
          pnew <- xyxy[ind.k, ]
          pnewl <- xyxy[ind.kk, ]
          if (debug.plots == "all") 
            points(pnew[1], pnew[2], col = "red")
          if (abs(angtan) > 1e+10) {
            if (verbose) 
              cat("kkk", kkk, "x=c case")
            pg.no <- sum(pg[, 1] < pnew[1])
            if (0 < pg.no) {
              cutp <- c(pnew[1], pg[pg.no, 2] + pg[pg.no, 
                                                   3] * (pnew[1] - pg[pg.no, 1]))
              pg <- rbind(pg[1:pg.no, ], c(cutp, angtan), 
                          c(cutp[1] + dxy, cutp[2] + angtan * dxy, 
                            NA))
            }
            else {
              if (verbose) 
                cat("!!! case degenerated UPPER polygon: pg.no==0")
              pg <- rbind(pg[1, ], c(pg[2, 1:2], NA))
            }
            pg.nol <- sum(pgl[, 1] >= pnewl[1])
            if (0 < pg.nol) {
              cutpl <- c(pnewl[1], pgl[pg.nol, 2] + pgl[pg.nol, 
                                                        3] * (pnewl[1] - pgl[pg.nol, 1]))
              pgl <- rbind(pgl[1:pg.nol, ], c(cutpl, angtan), 
                           c(cutpl[1] - dxy, cutpl[2] - angtan * dxy, 
                             NA))
            }
            else {
              if (verbose) 
                cat("!!! case degenerated LOWER polygon: pgl.no==0")
              pgl <- rbind(pgl[1, ], c(pgl[2, 1:2], NA))
            }
          }
          else {
            pg.inter <- pg[, 2] - angtan * pg[, 1]
            pnew.inter <- pnew[2] - angtan * pnew[1]
            pg.no <- sum(pg.inter < pnew.inter)
            if (is.na(pg[pg.no, 3])) 
              pg[pg.no, 3] <- -Inf
            cutp <- cut.p.sl.p.sl(pnew, ang[ia], pg[pg.no, 
                                                    1:2], pg[pg.no, 3])
            pg <- rbind(pg[1:pg.no, ], c(cutp, angtan), 
                        c(cutp[1] + dxy, cutp[2] + angtan * dxy, 
                          NA))
            pg.interl <- pgl[, 2] - angtan * pgl[, 1]
            pnew.interl <- pnewl[2] - angtan * pnewl[1]
            pg.nol <- sum(pg.interl > pnew.interl)
            if (is.na(pgl[pg.nol, 3])) 
              pgl[pg.nol, 3] <- Inf
            cutpl <- cut.p.sl.p.sl(pnewl, angtan, pgl[pg.nol, 
                                                      1:2], pgl[pg.nol, 3])
            pgl <- rbind(pgl[1:pg.nol, ], c(cutpl, angtan), 
                         c(cutpl[1] - dxy, cutpl[2] - angtan * dxy, 
                           NA))
          }
          if (debug.plots == "all") {
            points(pnew[1], pnew[2], col = "red")
            hx <- xyxy[ind.k, c(1, 1)]
            hy <- xyxy[ind.k, c(2, 2)]
            segments(hx, hy, c(10, -10), hy + ang[ia] * 
                       (c(10, -10) - hx), lty = 2)
            points(cutpl[1], cutpl[2], col = "red")
            hx <- xyxy[ind.kk, c(1, 1)]
            hy <- xyxy[ind.kk, c(2, 2)]
            segments(hx, hy, c(10, -10), hy + ang[ia] * 
                       (c(10, -10) - hx), lty = 2)
          }
        }
        if (2 < nrow(pg) && 2 < nrow(pgl)) {
          limit <- 1e-10
          idx <- c(TRUE, (abs(diff(pg[, 1])) > limit) | 
                     (abs(diff(pg[, 2])) > limit))
          if (any(idx == FALSE)) {
            pg <- pg[idx, ]
            pg[, 3] <- c(diff(pg[, 2])/diff(pg[, 1]), NA)
          }
          idx <- c((abs(diff(pgl[, 1])) > limit) | (abs(diff(pgl[, 
                                                                 2])) > limit), TRUE)
          if (any(idx == FALSE)) {
            pgl <- pgl[idx, ]
            pgl[, 3] <- c(diff(pgl[, 2])/diff(pgl[, 1]), 
                          NA)
          }
          pgl[, 2] <- pgl[, 2] - 1e-05
          pg <- pg[-nrow(pg), ][-1, , drop = FALSE]
          pgl <- pgl[-nrow(pgl), ][-1, , drop = FALSE]
          indl <- pos.to.pg(round(pgl, digits = 10), round(pg, 
                                                           digits = 10))
          indu <- pos.to.pg(round(pg, digits = 10), round(pgl, 
                                                          digits = 10), TRUE)
          sr <- sl <- NULL
          if (indu[(npg <- nrow(pg))] == "lower" & indl[1] == 
              "higher") {
            rnuml <- which(indl == "lower")[1] - 1
            rnumu <- npg + 1 - which(rev(indu == "higher"))[1]
            if (is.na(rnuml)) 
              rnuml <- sum(pg[rnumu, 1] < pgl[, 1])
            if (is.na(rnumu)) 
              rnumu <- sum(pg[, 1] < pgl[rnuml, 1])
            xyl <- pgl[rnuml, ]
            xyu <- pg[rnumu, ]
            sr <- cut.p.sl.p.sl(xyl[1:2], xyl[3], xyu[1:2], 
                                xyu[3])
          }
          if (indl[(npgl <- nrow(pgl))] == "higher" & indu[1] == 
              "lower") {
            lnuml <- npgl + 1 - which(rev(indl == "lower"))[1]
            lnumu <- which(indu == "higher")[1] - 1
            if (is.na(lnuml)) 
              lnuml <- sum(pg[lnumu, 1] < pgl[, 1])
            if (is.na(lnumu)) 
              lnumu <- sum(pg[, 1] < pgl[lnuml, 1])
            xyl <- pgl[lnuml, ]
            xyu <- pg[lnumu, ]
            sl <- cut.p.sl.p.sl(xyl[1:2], xyl[3], xyu[1:2], 
                                xyu[3])
          }
          pg <- rbind(pg[indu == "higher", 1:2, drop = FALSE], 
                      sr, pgl[indl == "lower", 1:2, drop = FALSE], 
                      sl)
          if (debug.plots == "all") 
            lines(rbind(pg, pg[1, ]), col = "red")
          if (!any(is.na(pg))) 
            pg <- pg[chull(pg[, 1], pg[, 2]), ]
        }
        else {
          if (2 < nrow(pgl)) {
            pg <- rbind(pg[2, 1:2], pgl[-c(1, length(pgl[, 
                                                         1])), 1:2])
          }
          else {
            pg <- rbind(pg[-c(1, length(pg[, 1])), 1:2], 
                        pgl[2, 1:2])
          }
        }
        if (verbose) 
          cat("END of computation of the directions")
        exp.dk <- cbind(pg[, 1] * xysd[1] + xym[1], pg[, 
                                                       2] * xysd[2] + xym[2])
      }
      else {
        exp.dk <- NULL
      }
      if (1 < kkk) 
        kkk <- kkk - 1
      if (verbose) 
        print("find polygon with depth a little bit lower than that of the bag")
      ia <- 1
      a <- angles[ia]
      xyt <- xyxy %*% c(cos(a), -sin(a))
      xyto <- order(xyt)
      ind.k <- xyto[kkk]
      cutp <- c(xyxy[ind.k, 1], -10)
      dxy <- diff(range(xyxy))
      pg <- rbind(c(cutp[1], -dxy, Inf), c(cutp[1], dxy, NA))
      ind.kk <- xyto[n + 1 - kkk]
      cutpl <- c(xyxy[ind.kk, 1], 10)
      pgl <- rbind(c(cutpl[1], dxy, -Inf), c(cutpl[1], -dxy, 
                                             NA))
      if (debug.plots == "all") {
        plot(xyxy, type = "p", bty = "n")
        text(xy, , 1:n, col = "blue")
        hx <- xy[ind.k, c(1, 1)]
        hy <- xy[ind.k, c(2, 2)]
        segments(hx, hy, c(10, -10), hy + ang[ia] * (c(10, 
                                                       -10) - hx), lty = 2)
        text(hx + rnorm(1, , 0.1), hy + rnorm(1, , 0.1), 
             ia)
      }
      if (verbose) 
        cat("start of computation of the directions: ", "kkk=", 
            kkk)
      for (ia in seq(angles)[-1]) {
        a <- angles[ia]
        angtan <- ang[ia]
        xyt <- xyxy %*% c(cos(a), -sin(a))
        xyto <- order(xyt)
        ind.k <- xyto[kkk]
        ind.kk <- xyto[n + 1 - kkk]
        pnew <- xyxy[ind.k, ]
        pnewl <- xyxy[ind.kk, ]
        if (debug.plots == "all") 
          points(pnew[1], pnew[2], col = "red")
        if (abs(angtan) > 1e+10) {
          if (verbose) 
            cat("kkk", kkk, "x=c case")
          pg.no <- sum(pg[, 1] < pnew[1])
          if (0 < pg.no) {
            cutp <- c(pnew[1], pg[pg.no, 2] + pg[pg.no, 
                                                 3] * (pnew[1] - pg[pg.no, 1]))
            pg <- rbind(pg[1:pg.no, ], c(cutp, angtan), 
                        c(cutp[1] + dxy, cutp[2] + angtan * dxy, 
                          NA))
          }
          else {
            if (verbose) 
              cat("!!! case degenerated UPPER polygon: pg.no==0")
            pg <- rbind(pg[1, ], c(pg[2, 1:2], NA))
          }
          pg.nol <- sum(pgl[, 1] >= pnewl[1])
          if (0 < pg.nol) {
            cutpl <- c(pnewl[1], pgl[pg.nol, 2] + pgl[pg.nol, 
                                                      3] * (pnewl[1] - pgl[pg.nol, 1]))
            pgl <- rbind(pgl[1:pg.nol, ], c(cutpl, angtan), 
                         c(cutpl[1] - dxy, cutpl[2] - angtan * dxy, 
                           NA))
          }
          else {
            if (verbose) 
              cat("!!! case degenerated LOWER polygon: pgl.no==0")
            pgl <- rbind(pgl[1, ], c(pgl[2, 1:2], NA))
          }
        }
        else {
          pg.inter <- pg[, 2] - angtan * pg[, 1]
          pnew.inter <- pnew[2] - angtan * pnew[1]
          pg.no <- sum(pg.inter < pnew.inter)
          if (is.na(pg[pg.no, 3])) 
            pg[pg.no, 3] <- -Inf
          cutp <- cut.p.sl.p.sl(pnew, ang[ia], pg[pg.no, 
                                                  1:2], pg[pg.no, 3])
          pg <- rbind(pg[1:pg.no, ], c(cutp, angtan), c(cutp[1] + 
                                                          dxy, cutp[2] + angtan * dxy, NA))
          pg.interl <- pgl[, 2] - angtan * pgl[, 1]
          pnew.interl <- pnewl[2] - angtan * pnewl[1]
          pg.nol <- sum(pg.interl > pnew.interl)
          if (is.na(pgl[pg.nol, 3])) 
            pgl[pg.nol, 3] <- Inf
          cutpl <- cut.p.sl.p.sl(pnewl, angtan, pgl[pg.nol, 
                                                    1:2], pgl[pg.nol, 3])
          pgl <- rbind(pgl[1:pg.nol, ], c(cutpl, angtan), 
                       c(cutpl[1] - dxy, cutpl[2] - angtan * dxy, 
                         NA))
        }
        if (debug.plots == "all") {
          points(pnew[1], pnew[2], col = "red")
          hx <- xyxy[ind.k, c(1, 1)]
          hy <- xyxy[ind.k, c(2, 2)]
          segments(hx, hy, c(10, -10), hy + ang[ia] * (c(10, 
                                                         -10) - hx), lty = 2)
          points(cutpl[1], cutpl[2], col = "red")
          hx <- xyxy[ind.kk, c(1, 1)]
          hy <- xyxy[ind.kk, c(2, 2)]
          segments(hx, hy, c(10, -10), hy + ang[ia] * (c(10, 
                                                         -10) - hx), lty = 2)
        }
      }
      if (2 < nrow(pg) && 2 < nrow(pgl)) {
        limit <- 1e-10
        idx <- c(TRUE, (abs(diff(pg[, 1])) > limit) | (abs(diff(pg[, 
                                                                   2])) > limit))
        if (any(idx == FALSE)) {
          pg <- pg[idx, ]
          pg[, 3] <- c(diff(pg[, 2])/diff(pg[, 1]), NA)
        }
        idx <- c((abs(diff(pgl[, 1])) > limit) | (abs(diff(pgl[, 
                                                               2])) > limit), TRUE)
        if (any(idx == FALSE)) {
          pgl <- pgl[idx, ]
          pgl[, 3] <- c(diff(pgl[, 2])/diff(pgl[, 1]), 
                        NA)
        }
        pgl[, 2] <- pgl[, 2] - 1e-05
        pg <- pg[-nrow(pg), ][-1, , drop = FALSE]
        pgl <- pgl[-nrow(pgl), ][-1, , drop = FALSE]
        indl <- pos.to.pg(round(pgl, digits = 10), round(pg, 
                                                         digits = 10))
        indu <- pos.to.pg(round(pg, digits = 10), round(pgl, 
                                                        digits = 10), TRUE)
        sr <- sl <- NULL
        if (indu[(npg <- nrow(pg))] == "lower" & indl[1] == 
            "higher") {
          rnuml <- which(indl == "lower")[1] - 1
          rnumu <- npg + 1 - which(rev(indu == "higher"))[1]
          if (is.na(rnuml)) 
            rnuml <- sum(pg[rnumu, 1] < pgl[, 1])
          if (is.na(rnumu)) 
            rnumu <- sum(pg[, 1] < pgl[rnuml, 1])
          xyl <- pgl[rnuml, ]
          xyu <- pg[rnumu, ]
          sr <- cut.p.sl.p.sl(xyl[1:2], xyl[3], xyu[1:2], 
                              xyu[3])
        }
        if (indl[(npgl <- nrow(pgl))] == "higher" & indu[1] == 
            "lower") {
          lnuml <- npgl + 1 - which(rev(indl == "lower"))[1]
          lnumu <- which(indu == "higher")[1] - 1
          if (is.na(lnuml)) 
            lnuml <- sum(pg[lnumu, 1] < pgl[, 1])
          if (is.na(lnumu)) 
            lnumu <- sum(pg[, 1] < pgl[lnuml, 1])
          xyl <- pgl[lnuml, ]
          xyu <- pg[lnumu, ]
          sl <- cut.p.sl.p.sl(xyl[1:2], xyl[3], xyu[1:2], 
                              xyu[3])
        }
        pg <- rbind(pg[indu == "higher", 1:2, drop = FALSE], 
                    sr, pgl[indl == "lower", 1:2, drop = FALSE], 
                    sl)
        if (debug.plots == "all") 
          lines(rbind(pg, pg[1, ]), col = "red")
        if (!any(is.na(pg))) 
          pg <- pg[chull(pg[, 1], pg[, 2]), ]
      }
      else {
        if (2 < nrow(pgl)) {
          pg <- rbind(pg[2, 1:2], pgl[-c(1, length(pgl[, 
                                                       1])), 1:2])
        }
        else {
          pg <- rbind(pg[-c(1, length(pg[, 1])), 1:2], 
                      pgl[2, 1:2])
        }
      }
      if (verbose) 
        cat("END of computation of the directions")
      exp.dk.1 <- cbind(pg[, 1] * xysd[1] + xym[1], pg[, 2] * 
                          xysd[2] + xym[2])
      if (is.null(exp.dk)) 
        exp.dk <- exp.dk.1
      if (verbose) 
        print("End of find hulls, method two")
    }
    if (nrow(d.k) == k.1 || nrow(d.k) == 1) 
      lambda <- 0
    else {
      ind <- sum(d.k[, 2] <= k.1)
      ind <- k.1
      ndk.1 <- d.k[ind, 1]
      ndk <- d.k[ind + 1, 1]
      lambda <- (n*alph - ndk)/(ndk.1 - ndk)
    }
    if (verbose) 
      cat("lambda", lambda)
    cut.on.pdk.1 <- find.cut.z.pg(exp.dk, exp.dk.1, center = center)
    cut.on.pdk <- find.cut.z.pg(exp.dk.1, exp.dk, center = center)
    h1 <- (1 - lambda) * exp.dk + lambda * cut.on.pdk.1
    h2 <- (1 - lambda) * cut.on.pdk + lambda * exp.dk.1
    h <- rbind(h1, h2)
    h <- h[!is.nan(h[, 1]) & !is.nan(h[, 2]), ]
    hull.bag <- h[chull(h[, 1], h[, 2]), ]
    if (verbose) 
      cat("bag completed:")
    if (debug.plots == "all") {
      lines(hull.bag, col = "red")
    }
    hull.loop <- cbind(hull.bag[, 1] - center[1], hull.bag[, 
                                                           2] - center[2])
    hull.loop <- factor * hull.loop
    hull.loop <- cbind(hull.loop[, 1] + center[1], hull.loop[, 
                                                             2] + center[2])
    hull.fullloop <-  hull.loop 
    if (verbose) 
      cat("loop computed")
    if (!very.large.data.set) {
      pxy.bag <- xydata[hdepth >= k, , drop = FALSE]
      pkt.cand <- xydata[hdepth == (k - 1), , drop = FALSE]
      pkt.not.bag <- xydata[hdepth < (k - 1), , drop = FALSE]
      if (0 < length(pkt.cand) && 0 < length(hull.bag)) {
        outside <- out.of.polygon(pkt.cand, hull.bag)
        if (sum(!outside) > 0) 
          pxy.bag <- rbind(pxy.bag, pkt.cand[!outside, 
                                             ])
        if (sum(outside) > 0) 
          pkt.not.bag <- rbind(pkt.not.bag, pkt.cand[outside, 
                                                     ])
      }
    }
    else {
      extr <- out.of.polygon(xydata, hull.bag)
      pxy.bag <- xydata[!extr, ]
      pkt.not.bag <- xydata[extr, , drop = FALSE]
    }
    if (length(pkt.not.bag) > 0) {
      extr <- out.of.polygon(pkt.not.bag, hull.loop)
      pxy.outlier <- pkt.not.bag[extr, , drop = FALSE]
      if (0 == length(pxy.outlier)) 
        pxy.outlier <- NULL
      pxy.outer <- pkt.not.bag[!extr, , drop = FALSE]
    }
    else {
      pxy.outer <- pxy.outlier <- NULL
    }
    if (verbose) 
      cat("points of bag, outer points and outlier identified")
    hull.loop <- rbind(pxy.outer, hull.bag)
    hull.loop <- hull.loop[chull(hull.loop[, 1], hull.loop[, 
                                                           2]), ]
    if (verbose) 
      cat("end of computation of loop")
    res <- list(center = center, hull.center = hull.center, hull.bag = hull.bag, 
                hull.loop = hull.loop, pxy.bag = pxy.bag, pxy.outer = if (length(pxy.outer) > 
                                                                          0) pxy.outer else NULL, pxy.outlier = if (length(pxy.outlier) > 
                                                                                                                    0) pxy.outlier else NULL, hdepths = hdepth, is.one.dim = is.one.dim, 
                prdata = prdata, xy = xy, xydata = xydata, hull.fullloop = hull.fullloop)
    if (verbose) 
      res <- c(res, list(exp.dk = exp.dk, exp.dk.1 = exp.dk.1, 
                         hdepth = hdepth))
    class(res) <- "bagplot"
    return(res)
  }
  
  
  p1.unscaled=data
  
  if(plot_the_dens==TRUE & nrow(p1.unscaled)<30){
    print('Too few datapoints to construct density plots. Set plot_the_dens=FALSE')
  }
  if(is.null(cluster_var)==FALSE &any(table(cluster_var)<30)) {
    print('Too few datapoints per cluster to construct proper densities. Set plot_the dens=FALSE or cluster_var=NULL')
  }
  
  m_p1=apply(p1.unscaled,2,mean)
  s_p1=apply(p1.unscaled,2,std)
  p1.scaled=scale(p1.unscaled,center = TRUE,scale=TRUE)
  
  if(is.null(colnames(p1.scaled))==TRUE){
    colnames(p1.scaled)=paste0('Var',1:length(p1.scaled[1,]))
    colnames(p1.unscaled)=paste0('Var',1:length(p1.scaled[1,]))
  }
  
  #PCA
  PCA_biplot.inputs=function(dat,r=2){#r=2
    n=nrow(dat)
    p=ncol(dat)
    scaled.dat=scale(dat,center=TRUE,scale=TRUE)
    p_comp=svd(scaled.dat)
    SIG=diag(p_comp$d)
    U_r=p_comp$u[,EVchoice]
    SIG_r=SIG[EVchoice,EVchoice]
    V_r=p_comp$v[,EVchoice]
    points=U_r%*%SIG_r
    axes=V_r
    #coordinate of a marker on the kth biplot axes of 1 unit of the kth variable
    e_k=diag(1,nrow=p)
    c_proj=matrix(0,nrow=r,ncol=p)
    for(k in 1:p){
      t1=(t(V_r)%*%e_k[,k])
      t2=(t(e_k[,k])%*%V_r%*%t(V_r)%*%e_k[,k])#coordinate for predicted value of X
      t3=t1/as.numeric(t2)
      c_proj[,k]=t3
    }
    #quality of biplot
    
    svd.out <- p_comp
    V.mat <- svd.out$v
    U.mat <- svd.out$u
    Sigma.mat <- diag(svd.out$d)
    Vr.before.rotate <- svd.out$v[, EVchoice]
    eigval <- svd.out$d^2
    lambda.mat <- diag(eigval)
    eigval.r <- eigval[EVchoice]
    lambda.r.mat <- diag(eigval.r)
    fit.predictivity.mat <- diag(diag(Vr.before.rotate %*% 
                                        lambda.r.mat %*% t(Vr.before.rotate))) %*% solve(diag(diag(V.mat %*% 
                                                                                                     lambda.mat %*% t(V.mat))))
    fit.predictivity <- round(diag(fit.predictivity.mat),digits = 3)
    names(fit.predictivity) <- dimnames(dat)[[2]]
    fit.quality <- paste0("Quality of display = ", round(((eigval[EVchoice[1]] + 
                                                             eigval[EVchoice[2]])/sum(eigval)) * 100, digits = 2), "%")
    
    
    
    #quality=sum(p_comp$d[EVchoice])/sum(p_comp$d)
    
    ret=list('points'=points,'axes'=axes,'len_1_unit'=c_proj,quality=fit.quality,axqual = fit.predictivity)
    return(ret)
  }
  
  
  PCA=PCA_biplot.inputs(p1.unscaled)
  rownames(PCA$axes)=paste0('Var',1:length(PCA$axes[,1]))
  
  p11=PCA$points #all
  bag_=compute.bagplot(p11[,1],p11[,2],alph=percentile_included)$hull.bag
  p1=cbind(as.numeric(bag_[,1]),as.numeric(bag_[,2]))#over the bagplot
  axes=PCA$axes
  rownames(axes)=paste0('Var',1:length(PCA$axes[,1]))
  nr.of.axes=length(PCA$axes[,1])
  m=rep(0,(nr.of.axes))
  for(i in 1:length(m)){
    m[i]=(axes[i,2])/(axes[i,1]-0)
  }
  order_m_start=order(m,decreasing =FALSE)
  
  get_semi_axes_ellipse=function(t,e,exy.){
    qq=matrix(0,nrow=length(e),ncol=length(e))
    for(ww in 1:length(e)){
      for(vv in 1:length(e)){
        qq[ww,vv]=sqrt((t[ww]-t[vv])^2+(e[ww]-e[vv])^2)
      }
    }
    aa=max(qq)/2
  
    ind=which(qq==max(qq),arr.ind = TRUE) 
    delta=(e[ind[1]]-e[ind[2]])/(t[ind[1]]-t[ind[2]])
    mid=cbind((t[ind[1]]+t[ind[2]])/2,(e[ind[1]]+e[ind[2]])/2)

    A=atan(delta)
    if(abs(A*180/pi-90)<0.001 | is.na(A) |abs(A*180/pi-270)<0.001 ){ 
      A=pi/2
      a_ends=unique(cbind(t[ind[,1]],e[ind[,2]]))
      mid=cbind((a_ends[1,1]+a_ends[2,1])/2,(a_ends[1,2]+a_ends[2,2])/2)
    }
    eigen_val=eigen(exy.$cov)$values
    area=sqrt(eigen_val[1]*eigen_val[2])*pi*2
    b=area/(aa*pi)
    ret=list(aa,b,mid,A)
    return(ret)
  }
  
  #ellipse
  exy.=ellipsoidhull(p1)
  exy.full=ellipsoidhull(p11) #
  p_exy.=predict.ellipsoid(exy.,d2=d2.perc)
  p_exy.full=predict.ellipsoid(exy.full)
  t=p_exy.[,1]
  e=p_exy.[,2]
  t.full=p_exy.full[,1]
  e.full=p_exy.full[,2]
  semi_axes.full=get_semi_axes_ellipse(t.full,e.full,exy.full)
  semi_axes=get_semi_axes_ellipse(t,e,exy.) 
  a=semi_axes[[1]]
  b=semi_axes[[2]]
  mid=semi_axes[[3]]
  A=semi_axes[[4]]
  a.full=semi_axes.full[[1]]
  b.full=semi_axes.full[[2]]
  mid.full=semi_axes.full[[3]]
  A.full=semi_axes.full[[4]]
  
  
  compute_intersection_line_conic=function(a,b,mid,m,A){

    cross_start_x=rep(NA,length(m))
    cross_end_x=rep(NA,length(m))
    cross_start_y=rep(NA,length(m))
    cross_end_y=rep(NA,length(m))
    for(i in 1:length(m)){
      C <- ellipseToConicMatrix(c(a,b),c(mid),A)
      if(m[i]>-Inf & m[i]<Inf){
        l= c(m[i],-1,0)
        p_Cl <- intersectConicLine(C,l)
      }
      
      if(m[i]==Inf | m[i]==-Inf){ l=c(1,0,0)
      p_Cl <- intersectConicLine(C,l)
      }
      if(m[i]>=0){
        cross_end_x[i]=max(p_Cl[1,1],p_Cl[1,2])
        cross_start_x[i]=min(p_Cl[1,1],p_Cl[1,2])
        cross_end_y[i]=max(p_Cl[2,1],p_Cl[2,2])
        cross_start_y[i]=min(p_Cl[2,1],p_Cl[2,2])
      }
      if(m[i]<0){
        cross_end_x[i]=max(p_Cl[1,1],p_Cl[1,2])
        cross_start_x[i]=min(p_Cl[1,1],p_Cl[1,2])
        cross_end_y[i]=min(p_Cl[2,1],p_Cl[2,2])
        cross_start_y[i]=max(p_Cl[2,1],p_Cl[2,2])
      }

    }
    ret=list(cross_start_x,cross_end_x,cross_start_y,cross_end_y)
    return(ret)
  }
  
  
  
  #intersections
  ####################################################

  cross_start_x.full=NULL
  cross_start_y.full=NULL
  cross_end_x.full=NULL
  cross_end_y.full=NULL
  
  proj_list=list()# x an y coordinates of the projections
  proj=matrix(0,length(p11[,1]),2)
  for(i in 1:length(m)){
    for(j in 1:length(p11[,1])){
      proj[j,]=t(axes[i,]%*%((t(axes[i,])%*%p11[j,])/(t(axes[i,])%*%axes[i,])))
    }
    proj_list=append(proj_list,list(proj),after=length(proj_list))
    
  }
  for(i in 1:length(m)){
    if(m[i]>=0){
      cross_start_x.full=append(cross_start_x.full,min(proj_list[[i]][,1]),after = length(cross_start_x.full))
      cross_start_y.full=append(cross_start_y.full,min(proj_list[[i]][,2]),after = length(cross_start_y.full))
      cross_end_x.full=append(cross_end_x.full,max(proj_list[[i]][,1]),after = length(cross_end_x.full))
      cross_end_y.full=append(cross_end_y.full,max(proj_list[[i]][,2]),after=length(cross_end_y.full))
    }
    if(m[i]<0){
      cross_start_x.full=append(cross_start_x.full,min(proj_list[[i]][,1]),after = length(cross_start_x.full))
      cross_start_y.full=append(cross_start_y.full,max(proj_list[[i]][,2]),after = length(cross_start_y.full))
      cross_end_x.full=append(cross_end_x.full,max(proj_list[[i]][,1]),after = length(cross_end_x.full))
      cross_end_y.full=append(cross_end_y.full,min(proj_list[[i]][,2]),after=length(cross_end_y.full))
    }
  }
  
  comp_int=compute_intersection_line_conic(a,b,mid,m,A)
  cross_start_x=comp_int[[1]] #for percentile
  cross_end_x=comp_int[[2]]
  cross_start_y=comp_int[[3]]
  cross_end_y=comp_int[[4]]
  
  move=function(x_start,x_end,y_start,y_end,m,d){ 
    x=c(x_start,x_end)
    y=c(y_start,y_end)
    aa=d/sqrt((1/m)^2+1)+x
    bb=-1/m*(aa-x)+y
    new=cbind(aa,bb)
    return(new)
  }
  #these are all ordered
  x=cbind(cross_start_x,cross_end_x)
  x.full=cbind(cross_start_x.full,cross_end_x.full)
  y=cbind(cross_start_y,cross_end_y)
  y.full=cbind(cross_start_y.full,cross_end_y.full)
  
  
  # Given three colinear points p, q, r, the function checks if  
  # point q lies on line segment 'pr'  
  onSegment=function(p, q, r){
    if ( (q[1] <= max(p[1], r[1])) & (q[1] >= min(p[1], r[1])) & 
         (q[2] <= max(p[2], r[2])) & (q[2] >= min(p[2], r[2])))
      return (TRUE)
    else
      return (FALSE)
  }
  
  orientation=function(p, q, r){ 
    # to find the orientation of an ordered triplet (p,q,r) 
    # function returns the following values: 
    # 0 : Colinear points 
    # 1 : Clockwise points 
    # 2 : Counterclockwise 

    val = ((q[2] - p[2]) * (r[1] - q[1])) - ((q[1] - p[1]) * (r[2] - q[2])) 
    if (val > 0) 
      # Clockwise orientation 
      return (1)
    if (val < 0) 
      # Counterclockwise orientation 
      return (2)
    else
      # Colinear orientation 
      return (0)
  }
  # The main function that returns true if  
  # the line segment 'p1q1' & 'p2q2' intersect. 
  doIntersect=function(p1,q1,p2,q2){ 
    # Find the 4 orientations required for  
    # the general & special cases 
    o1 = orientation(p1, q1, p2) 
    o2 = orientation(p1, q1, q2) 
    o3 = orientation(p2, q2, p1) 
    o4 = orientation(p2, q2, q1) 
    # General case 
    if ((o1 != o2) & (o3 != o4)){
      return (TRUE)
    }
    
    # Special Cases 
    # p1 , q1 & p2 are colinear & p2 lies on segment p1q1 
    if ((o1 == 0) & onSegment(p1, p2, q1)) 
      return (TRUE)
    # p1 , q1 & q2 are colinear & q2 lies on segment p1q1 
    if ((o2 == 0) & onSegment(p1, q2, q1)) 
      return (TRUE)
    # p2 , q2 & p1 are colinear & p1 lies on segment p2q2 
    if ((o3 == 0) & onSegment(p2, p1, q2)) 
      return (TRUE)
    # p2 , q2 & q1 are colinear & q1 lies on segment p2q2 
    if ((o4 == 0) & onSegment(p2, q1, q2))
      return (TRUE)
    # If none of the cases 
    return (FALSE)
  }
  #doIntersect(p1, q1, p2, q2); line 1 is p1q1 , line 2 is p2q2
  
  
  storehouse=array(0,dim=c(2,2,length(m)))
  storehouse0=array(0,dim=c(2,2,length(m)))
  dimnames(storehouse)=list(c('',''),c('',''),colnames(p1.unscaled))
  move_lines_out_of_ellipse=function(m,x,y,x.full,y.full,
                                     e,t,a,b,a.full,b.full){
    
    
    order_m_start=order(m,decreasing = FALSE)
    dimnames(storehouse0)=list(c('',''),c('',''),colnames(p1.unscaled)[order_m_start])
    m0=m[order_m_start]
    x0=x[order_m_start,]
    y0=y[order_m_start,]
    x0.full=x.full[order_m_start,]
    y0.full=y.full[order_m_start,]

    if(use_percentile_for_plotting==TRUE){
      a_use=a 
      b_use=b
      mid_use=mid
      A_use=A
    }else{
      a_use=a.full 
      b_use=b.full
      mid_use=mid.full
      A_use=A.full
    }
    if(tan(A)>1 & tan(A)<(-1)){
      a_use=a.temp
      a_use=b_use
      b_use=a.temp 
    }

    
    d_list0=NULL
    for (i in 1:length(m)){
      d0=max(a_use,b_use)+0.3*min(a_use,b_use)
      
      if(i%%2==0 ){
        d0=d0
      }
      if(i%%2!=0){
        d0=-d0
      }
      s_1=sign(d0)
      new_=move(x0.full[i,1],x0.full[i,2],y0.full[i,1],y0.full[i,2],m0[i],d0)
      storehouse0[,,i]=new_
      
      d_list0=append(d_list0,d0,length(d_list0))

      if(i>1){
        d0=max(a_use,b_use)+0.3*min(a_use,b_use)
        d0=s_1*d0
        #for every previous line
        #reset the d value
        obj=rep(0,length(m0))
        for(jj in 1:(i-1)){
          obj[jj]=doIntersect(p1=storehouse0[1,,jj], q1=storehouse0[2,,jj], p2=storehouse0[1,,i], q2=storehouse0[2,,i])
        }
        counter=0
        while(any(obj==TRUE)){ 
          counter=counter+1
          
          d_list0=d_list0[-i]
          s_=sign(d0)

          
          if(counter%%2==0){
            d0=abs(d0) 
            d0=d0+move_lines_dist
            d0=-s_*d0 
          }

          if(counter%%2!=0){
            d0=-d0
          }
          
      
          new_=move(x0.full[i,1],x0.full[i,2],y0.full[i,1],y0.full[i,2],m0[i],d0)
          storehouse0[,,i]=new_
 
          
          d_list0=append(d_list0,d0,length(d_list0))
          obj=rep(0,length(i-1))
          for(jj in 1:(i-1)){
            obj[jj]=doIntersect(p1=storehouse0[1,,jj], q1=storehouse0[2,,jj], p2=storehouse0[1,,i], q2=storehouse0[2,,i])
          }
        }
        
      }
    }

    d_list=d_list0[order(order_m_start)]
    assign("d_list",d_list,envir = globalenv())
    assign("storehouse0",storehouse0,envir = globalenv())
    assign("order_m_start",order_m_start,envir = globalenv())
    storehouse=storehouse0[,,order(order_m_start)]
   
  }
  storehouse=move_lines_out_of_ellipse(m,x,y,x.full,y.full, 
                                       e,t,a,b,a.full,b.full)
  
  #project onto axes

  rotate <- function(point, theta, degree = F) {
    if (degree) theta <- theta * pi / 180
    rotate.matrix <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), byrow = T, nrow = 2)
    rotate.point <- rotate.matrix %*% point
    rotate.point
  }
  proj_list=list()
  proj=matrix(0,length(p11[,1]),2)
  for(i in 1:length(m)){
    for(j in 1:length(p11[,1])){
      proj[j,]=t(axes[i,]%*%((t(axes[i,])%*%p11[j,])/(t(axes[i,])%*%axes[i,])))
    }
    proj_list=append(proj_list,list(proj),after=length(proj_list))
  }
  proj_list_hor=list() 
  proj_hor=matrix(0,length(p11[,1]),2)
  for(i in 1:length(m)){
    for(j in 1:length(p11[,1])){
      proj_hor[j,]=t(rotate (proj_list[[i]][j,], -atan(m[i]), degree = F))
    }
    proj_list_hor=append(proj_list_hor,list(proj_hor),after=length(proj_list_hor))
  }
  proj_onto_var=matrix(0,length(p11[,1]),length(m))
  for(i in 1:length(m)){
    proj_onto_var[,i]=proj_list_hor[[i]][,1] 
  }
  colnames(proj_onto_var)=paste0('Var',1:length(PCA$axes[,1]))
  
  
  
  #densities
  rotate <- function(point, theta, degree = F) {
    if (degree) theta <- theta * pi / 180 #convert to radians
    rotate.matrix <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), byrow = T, nrow = 2)
    rotate.point <- rotate.matrix %*% point
    rotate.point
  }
  
  nr_of_breaks=round(sqrt(length(p1.unscaled[,1])),0)
  if(is.null(cluster_var)==FALSE){
    grps=sort(unique(cluster_var),decreasing=FALSE)}
  if(is.null(cluster_var)==TRUE){
    grps=1}
  D_x=list()
  D_y=list()
  for(i in 1:length(m)){
    for(vv in 1:length(grps)){
      #generate densities
      #rotate them
      #move them out towards their lines
      rotated_x=NULL
      rotated_y=NULL
      
      if(length(grps)>1){
        if(counts_or_densities=='densities'){
          if(sum(cluster_var==grps[vv])<=1){
            stop('Too few points per cluster to produce a distribution by cluster. Either set plot_the_dens=FALSE or cluster_var=NULL')}
          size_of_br=((max(proj_onto_var[,i][cluster_var==grps[vv]])-min(proj_onto_var[,i][cluster_var==grps[vv]])))/nr_of_breaks
          br=seq(min(proj_onto_var[,i][cluster_var==grps[vv]]),max(proj_onto_var[,i][cluster_var==grps[vv]]),by=size_of_br)
          den_=hist(proj_onto_var[,i][cluster_var==grps[vv]],breaks=br,plot=FALSE)$density*(fact_)
          mids_=hist(proj_onto_var[,i][cluster_var==grps[vv]],breaks=br,plot=FALSE)$mids
          den=c(0,den_,0)
          size_of_jump=mids_[2]-mids_[1]
          mids=c(mids_[1]-size_of_jump,mids_,mids_[length(mids_)]+size_of_jump)
          if(smoothing_method=='smooth.spline'){
            smoothingSpline = smooth.spline(mids, den, spar=smooth_par)
            #to make sure the smoothing par for the first distribution is used for all the densities
            if(is.null(smooth_par)& i==1 & vv==1){
              smooth_par=smoothingSpline$spar
            }
            mids=smoothingSpline$x
            den=smoothingSpline$y
          }
          
        }
        
        if(counts_or_densities=='counts'){
          if(sum(cluster_var==grps[vv])<=1){
            stop('Too few points per cluster to produce a distribution by cluster. Either set plot_the_dens=FALSE or cluster_var=NULL')}
          size_of_br=((max(proj_onto_var[,i][cluster_var==grps[vv]])-min(proj_onto_var[,i][cluster_var==grps[vv]])))/nr_of_breaks
          br=seq(min(proj_onto_var[,i][cluster_var==grps[vv]]),max(proj_onto_var[,i][cluster_var==grps[vv]]),by=size_of_br)
          den_=hist(proj_onto_var[,i][cluster_var==grps[vv]],breaks=br,plot=FALSE)$counts*fact_/10
          mids_=hist(proj_onto_var[,i][cluster_var==grps[vv]],breaks=br,plot=FALSE)$mids
          den=c(0,den_,0)
          size_of_jump=mids_[2]-mids_[1]
          mids=c(mids_[1]-size_of_jump,mids_,mids_[length(mids_)]+size_of_jump)
          if(smoothing_method=='smooth.spline'){
            smoothingSpline = smooth.spline(mids, den, spar=smooth_par)
            #to make sure the smoothing par for the first distribution is used for all the densities
            if(is.null(smooth_par)& i==1 & vv==1){
              smooth_par=smoothingSpline$spar
            }#
            mids=smoothingSpline$x
            den=smoothingSpline$y
          }
        }
      }
      if(length(grps)<=1){
        if(counts_or_densities=='densities'){
          size_of_br=((max(proj_onto_var[,i])-min(proj_onto_var[,i])))/nr_of_breaks
          br=seq(min(proj_onto_var[,i]),max(proj_onto_var[,i]),by=size_of_br)
          den_=hist(proj_onto_var[,i],breaks=br,plot=FALSE)$density*(fact_)
          mids_=hist(proj_onto_var[,i],breaks=br,plot=FALSE)$mids
          den=c(0,den_,0)
          size_of_jump=mids_[2]-mids_[1]
          mids=c(mids_[1]-size_of_jump,mids_,mids_[length(mids_)]+size_of_jump)
          if(smoothing_method=='smooth.spline'){
            smoothingSpline = smooth.spline(mids, den, spar=smooth_par)
            #to make sure the smoothing par for the first distribution is used for all the densities
            if(is.null(smooth_par)& i==1 & vv==1){
              smooth_par=smoothingSpline$spar
            }#
            mids=smoothingSpline$x
            den=smoothingSpline$y
          }
        }
        if(counts_or_densities=='counts'){
          size_of_br=((max(proj_onto_var[,i])-min(proj_onto_var[,i])))/nr_of_breaks
          br=seq(min(proj_onto_var[,i]),max(proj_onto_var[,i]),by=size_of_br)
          den_=hist(proj_onto_var[,i],breaks=br,plot=FALSE)$counts*(fact_)/10
          mids_=hist(proj_onto_var[,i],breaks=br,plot=FALSE)$mids
          den=c(0,den_,0)
          size_of_jump=mids_[2]-mids_[1]
          mids=c(mids_[1]-size_of_jump,mids_,mids_[length(mids_)]+size_of_jump)
          if(smoothing_method=='smooth.spline'){
            smoothingSpline = smooth.spline(mids, den, spar=smooth_par)
            #to make sure the smoothing par for the first distribution is used for all the densities
            if(is.null(smooth_par)& i==1 & vv==1){
              smooth_par=smoothingSpline$spar
            }#
            mids=smoothingSpline$x
            den=smoothingSpline$y
          }
        }
      }
      
      
      rot_x=rep(0,length(den))
      rot_y=rep(0,length(den))
      for(jj in 1:length(den)){
        rotated=rotate(c(mids[jj],den[jj]),atan(m[i]))
        rot_x[jj]=rotated[1,]
        rot_y[jj]=rotated[2,]
      }
      moved_rotated_x=d_list[i]/sqrt((1/m[i])^2+1)+rot_x
      moved_rotated_y=-1/m[i]*(moved_rotated_x-rot_x)+rot_y
      D_x=append(D_x,list(moved_rotated_x),after=length(D_x))
      D_y=append(D_y,list(moved_rotated_y),after=length(D_y))
      
    }
    assign('D_x',D_x,envir=globalenv())
    assign('D_y',D_y,envir=globalenv())
  }
  
  ticks=array(0,dim = c(length(c(-1:-100,0,1:100)),3,ncol(p1.unscaled)))
  #create tick marks
  for(i in 1:ncol(p1.unscaled)){
    x_=c(rev((-1:-100)/(draw_more_ticks)*PCA$len_1_unit[1,i]),0,(1:100)/(draw_more_ticks)*PCA$len_1_unit[1,i])
    y_=c(rev((-1:-100)/(draw_more_ticks)*PCA$len_1_unit[2,i]),0,(1:100)/(draw_more_ticks)*PCA$len_1_unit[2,i])
    tick=c(rev(-1:-100)/(draw_more_ticks),0,1:100/(draw_more_ticks))
    ticks[,,i]=cbind(x_,y_,tick)
  }
  
  
  tick_x=list()
  tick_y=list()
  for(i in 1:length(m)){
    #move them out towards their lines
    moved_rotated_x=d_list[i]/sqrt((1/m[i])^2+1)+ticks[,1,i]
    moved_rotated_y=-1/m[i]*(moved_rotated_x-ticks[,1,i])+ticks[,2,i]
    tick_x=append(tick_x,list(moved_rotated_x),after=length(tick_x))
    tick_y=append(tick_y,list(moved_rotated_y),after=length(tick_y))
  }
  I_=list()
  for(i in 1:length(m)){
    I_=append(I_,list(tick_x[[i]]<=max(storehouse[,1,i]) & tick_x[[i]]>=min(storehouse[,1,i])&(tick_y[[i]]<=max(storehouse[,2,i]) & tick_y[[i]]>=min(storehouse[,2,i]))),length(I_))
  }
  
  find_proj_line=function(x1,y1,x2,y2,m){
    #p1.unscaled(x,y)=x1 and y1
    #variable line (x,y)= x2 and y2 and m
    a1=-(-1/m)
    b1=1
    c1=(-1/m)*(-x1)+y1
    
    a2=-m
    b2=1
    c2=m*(-x2)+y2
    
    int_x=(c1*b2-c2*b1)/(a1*b2-a2*b1)
    int_y=(a2*c1-a1*c2)/(a2*b1-a1*b2)
    int=c(int_x,int_y)
    return(int)
  }
  
  
  if(use_convex_hull==TRUE){
    if(use_percentile_for_plotting==TRUE){
      hpts <- chull(p1)
      hpts <- c(hpts, hpts[1])
      hpts=p1[hpts, ]
    }
    if(use_percentile_for_plotting==FALSE){
      hpts <- chull(p11)
      hpts <- c(hpts, hpts[1])
      hpts=p11[hpts, ]
      
    }
  }
  
  
  #plotly
  #get dataset ready
  p11=as.data.frame(p11)
  colnames(p11)<-c('point_x','point_y')
  if(is.null(row.names(p1.unscaled))){
    names_m=paste0(1:length(p1.unscaled[,1]))
  }else {
    rownames(p11)=rownames(p1.unscaled)
  }
  
  if(is.null(colnames(p1.unscaled))){
    c_names=paste0('Var',1:length(m))
  }else {
    c_names=colnames(p1.unscaled)
  }
  
  storehouse<-as.data.frame(storehouse)
  for(i in 1:(length(m))*2){
    colnames(storehouse)[i-1]=paste0('x',i/2)
    colnames(storehouse)[i]=paste0('y',i/2)
  } 
  ticks=as.data.frame(ticks[,3,1])
  colnames(ticks)=c('v')
  I_=array(I_)
  dimnames(I_)=list(paste(colnames(p1.unscaled)))
  
  
  see_projections=function(name=NULL){
    proj_intersect=matrix(0,nrow=length(m),ncol=2)
    if(is.numeric(name) & !is.null(name)){
      val=proj_onto_var[name,]
      for(i in 1:length(m)){
        proj_intersect[i,]=find_proj_line(PCA$points[name,1],PCA$points[name,2],
                                          storehouse[1,(i*2)-1],storehouse[1,i*2],m[i])
      }
    }
    if(!is.numeric(name) & !is.null(name)){
      val=proj_onto_var[paste(name),]
      for(i in 1:length(m)){
        proj_intersect[i,]=find_proj_line(PCA$points[paste(name),1],PCA$points[paste(name),2],
                                          storehouse[1,i*2-1],storehouse[1,i*2],m[i])
        
      }
    }
    return(proj_intersect)
  }
  
  storehouse_seq=matrix(0,nrow=20,ncol=length(storehouse[1,]))
  storehouse_seq=as.data.frame(storehouse_seq)
  for(i in 1:(length(m))*2){
    colnames(storehouse_seq)[i-1]=paste0('x',i/2)
    colnames(storehouse_seq)[i]=paste0('y',i/2)
  } 
  for(i in 1:(length(m)*2)){
    storehouse_seq[,i]=seq(from=storehouse[1,i],storehouse[2,i],length.out = 20)
  }
  
  
  proj_lines=array(0,dim=c(length(m),2,length(p1.unscaled[,1])))
  for(ii in 1:length(p1.unscaled[,1])){
    proj_lines[,,ii]=see_projections(name=(1:length(p1.unscaled[,1]))[ii])

  }
  dimnames(proj_lines)=list(c(paste(colnames(p1.unscaled))),c('x','y'),c(paste(rownames(p1.unscaled))))
  
  
  if(is.null(cluster_var)==FALSE){
    if(use_convex_hull==TRUE){
      dataset=list(points=cbind(p11,cluster_var),lines=storehouse,D_x=D_x,D_y=D_y,t.full=t.full,e.full=e.full,
                   t=t,e=e,projected_val=proj_onto_var,proj_lines=proj_lines,
                   tick_x=tick_x,tick_y=tick_y,ticks=ticks,I_=I_,hpts=hpts,s_p1=s_p1,m_p1=m_p1,storehouse_seq=storehouse_seq)
    }
    if(use_convex_hull==FALSE){
      dataset=list(points=cbind(p11,cluster_var),lines=storehouse,D_x=D_x,D_y=D_y,t.full=t.full,e.full=e.full,
                   t=t,e=e,projected_val=proj_onto_var,proj_lines=proj_lines,
                   tick_x=tick_x,tick_y=tick_y,ticks=ticks,I_=I_,s_p1=s_p1,m_p1=m_p1,storehouse_seq=storehouse_seq)
    }
    
  }
  if(is.null(cluster_var)==TRUE){
    if(use_convex_hull==TRUE){
      cluster_var=0
      dataset=list(points=cbind(p11,cluster_var),lines=storehouse,D_x=D_x,D_y=D_y,t.full=t.full,e.full=e.full,
                   t=t,e=e,projected_val=proj_onto_var,proj_lines=proj_lines,
                   tick_x=tick_x,tick_y=tick_y,ticks=ticks,I_=I_,hpts=hpts,s_p1=s_p1,m_p1=m_p1,storehouse_seq=storehouse_seq)
    }
    if(use_convex_hull==FALSE){
      cluster_var=0
      dataset=list(points=cbind(p11,cluster_var),lines=storehouse,D_x=D_x,D_y=D_y,t.full=t.full,e.full=e.full,
                   t=t,e=e,projected_val=proj_onto_var,proj_lines=proj_lines,
                   tick_x=tick_x,tick_y=tick_y,ticks=ticks,I_=I_,s_p1=s_p1,m_p1=m_p1,storehouse_seq=storehouse_seq)
    }
    
  }
  
  #Add densities here to the dataset
  
  if(gray_scale==TRUE){
    col_m=gray.colors(length(m)+length(grps), start = 0.45, end = 0.8, gamma = 2.2, alpha = NULL)
  }
  else{
    pal = c(
      rgb(57 / 225, 106 / 225, 177 / 225),
      rgb(218 / 225, 124 / 225, 48 / 225),
      rgb(62 / 225, 150 / 225, 81 / 225),
      rgb(204 / 225, 37 / 225, 41 / 225),
      rgb(83 / 225, 81 / 225, 84 / 225),
      rgb(107 / 225, 76 / 225, 154 / 225),
      rgb(146 / 225, 36 / 225, 40 / 225),
      rgb(148 / 225, 139 / 225, 61 / 225)
    )
    pal2=c(
      rgb(114 / 225, 147 / 225, 203 / 225),
      rgb(225 / 225, 151 / 225, 76 / 225),
      rgb(132 / 225,186 / 225, 91 / 225),
      rgb(211 / 225, 94 / 225, 96 / 225),
      rgb(128 / 225, 133 / 225, 133 / 225),
      rgb(144/ 225, 103 / 225, 167 / 225),
      rgb(171 / 225, 104 / 225, 87 / 225),
      rgb(204 / 225, 194/ 225, 16 / 225)
      
    )
    

    col_m=c(pal,pal2)[1:(length(grps))]
    
  }
  
  
  p <- plot_ly()
  
  #DATAPOINTS
  col_seq=seq(1,(length(m)+length(grps)),length.out = length(grps))
  if(plot_all_the_data==TRUE){
    if(length(cluster_var)==1)leg_entree='data'
    if(length(cluster_var)>1) leg_entree=paste(grps)
    
    if(is.null(rownames(p1.unscaled))){
      if(length(cluster_var)==1){
        p= add_markers(p,data=dataset, x = ~round(dataset$points[,1],3), y = ~round(dataset$points[,2],3),
                       color=I(col_m[1]),colors=col_m,
                       type='scatter',mode='markers',
                       marker=list(size=4),showlegend=TRUE,name='data',
                       hoverinfo='text',text=~paste(1:length(dataset$points[,1])))
      }else{
        p= add_markers(p,data=dataset, x = ~round(dataset$points[,1],3), y = ~round(dataset$points[,2],3),color=~as.factor(dataset$points[,3]),
                       colors=col_m,
                       type='scatter',mode='markers',
                       marker=list(size=4),showlegend=TRUE,name=paste(unique(cluster_var)),
                       hoverinfo='text',text=~paste(1:length(dataset$points[,1])))
      }
    }
    else{
      if(length(cluster_var)==1){
        p= add_markers(p,data=dataset, x = ~round(dataset$points[,1],3), y = ~round(dataset$points[,2],3),
                       color=I(col_m[1]),colors=col_m,
                       type='scatter',mode='markers',
                       marker=list(size=4),showlegend=TRUE,name='data',
                       hoverinfo='text',text=~paste(rownames(p1.unscaled)))
      }
      else 
        p= add_markers(p,data=dataset, x = ~round(dataset$points[,1],3), y = ~round(dataset$points[,2],3),color=~as.factor(dataset$points[,3]),
                       colors=col_m,
                       type='scatter',mode='markers',
                       marker=list(size=4),showlegend=TRUE,
                       hoverinfo='text',text=~paste(rownames(p1.unscaled)))
    }                 
  }
  
  
  #VARIABLE LINES
  m_temp=round(m,3)
  for(i in 1:length(m)){
    x = paste0("x",i)
    y = paste0("y",i)
    te = paste0(c_names[i]," (",PCA$axqual[i]*100,"%)")
    n = paste0(c_names[i]," (",PCA$axqual[i]*100,"%)")
    
    p <- add_trace(p, 
                   x = dataset$storehouse_seq[[x]], 
                   y = dataset$storehouse_seq[[y]], 
                   type = 'scatter', 
                   mode = 'lines', 
                   line = list( width = 1.4,color="darkolivegreen"),colors=col_m,name=n,hoverinfo='text',text=paste(te)) #list(color = (brewer.pal(n = 3, name = "RdBu"))
  }
  
  
  
  #DENSITY
  if(plot_the_dens==TRUE){

    for(i in 1:(length(m)*length(grps))){
      x = i
      y = i
      if(length(grps)>1){
        if(!is.null(colnames(p1.unscaled))){#if there is colnames
          n=paste0(unique(cluster_var)[rep((1:length(grps)),length(m))[i]],', ',colnames(p1.unscaled)[ceiling((1:(length(grps)*length(m)))/length(grps))][i])
        }
        if(is.null(colnames(p1.unscaled))){
          n=paste0(unique(cluster_var)[rep((1:length(grps)),length(m))[i]],', ','Var',ceiling(i/length(grps)))
        }
      }
      if(length(grps)==1){
        if(!is.null(colnames(p1.unscaled))){
          n=paste0(colnames(p1.unscaled)[i],' density')
        }
        if(!is.null(colnames(p1.unscaled))){
          n=paste0('Var',ceiling(i/length(grps)),' density')
        }
        
      }
      
      if(length(length(grps)>1)){
        C=col_m[rep(1:length(grps),length(m))[i]]
      }
      if(length(grps)==1) {C=col_m[1]}
      
      te= paste0(c_names[ceiling(i/length(grps))],', m=',m_temp[ceiling(i/length(grps))])
      p <- add_trace(p, 
                     x = dataset$D_x[[x]], 
                     y = dataset$D_y[[y]], 
                     type = 'scatter', 
                     mode = 'lines',
                     showlegend=TRUE,line = list(width =1.5 ,dash='dot',color=paste(C)),name=n ,hoverinfo='text',text=n)#
    }
    # }
  }
  
  #ELLIPSE
  if(plot_ellipse==TRUE){
    p<-add_trace(p,dataset,x=~t.full,y=~e.full,type='scatter',mode='lines',
                 line = list(width = 1.25,color =gray(0.6),dash='dash'),showlegend=FALSE,hoverinfo='text')
  }
  if(plot_percentile_ellipse==TRUE){
    p<-add_trace(p,dataset,x=~t,y=~e,type='scatter',mode='lines',
                 line = list(width = 1.25,color =gray(0.6),dash='dash'),showlegend=FALSE,hoverinfo='text')
  }
  #LAYOUT
  p<-layout(p,xaxis=list(title='',showticklabels = FALSE,zeroline=FALSE,showgrid = FALSE),
            yaxis=list(title=PCA$quality,showticklabels = FALSE,zeroline=FALSE,scaleanchor={'x'},scaleratio=1,showgrid = FALSE), 
            title=list('Biplot'),plot_bgcolor='rgba(0,0,0,0)',
            annotations=list((color=col_m[3])))
  
  
  #TICK MARKS
  tick_list=list()
  tick_position_x=list()
  tick_position_y=list()
  predicted_val=rep(0,length(m))
  for(i in 1:length(m)){
    ind=i
    
    x=dataset$tick_x[[ind]]
    x_n=x[dataset$I_[[ind]]]
    y=dataset$tick_y[[ind]]
    y_n=y[dataset$I_[[ind]]]
    C_=col_m[-col_seq][i]
    if(ticks_in_original_units==TRUE){
      tick=dataset$ticks*dataset$s_p1[i]+dataset$m_p1[i]
    }
    if(ticks_in_original_units==FALSE){
      tick=dataset$ticks
    }
    
    tick_n=round(tick[I_[[ind]],1],round_ticks_to)
    tick_list=append(tick_list,list(tick_n),length(tick_list))
    tick_position_x=append(tick_position_x,list(x_n),after=length(tick_position_x))
    tick_position_y=append(tick_position_y,list(y_n),length(tick_position_y))
    
    p<-add_text(p,x=x_n,y=y_n,text=tick_n,textposition='top right',showlegend=FALSE,size=I(tick_size),hoverinfo='text')
    p<-add_markers(p,x=x_n,y=y_n,size=I(6),showlegend=FALSE,hoverinfo='text',marker=list(color=I(gray(0.3)),
                                                                                         line = list(width = 0)))
   
    if(is.numeric(name)==TRUE& !is.null(name))
      predicted_val[i]=(proj_lines[i,1,name]-x_n[1])*(tick_n[length(tick_n)]-tick_n[1])/(x_n[length(x_n)]-x_n[1])+tick_n[1]
    if(is.numeric(name)==FALSE & !is.null(name)){
      predicted_val[i]=(proj_lines[i,1,paste(name)]-x_n[1])*(tick_n[length(tick_n)]-tick_n[1])/(x_n[length(x_n)]-x_n[1])+tick_n[1]
    }
    if(is.null(name)){
      predicted_val=NULL
    }
  }
  
  for(i in 1:length(m)){
    x=paste0('x',i)
    y=paste0('y',i)
    ind=i
    n=paste0(c_names[i])
    if(plotqual == TRUE) texpos <- "top center"
    if(plotqual == FALSE) texpos <- "bottom center"
    tick_position_x[[i]][which.max(tick_list[[i]])]
    dist_first=sqrt((tick_position_x[[i]][which.max(tick_list[[i]])]-dataset$lines[[x]][1])^2+(tick_position_y[[i]][which.max(tick_list[[i]])]-dataset$lines[[y]][1])^2)#######################################!!!!!!!!!!!!!!!!
    dist_last=sqrt((tick_position_x[[i]][which.max(tick_list[[i]])]-dataset$lines[[x]][2])^2+(tick_position_y[[i]][which.max(tick_list[[i]])]-dataset$lines[[y]][2])^2)
    if(min(dist_first,dist_last)==dist_first){
      if(var_names==TRUE){
        p=add_text(p,x=dataset$lines[[x]][1],
                   y=dataset$lines[[y]][1],text=paste(n),textposition=texpos,textfont=list(size=variable_name_size),showlegend=FALSE,hoverinfo='text')
      }
      if(var_names==FALSE){
        p=add_markers(p,type='scatter',dataset$lines[[x]][1],dataset$lines[[y]][1],mode='markers',marker = list(symbol=(3),color=(paste(col_m[3])),size=8),showlegend=FALSE,hoverinfo='text')
        
      }
      if(plotqual==TRUE){
        p=add_text(p,x=dataset$lines[[x]][1],
                   y=dataset$lines[[y]][1],text=paste0("(",PCA$axqual[i]*100,"%)"),textposition="bottom center",textfont=list(size=variable_name_size*0.7),showlegend=FALSE,hoverinfo='text')
      }
      
    }
    
    if(min(dist_first,dist_last)==dist_last){
      if(var_names==TRUE){
        p=add_text(p,x=dataset$lines[[x]][2],
                   y=dataset$lines[[y]][2],text=paste(n),textposition=texpos,textfont=list(size=variable_name_size),showlegend=FALSE,hoverinfo='text')#,textpositionsrc
        
      }
      if(var_names==FALSE){
        p=add_markers(p,type='scatter',dataset$lines[[x]][2],dataset$lines[[y]][2],mode='markers',marker = list(symbol=(3),color=(paste(col_m[3])),size=8),showlegend=FALSE,hoverinfo='text')
      }
      if(plotqual==TRUE){
        p=add_text(p,x=dataset$lines[[x]][2],
                   y=dataset$lines[[y]][2],text=paste0("(",PCA$axqual[i]*100,"%)"),textposition="bottom center",textfont=list(size=variable_name_size*0.7),showlegend=FALSE,hoverinfo='text')
      }
    }
    
  }
  
  #PROJECTION
  if(is.null(name)==FALSE){
    
    if(is.numeric(name)==FALSE){
      
      val=predicted_val
      
      p<-add_segments(p,dataset$proj_lines[,1,paste(name)],dataset$proj_lines[,2,paste(name)],
                      xend = rep(dataset$points[paste(name),1],length(m)),
                      yend=rep(dataset$points[paste(name),2],length(m)),line=list(width=0.7,dash='dot',color='black'),showlegend=TRUE)
      p<-add_markers(p,dataset$proj_lines[,1,paste(name)],dataset$proj_lines[,2,paste(name)],type='scatter',
                     mode='markers',marker=list(color='black',symbols=c('x')),
                     text = c(paste(round(val,round_ticks_to))),
                     hoverinfo = 'text',showlegend=TRUE,name='predicted values')
      p<-add_markers(p,dataset$points[paste(name),1],dataset$points[paste(name),2],marker=list(color='black'),
                     hoverinfo='text',text=~paste(name),name='chosen value')
    }
    else {
      
      val=predicted_val
      
      p<-add_segments(p,dataset$proj_lines[,1,name],dataset$proj_lines[,2,name],
                      xend = rep(dataset$points[name,1],length(m)),
                      yend=rep(dataset$points[name,2],length(m)),line=list(width=0.7,dash='dot',color='black'),showlegend=TRUE)
      p<-add_markers(p,dataset$proj_lines[,1,name],dataset$proj_lines[,2,name],type='scatter',
                     mode='markers',marker=list(color='black',symbols=c('x')),
                     text = c(paste(round(val,round_ticks_to))),
                     hoverinfo = 'text',showlegend=TRUE,name='predicted values')
      p<-add_markers(p,dataset$points[name,1],dataset$points[name,2],marker=list(color='black'),
                     hoverinfo='text',text=~paste(name))
    }
  }
  
  #CONVEX HULL
  if(use_convex_hull==TRUE){
    if(use_percentile_for_plotting==TRUE){
      p=add_paths(p,dataset$hpts[,1],dataset$hpts[,2],line=list(width = 1.25,color =gray(0.7),dash='dash'),showlegend=FALSE)
    }
    if(use_percentile_for_plotting==FALSE){
      p=add_paths(p,dataset$hpts[,1],dataset$hpts[,2],line=list(width = 1.25,color =gray(0.7),dash='dash'),showlegend=FALSE)
    }
  }
  if(is.null(colnames(p1.unscaled))==FALSE){
    if(!is.null(name)){
      names(val)=colnames(p1.unscaled)
    }
  }
  if(is.null(rownames(p1.unscaled))==FALSE){
    rownames(axes)=colnames(p1.unscaled)
    rownames(p11)=rownames(p1.unscaled)
  }
  if(is.null(name)){
    val=NULL
  }
  ret=list(points=p11,p=p,axes=axes,predicted_values=val,length_of_one_unit=PCA$len_1_unit,quality_of_biplot=PCA$quality)#
  
  return(ret)
}