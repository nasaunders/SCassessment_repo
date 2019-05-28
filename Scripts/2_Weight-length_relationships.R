#-*- coding: latin-9 -*-

### File: 2_Weight-length_relationships.R
### Time-stamp: <2019-01-30 14:56:11 yreecht>
###
### Created: 12/09/2018	13:14:33
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################


res <- by(data = LWdata, INDICES = as.list(LWdata[ , "species", drop = FALSE]),
          FUN = function(x)
          {
              glmLW <- glm(log(weight) ~ log(size), data = x,
                           family = gaussian)

              newdata <- data.frame(size = seq(from = min(x$size, na.rm = TRUE),
                                               to = max(x$size, na.rm = TRUE),
                                               length.out = 10))

              glmPred <- predict.glm(object = glmLW, newdata = newdata, se.fit = TRUE)

              .glmRes[[unique(as.character(x$species))]] <<-
                  list(species = unique(as.character(x$species)),
                       glmObj = glmLW,
                       glmPred = glmPred)

              ## Back calculation of effective N used for the estimates of se.fit:
              Neff <- summary(glmLW)$dispersion / (glmPred$se.fit ^ 2)

              ## Correction of se.fit for back calculation of the log transformation (Cox method):
              se2BT <- sqrt(summary(glmLW)$dispersion / Neff +
                                          summary(glmLW)$dispersion ^ 2 / (2 * (Neff - 1)))

              sd2BT <- sqrt(summary(glmLW)$dispersion + summary(glmLW)$dispersion^2 / 2)

              a <- exp(glmLW$coefficients[1] + summary(glmLW)$dispersion / 2)
              b <- glmLW$coefficients[2]

              ## Check N effectifs:
              if (! all.equal(sqrt(summary(glmLW)$dispersion / Neff),
                              glmPred$se.fit)) warning("Weight-Length:",
                                                       " something wrong in the estimates of N effectif (",
                                                       round(Neff, 1),
                                                       ")")


              ## Saving graphics:
              for (devType in getOption("surveyPlot.dev"))
              {

                  if (devType %in% c("X11", "pdf", "jpg", "jpeg", "png"))
                  {
                      width <- c(png = 1000, jpg = 1000, jpeg = 1000, X11 = 10, pdf = 10)
                      dev <- openDev(device = devType,
                                     directory = ResultsPath,
                                     filename = paste0("Weight-Length_",
                                                       gsub("[[:blank:]]+", "_", unique(x$species))),
                                     width = width,
                                     height = width * 7 / 10,
                                     pointsize = c(png = 20, jpg = 20, jpeg = 20, X11 = 12, pdf = 12),
                                     counter = FALSE,
                                     verbose=FALSE)
                  }else{
                      warning("Device \"", devType, "\" not supported")
                      next()
                  }

                  par(mfrow = c(2, 3), oma = c(0, 0, 3.7, 0), mar = c(4.5, 4.5, 2, 1) + 0.1)
                  plot(glmLW)
                  mtext(text = bquote(italic(.(unique(as.character(x$species))))),
                        outer = TRUE, line = 1.8, cex = 1.2)

                  ## X11(width = 10)
                  ## par(mfrow = c(1, 2), oma = c(0, 0, 1.5, 0), mar = c(4.5, 4.5, 1, 1) + 0.1)

                  plot(log(weight) ~ log(size), data = x, col = "red")
                  lines(x = log(newdata$size), y = glmPred$fit, col = "blue", lwd = 2)
                  lines(x = log(newdata$size), y = glmPred$fit +
                                                   qt(p = 0.025, df = Neff - 1) * glmPred$se.fit,
                        col = "blue", lwd = 2, lty = 2)
                  lines(x = log(newdata$size), y = glmPred$fit +
                                                   qt(p = 0.975, df = Neff - 1) * glmPred$se.fit,
                        col = "blue", lwd = 2, lty = 2)
                  lines(x = log(newdata$size), y = glmPred$fit +
                                                   qnorm(p = 0.025) * glmPred$residual.scale,
                        col = "green", lwd = 2, lty = 2)
                  lines(x = log(newdata$size), y = glmPred$fit +
                                                   qnorm(p = 0.975) * glmPred$residual.scale,
                        col = "green", lwd = 2, lty = 2)

                  plot(weight ~ size, data = x, col = "blue", pch = 16)

                  lines(x = newdata$size, y = exp(glmPred$fit + summary(glmLW)$dispersion / 2), col = "red", lwd = 2)
                  lines(x = newdata$size, y = exp(glmPred$fit +
                                                  summary(glmLW)$dispersion / 2 +
                                                                qt(p = 0.025, df = Neff - 1) * se2BT),
                        col = "red", lwd = 2, lty = 2)
                  lines(x = newdata$size, y = exp(glmPred$fit +
                                                  summary(glmLW)$dispersion / 2 +
                                                                qt(p = 0.975, df = Neff - 1) * se2BT),
                        col = "red", lwd = 2, lty = 2)
                  lines(x = newdata$size, y = exp(glmPred$fit + qnorm(p = 0.025) * sd2BT),
                        col = "green", lwd = 2, lty = 2)
                  lines(x = newdata$size, y = exp(glmPred$fit + qnorm(p = 0.975) * sd2BT),
                        col = "green", lwd = 2, lty = 2)
                  ## mtext(text = bquote(italic(.(unique(as.character(x$species))))), outer = TRUE, line = 0)

                  mtext(text = bquote(expr = "W ="~.(format(a, digits = 3)) %*% L ^ .(format(b, digits = 3))),
                        side = 3, line = -1.9, col = "red", font = 2)

                  if (devType != "X11") dev.off()
              }

              ## Saving the coefficients:
              write.csv(matrix(data = c(a, b),
                               ncol = 2,
                               dimnames = list(unique(as.character(x$species)),
                                               c("a (bias corrected)", "b"))),
                        file = file.path(ResultsPath,
                                         paste0("Weight-Length_",
                                                gsub("[[:blank:]]+", "_", unique(x$species)),
                                                ".csv")),
                        row.names = TRUE)

              return(matrix(data = c(a, b),
                            ncol = 2,
                            dimnames = list(unique(as.character(x$species)),
                                            c("a", "b"))))
          },
          simplify = FALSE)

## Binding rows together:
LWparamCalc <- do.call(rbind, res)







### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
