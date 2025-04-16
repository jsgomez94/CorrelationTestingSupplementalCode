
## Colorblind friendly palette. Required to use these colors for assignments.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", 
               "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7")
plot(1:8, col = cbPalette, pch = 19)

method_names <- c(
  "Cai & Zhang (2016)",
  "Zheng et al. (2019) sum",
  "Zheng et al. (2019) max",
  "Permutation Test")
n_vals <- seq(50, 600, by = 50)


par(cex.main = 1.5, cex.lab = 1)
layout(mat = matrix(c(1, 2), ncol = 2), widths = c(1,0.5))
# table <- read.csv(paste0("42_Results_average.csv"))
table <- read.csv(paste0("42_Results_average_400.csv"))

plot(
  x = n_vals, y = table[1, -1],
  ylim = c(0,1), xlim = c(0, 600),
  main = "",
  pch = 0, col = cbPalette[4],
  cex.main = 2, cex.lab = 1.5,
  xlab = "Total Sample Size (N)",
  ylab = "Proportion of Rejection")
lines(
  x = n_vals, y = table[1, -1],
  lty = 2, col = cbPalette[4])


abline(h = 0.05, col = "red", lty = 3)
points(
  x = n_vals, y = rep(0.05, length(n_vals)),
  pch =".", col = "red")


points(
  x = n_vals, y = table[2, -1],
  pch = 8, col = cbPalette[6])
lines(
  x = n_vals, y = table[2, -1],
  lty = 3, col = cbPalette[6]
)

points(
  x = n_vals, y = table[3, -1],
  pch = 10, col = cbPalette[7])
lines(
  x = n_vals, y = table[3, -1],
  lty = 4, col = cbPalette[7]
)

points(
  x = n_vals, y = table[4, -1],
  pch = 19, col = "#000000")
lines(
  x = n_vals, y = table[4, -1],
  lty = 1, col = "#000000"
) 

plot.new()

legend(
  "center", legend = , c(method_names, "Target level 0.05"), col = c(cbPalette[c(4, 6, 7)], "#000000", "red"),
  lty = c(2, 3, 4, 1, 3), pch = c(0, 8, 10, 19, 19), cex = 1)

mtext(
  text = "Rejection Rate Under Null: Gene Expression Data",
  side = 3, outer = TRUE, cex = 2, line = -2.5)
