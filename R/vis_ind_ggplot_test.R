# 
# vF <- rbind(c(0.4, 0.2, 0.1, 0.1, 0.1, 0.1),
#             c(0.2, 0.4, 0.1, 0.3, 0, 0),
#             c(0.05, 0.05, 0.05, 0.85, 0, 0),
#             c(0.3, 0.4, 0.1, 0.2, 0, 0),
#             c(0.1, 0.3, 0.4, 0.2, 0, 0))
# 
# 
# if (is.na(baseCol)) {
#   gg_color_hue6 <- hcl(h = seq(15, 375, length = 7), l=65, c=100)[1:6]
#   baseCol <- c(gg_color_hue6[3], gg_color_hue6[5], gg_color_hue6[2], gg_color_hue6[1], gg_color_hue6[4], gg_color_hue6[6]);
# }
# 
# numBases <- 5;
# 
# centerBase <- (1 + numBases) / 2;
# 
# v1 <- vF[1,1:6];
# V2 <- vF[2:(numBases),1:4];
# A <- matrix(0, numBases, 4);
# B <- matrix(0, 4, 4);
# 
# if (trDir == TRUE) {
#   v3 <- vF[(numBases + 1),1:2];
# }
# 
# for (l in 1:numBases) {
#   if (l < centerBase) {
#     A[l, ] <- V2[l, ];
#   } else if (l > centerBase) {
#     A[l, ] <- V2[l - 1, ];
#   }
# }
# A[centerBase,2] <- sum(v1[1:3]);
# A[centerBase,4] <- sum(v1[4:6]);
# 
# B[2, c(1, 3, 4)] <- v1[1:3] / sum(v1[1:3]);
# B[4, c(1, 2, 3)] <- v1[4:6] / sum(v1[4:6]);
# 
# 
# x_start <- c()
# x_end <- c()
# y_start <- c()
# y_end <- c()
# text_x <- c()
# text_y <- c()
# text_lab <- c()
# text_col <- c()
# nucBases <- c()
# num2base <- c("A", "C", "G", "T");
# 
# tempStartX <- 0;
# for (i in 1:numBases) {
#   x_start <- c(x_start, tempStartX + c(0, cumsum(A[i,1:3])))
#   x_end <- c(x_end, tempStartX + cumsum(A[i,1:4]));
#   y_start <- c(y_start, rep(0, 4))
#   y_end <- c(y_end, rep(1, 4))
#   nucBases <- c(nucBases, c("A", "C", "G", "T")) 
#   for (j in 1:4) {
#     tempPos <- c(0, cumsum(A[i,1:4]))
#     if (A[i,j] > 0.25) {
#       text_x <- c(text_x, tempStartX + 0.5 * (tempPos[j] + tempPos[j + 1]))
#       text_y <- c(text_y, 0.5 * (0 + 1))
#       text_lab <- c(text_lab, num2base[j])
#       text_col <- c(text_col, "w")
#     }
#   }
#   tempStartX <- tempStartX + 1.25
# }
# 
# tempStartX <- (centerBase - 1) * 1.25;
# x_start <- c(x_start, rep(tempStartX, 4))
# x_end <- c(x_end, rep(tempStartX + A[centerBase, 2], 4))
# y_start <- c(y_start, 2 + c(0, cumsum(B[2,1:3])))
# y_end <- c(y_end, 2 + cumsum(B[2,1:4]))
# nucBases <- c(nucBases, c("A", "C", "G", "T"))
# 
# tempPos <- c(0, cumsum(B[2,1:4]))
# for (j in 1:4) {
#   if (B[2,j] > 0.25) {
#     text_x <- c(text_x, tempStartX + 0.5 * A[centerBase, 2])
#     text_y <- c(text_y, 2 + 0.5 * (tempPos[j] + tempPos[j + 1]))
#     text_lab <- c(text_lab, num2base[j])
#     text_col <- c(text_col, "w")
#   }
# }
# 
# 
# tempStartX <- tempStartX + A[centerBase, 2]
# x_start <- c(x_start, rep(tempStartX, 4))
# x_end <- c(x_end, rep(tempStartX + A[centerBase, 4], 4))
# y_start <- c(y_start, 2 + c(0, cumsum(B[4,1:3])))
# y_end <- c(y_end, 2 + cumsum(B[4,1:4]))
# nucBases <- c(nucBases, c("A", "C", "G", "T"))
# 
# tempPos <- c(0, cumsum(B[4,1:4]))
# for (j in 1:4) {
#   if (B[4,j] > 0.25) {
#     text_x <- c(text_x, tempStartX + 0.5 * A[centerBase, 4])
#     text_y <- c(text_y, 2 + 0.5 * (tempPos[j] + tempPos[j + 1]))
#     text_lab <- c(text_lab, num2base[j])
#     text_col <- c(text_col, "w")
#   }
# }
# 
# rect_data <- data.frame(x_start = x_start, x_end = x_end, y_start = y_start, y_end = y_end, nucBases = nucBases)
# 
# text_data <- data.frame(x = text_x, y = text_y, label = text_lab, text_col = text_col)
# ##########
# # draw arrow
# xs <- c(1 /3, 2 / 3, 2 / 3, 5 / 6, 1 / 2, 1 / 6, 1 / 3, 1 / 3) + (centerBase - 1) * 1.25
# ys <- c(1 / 4, 1 / 4, 1 / 2, 1 / 2, 3 / 4, 1 / 2, 1 / 2, 1 / 4) + 1
# vs <- rep("arrow", length(xs))
# 
# arrow_poly <- data.frame(x = xs, y = ys, v = vs)
# 
# ggplot() + 
# geom_rect(data = df, aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end, fill = nucBases)) + 
# geom_polygon(data = arrow_poly, aes(x = x, y = y, fill = v))  +
# geom_text(data = text_data, aes(label = label, x = x, y = y, colour = text_col), size = 5) + 
# scale_colour_manual(values = c("#FFFFFF")) + 
# scale_fill_manual(values = c("A" = baseCol[1], "C" = baseCol[2], "G" = baseCol[3],
#                             "T" = baseCol[4], "arrow" = "#A8A8A8")) + 
# guides(fill=FALSE) +
# guides(colour=FALSE) +
# guides(size=FALSE) +
# theme(axis.text = element_blank(),
#       axis.ticks = element_blank(),
#       panel.background = element_blank(),
#       panel.grid = element_blank(), 
#       axis.title = element_blank())
#   
#   
# 
# 
# gg2 <- geom_polygon(data = arrow_poly, aes(x = x, y = y, fill = v)) + scale_fill_manual(values = c("#A8A8A8"))
# 
# 
# 
# 
# 
# xs <- c()
# ys <- c()
# types <- c()
# tempStartX <- 0
# for (i in 1:numBases) {
#   xs <- c(xs, rep(c(tempStartX, tempStartX + cumsum(A[i,1:4])), times = c(2, 4, 4, 4, 2)))
#   ys <- c(ys, rep(c(0, 1, 1, 0), 4))
#   types <- c(types, rep(c("A", "C", "G", "T"), each = 4))
#   tempStartX <- tempStartX + 1.25
# }
# 
# tempStartX <- (centerBase - 1) * 1.25;
# xs <- c(xs, rep(tempStartX + c(0, A[centerBase, 2], A[centerBase, 2], 0), 3))
# ys <- c(ys, rep(c(2, 2 + cumsum(B[2, c(1, 3, 4)])), times = c(2, 4, 4, 2)))
# types <- c(types, rep(c("A_2", "G_2", "T_2"), each = 4))
# 
# tempStartX <- tempStartX + A[centerBase, 2]
# xs <- c(xs, rep(tempStartX + c(0, A[centerBase, 4], A[centerBase, 4], 0), 3))
# ys <- c(ys, rep(c(2, 2 + cumsum(B[4, 1:3])), times = c(2, 4, 4, 2)))
# types <- c(types, rep(c("A", "C", "G"), each = 4))
# 
# 
# 
# poly_data <- data.frame(x = xs, y = ys, type = types)
# 
# 
# ggplot() + geom_polygon(data = poly_data, aes(x = x, y = y, group = type, fill = type))
# 
# 
# x_start <- c()
# x_end <- c()
# y_start <- c()
# y_end <- c()
# nucBases <- c()










