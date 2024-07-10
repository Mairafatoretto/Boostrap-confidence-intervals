#################################################################
##########Figure 4 and  5########################################
################################################################

############Figure 4#############################################
# Install packages if not already installed
if (!requireNamespace("ggplot2", quietly = TRUE)) {
   install.packages("ggplot2")
}
if (!requireNamespace("ggpubr", quietly = TRUE)) {
   install.packages("ggpubr")
}
if (!requireNamespace("extrafont", quietly = TRUE)) {
   install.packages("extrafont")
}

# Load necessary packages
library(ggplot2)
library(ggpubr)
library(extrafont)

# Define variables and settings
tam_leg <- 8
tam_val <- 6
tam_tit <- 12

# Import available fonts from the system
font_import()

# Load all fonts
loadfonts()

# Set a default font for all plots
fonte <- "ArialMT"  # Choose a font available on your system

# Set the font in ggplot theme
theme_set(theme_minimal(base_family = fonte))

# Read data
coverage_parameters <- read.table('coverage_parameters.txt', header = TRUE, dec = '.', sep = '\t')
coverage_components <- read.table('coverage_components.txt', header = TRUE, dec = '.', sep = '\t')

# Sort data as needed
ord <- c("BL1", "BL2", "BL3", "BL4", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9",
         "F1.0", "F1.1", "F1.2", "F1.3", "F1.4")

coverage_parameters$parameter <- factor(coverage_parameters$parameter, levels = ord)

# Divide data for different plots
coverage_parameters1 <- coverage_parameters[1:72, ]
coverage_parameters2 <- coverage_parameters[75:144, ]
coverage_components1 <- coverage_components[1:16, ]
coverage_components2 <- coverage_components[17:32, ]

# Define theme adjustments with the chosen font
theme_custom <- list(
   panel.background = element_rect(fill = "white"),
   plot.title = element_text(hjust = 0.8, size = tam_tit, vjust = 3),
   text = element_text(family = fonte),
   panel.grid.major = element_line(color = "lightgrey", size = 0.5),
   panel.grid.minor = element_blank(),
   axis.title.x = element_text(size = 25),
   axis.text.x = element_text(angle = 320, size = 25),
   axis.title.y = element_text(size = 25),
   axis.text.y = element_text(size = 20),
   legend.key.height = unit(0.8, "line")
)


# PLOT 1 (p1)
p1 <- ggplot(coverage_parameters1, aes(x = parameter, y = value, group = interaction(estimation, interval))) +
   geom_line(aes(color = estimation, linetype = interval), size = 0.8) +
   geom_point(aes(color = estimation, linetype = interval)) +
   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 0.95, by = 0.3)) +
   geom_hline(yintercept = 0.95, linetype = "longdash", size = 0.5, color = "grey") +
   xlab("Parameters") +
   scale_color_manual(values = c('red', 'blue3'), labels = c("Method 2", "Method 3"), name = "Bootstrap Method   ") +
   scale_linetype_manual(values = c("solid", "dashed"), labels = c("CI", "HDI"), name = "Interval") +
   scale_x_discrete(labels = c(expression(gamma[1]), expression(gamma[2]), expression(gamma[3]), expression(gamma[4]),
                               expression(beta[0.2]), expression(beta[0.3]), expression(beta[0.4]), expression(beta[0.5]),
                               expression(beta[0.6]), expression(beta[0.7]), expression(beta[0.8]), expression(beta[0.9]),
                               expression(beta[0.10]), expression(beta[0.11]), expression(beta[0.12]), expression(beta[0.13]),
                               expression(beta[0.14]), expression(beta[1]))) +
   theme_minimal() +ylab("Coverage rate") +
   xlab("Parameters") +
   theme(plot.margin = margin(15, 10, 10, 10),
         legend.position = "none",
         panel.background = element_rect(fill = "white", colour = "black"),
         plot.title = element_text(hjust = 0.8, size = tam_tit, vjust = 3),
         text = element_text(family = fonte),
         panel.grid.major = element_line(color = "lightgrey", size = 0.5),
         panel.grid.minor = element_blank(),
         axis.title.x = element_text(size = 12),
         axis.text.x = element_text(angle = 320, size = 10),
         axis.title.y = element_text(size = 12),
         axis.text.y = element_text(size = 10),
         legend.key.height = unit(0.8, "line")
   )

# PLOT 2 (p2)
p2 <- ggplot(coverage_parameters2, aes(x = parameter, y = value, group = interaction(estimation, interval))) +
   geom_line(aes(color = estimation, linetype = interval), size = 0.8) +
   geom_point(aes(color = estimation, linetype = interval)) +ylab("Coverage rate") +
   xlab("Parameters") +
   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 0.95, by = 0.3)) +
   geom_hline(yintercept = 0.95, linetype = "longdash", size = 0.5, color = "grey") +
   xlab("Parameters") +
   scale_color_manual(values = c('red', 'blue3'), labels = c("Method 2", "Method 3"), name = "Bootstrap Method   ") +
   scale_linetype_manual(values = c("solid", "dashed"), labels = c("CI", "HDI"), name = "Interval") +
   scale_x_discrete(labels = c(expression(gamma[1]), expression(gamma[2]), expression(gamma[3]), expression(gamma[4]),
                               expression(beta[0.2]), expression(beta[0.3]), expression(beta[0.4]), expression(beta[0.5]),
                               expression(beta[0.6]), expression(beta[0.7]), expression(beta[0.8]), expression(beta[0.9]),
                               expression(beta[0.10]), expression(beta[0.11]), expression(beta[0.12]), expression(beta[0.13]),
                               expression(beta[0.14]), expression(beta[1]))) +
   theme_minimal() +
   theme(plot.margin = margin(15, 10, 10, 10),
         legend.position = "none",
         panel.background = element_rect(fill = "white", colour = "black"),
         plot.title = element_text(hjust = 0.8, size = tam_tit, vjust = 3),
         text = element_text(family = fonte),
         panel.grid.major = element_line(color = "lightgrey", size = 0.5),
         panel.grid.minor = element_blank(),
         axis.title.x = element_text(size = 12),
         axis.text.x = element_text(angle = 320, size = 10),
         axis.title.y = element_text(size = 12),
         axis.text.y = element_text(size = 10),
         legend.key.height = unit(0.8, "line")
   )

# PLOT 3 (p3)
p3 <- ggplot(coverage_components1, aes(x = component, y = value, group = interaction(estimation, interval))) +
   geom_line(aes(color = estimation, linetype = interval), size = 0.8) +
   geom_point(aes(color = estimation, linetype = interval)) +
   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 0.95, by = 0.3)) +
   geom_hline(yintercept = 0.95, linetype = "longdash", size = 0.5, color = "grey") +
   ylab("Coverage rate") +
   xlab("Parameters") +
   scale_color_manual(values = c('red', 'blue3'), labels = c("Method 2", "Method 3"), name = "Bootstrap Method   ") +
   scale_linetype_manual(values = c("solid", "dashed"), labels = c("CI", "HDI"), name = "Interval") +
   scale_x_discrete(labels = c(expression(sigma[O]), expression(sigma[I]), expression(sigma[S]), expression(sigma[IS]))) +
   theme_minimal() +
   theme(plot.margin = margin(20, 10, 20, 10),
         legend.position = "none",
         panel.background = element_rect(fill = "white", colour = "black"),
         plot.title = element_text(hjust = 0.8, size = tam_tit, vjust = 3),
         text = element_text(family = fonte),
         panel.grid.major = element_line(color = "lightgrey", size = 0.5),
         panel.grid.minor = element_blank(),
         axis.title.x = element_text(size = 12),
         axis.text.x = element_text(size = 10),
         axis.title.y = element_text(size = 12),
         axis.text.y = element_text(size = 10),
         legend.key.height = unit(0.8, "line")
   )

# PLOT 4 (p4)
p4 <- ggplot(coverage_components2, aes(x = component, y = value, group = interaction(estimation, interval))) +
   geom_line(aes(color = estimation, linetype = interval), size = 0.8) +
   geom_point(aes(color = estimation, linetype = interval)) +
   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 0.95, by = 0.3)) +
   geom_hline(yintercept = 0.95, linetype = "longdash", size = 0.5, color = "grey")+
   ylab("Coverage rate") +
   xlab("Parameters") +
   scale_color_manual(values = c('red', 'blue3'), labels = c("Method 2", "Method 3"), name = "Bootstrap Method   ") +
   scale_linetype_manual(values = c("solid", "dashed"), labels = c("CI", "HDI"), name = "Interval") +
   scale_x_discrete(labels = c(expression(sigma[O]), expression(sigma[I]), expression(sigma[S]), expression(sigma[IS]))) +
   theme_minimal() +
   theme(plot.margin = margin(20, 10, 20, 10),
         legend.position = "none",
         panel.background = element_rect(fill = "white", colour = "black"),
         plot.title = element_text(hjust = 0.8, size = tam_tit, vjust = 3),
         text = element_text(family = fonte),
         panel.grid.major = element_line(color = "lightgrey", size = 0.5),
         panel.grid.minor = element_blank(),
         axis.title.x = element_text(size = 12),
         axis.text.x = element_text(size = 10),
         axis.title.y = element_text(size = 12),
         axis.text.y = element_text(size = 10),
         legend.key.height = unit(0.8, "line")
   )


# Combine plots using ggarrange
combined_plot <- ggarrange(
   p1, p2, p3, p4,
   ncol = 2,  # Number of columns
   nrow = 2,  # Number of rows
   common.legend = TRUE,  # Add a common legend
   legend = "right",  # Position of the common legend
   labels = c("(a)", "(b)", "(c)", "(d)"), # Labels for subplots
   hjust = (-7.2),
   vjust = c(1.1, 1.1, 1.5, 1.5),  # Vertical justification adjustments
   align = c("none", "h", "v", "hv"),  # Alignment of plots
   font.label = list(size = 12, color = "black", face = "plain")  # Label font settings
)

print(combined_plot)


########FIGURE 5################################################

require(tidyr)

# Carregue os pacotes necessários
library(ggplot2)
library(ggpubr)


total0 <- load("ML_1000_M3.RData")
total0 <- coverage[1:100]

total1 <- load("REML_1000_M3.RData")
total1 <- coverage[1:100]

#calculating width of percentile confidence intervals
REML <- lapply(1:100, function(i){  
   list(
      
      "NewID1" <-  total1[[i]]$qtl[,2][1]-total1[[i]]$qtl[,1][1],
      "NewID2"  <-    total1[[i]]$qtl[,2][2]-total1[[i]]$qtl[,1][2],
      "NewID3"  <-    total1[[i]]$qtl[,2][3]-total1[[i]]$qtl[,1][3],
      "NewID4"  <-   total1[[i]]$qtl[,2][4]-total1[[i]]$qtl[,1][4],
      "FUNG2"  <-    total1[[i]]$qtl[,2][5]-total1[[i]]$qtl[,1][5],
      "FUNG3"  <-    total1[[i]]$qtl[,2][6]-total1[[i]]$qtl[,1][6],
      "FUNG4"  <-    total1[[i]]$qtl[,2][7]-total1[[i]]$qtl[,1][7],
      "FUNG5"  <-    total1[[i]]$qtl[,2][8]-total1[[i]]$qtl[,1][8],
      "FUNG6"  <-    total1[[i]]$qtl[,2][9]-total1[[i]]$qtl[,1][9],
      "FUNG7"  <-    total1[[i]]$qtl[,2][10]-total1[[i]]$qtl[,1][10],
      "FUNG8"  <-    total1[[i]]$qtl[,2][11]-total1[[i]]$qtl[,1][11],
      "FUNG9"  <-    total1[[i]]$qtl[,2][12]-total1[[i]]$qtl[,1][12],
      "FUNG10"  <-    total1[[i]]$qtl[,2][13]-total1[[i]]$qtl[,1][13],
      "FUNG11"  <-    total1[[i]]$qtl[,2][14]-total1[[i]]$qtl[,1][14],
      "FUNG12"  <-    total1[[i]]$qtl[,2][15]-total1[[i]]$qtl[,1][15],
      "FUNG13"  <-    total1[[i]]$qtl[,2][16]-total1[[i]]$qtl[,1][16],
      "FUNG14"  <-    total1[[i]]$qtl[,2][17]-total1[[i]]$qtl[,1][7],
      "DOSE"  <-   total1[[i]]$qtl[,2][18]-total1[[i]]$qtl[,1][18],
      "obs"  <-    total1[[i]]$qtl[,2][19]-total1[[i]]$qtl[,1][19],
      "inter"  <-    total1[[i]]$qtl[,2][20]-total1[[i]]$qtl[,1][20],
      "slope"  <-    total1[[i]]$qtl[,2][21]-total1[[i]]$qtl[,1][21],
      "corr"  <-   total1[[i]]$qtl[,2][22]-total1[[i]]$qtl[,1][22]
   ) 
}) 




REML <- lapply(as.data.frame(do.call("rbind", REML)), as.numeric)
REML <- data.frame(matrix(unlist(REML),nrow=100))
REML$method <- "REML"


ML <- lapply(1:100, function(i){  
   list(
      
      "NewID1" <-  total0[[i]]$qtl[,2][1]-total0[[i]]$qtl[,1][1],
      "NewID2"  <-    total0[[i]]$qtl[,2][2]-total0[[i]]$qtl[,1][2],
      "NewID3"  <-    total0[[i]]$qtl[,2][3]-total0[[i]]$qtl[,1][3],
      "NewID4"  <-   total0[[i]]$qtl[,2][4]-total0[[i]]$qtl[,1][4],
      "FUNG2"  <-    total0[[i]]$qtl[,2][5]-total0[[i]]$qtl[,1][5],
      "FUNG3"  <-    total0[[i]]$qtl[,2][6]-total0[[i]]$qtl[,1][6],
      "FUNG4"  <-    total0[[i]]$qtl[,2][7]-total0[[i]]$qtl[,1][7],
      "FUNG5"  <-    total0[[i]]$qtl[,2][8]-total0[[i]]$qtl[,1][8],
      "FUNG6"  <-    total0[[i]]$qtl[,2][9]-total0[[i]]$qtl[,1][9],
      "FUNG7"  <-    total0[[i]]$qtl[,2][10]-total0[[i]]$qtl[,1][10],
      "FUNG8"  <-    total0[[i]]$qtl[,2][11]-total0[[i]]$qtl[,1][11],
      "FUNG9"  <-    total0[[i]]$qtl[,2][12]-total0[[i]]$qtl[,1][12],
      "FUNG10"  <-    total0[[i]]$qtl[,2][13]-total0[[i]]$qtl[,1][13],
      "FUNG11"  <-    total0[[i]]$qtl[,2][14]-total0[[i]]$qtl[,1][14],
      "FUNG12"  <-    total0[[i]]$qtl[,2][15]-total0[[i]]$qtl[,1][15],
      "FUNG13"  <-    total0[[i]]$qtl[,2][16]-total0[[i]]$qtl[,1][16],
      "FUNG14"  <-    total0[[i]]$qtl[,2][17]-total0[[i]]$qtl[,1][7],
      "DOSE"  <-   total0[[i]]$qtl[,2][18]-total0[[i]]$qtl[,1][18],
      "obs"  <-    total0[[i]]$qtl[,2][19]-total0[[i]]$qtl[,1][19],
      "inter"  <-    total0[[i]]$qtl[,2][20]-total0[[i]]$qtl[,1][20],
      "slope"  <-    total0[[i]]$qtl[,2][21]-total0[[i]]$qtl[,1][21],
      "corr"  <-   total0[[i]]$qtl[,2][22]-total0[[i]]$qtl[,1][22]
   ) 
}) 


ML <- lapply(as.data.frame(do.call("rbind", ML)), as.numeric)
ML <- data.frame(matrix(unlist(ML),nrow=100))
ML$method <- "ML"


longer_data <- rbind(ML,REML) %>%
   pivot_longer(X1:X22, names_to = "parameters", values_to = "values")




# Ajuste para garantir que `longer_data` esteja como um tibble
library(tibble)
longer_data <- as_tibble(longer_data)

# Definição de parâmetros para o tamanho do texto e da fonte
tam_tit <- 15
tam_leg <- 12
tam_val <- 9
fonte <- 'Times'

# Definir a ordem das variáveis
ord <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9",
         "X10","X11","X12","X13","X14","X15","X16","X17",
         "X18","X19","X20","X21","X22")
longer_data$parameters <- factor(longer_data$parameters, levels = ord)

# Criar os gráficos individuais
# Gráfico 5a - Width of confidence intervals
plot_5a <- ggplot(longer_data, aes(x = parameters, y = values)) +
   geom_boxplot(aes(fill = factor(method)), alpha = 0.8) +
   ylab("Width of confidence intervals\n") +
   xlab("\nParameters") +
   scale_fill_manual(values = c('deeppink', 'cornflowerblue'),
                     labels = c("ML ", "REML"),
                     name = "Estimation \n method") +
   scale_x_discrete(labels = c(expression(gamma[1]), expression(gamma[2]), expression(gamma[3]), expression(gamma[4]),
                               expression(beta[0.2]), expression(beta[0.3]), expression(beta[0.4]), expression(beta[0.5]),
                               expression(beta[0.6]), expression(beta[0.7]), expression(beta[0.8]), expression(beta[0.9]),
                               expression(beta[0.1[0]]), expression(beta[0.11]), expression(beta[0.12]), expression(beta[0.13]),
                               expression(beta[0.14]), expression(beta[1]), expression(sigma[O]), expression(sigma[I]),
                               expression(sigma[S]), expression(sigma[IS]))) +
   theme(plot.margin = margin(20, 10, 10, 10),,legend.key = element_rect(fill = "white", colour = "white"),
         legend.key.size = unit(1.5, "cm"),
         panel.background = element_rect(fill = "white", colour = "black"),
         plot.title = element_text(hjust = 0.5, size = tam_tit, face = "bold", vjust = 3),
         text = element_text(family = fonte), panel.grid = element_blank(),
         legend.title = element_text(size = tam_leg),
         legend.text = element_text(size = tam_val),
         axis.title.x = element_text(size = tam_leg),
         axis.text.x = element_text(angle = 320, size = tam_leg),
         axis.title.y = element_text(size = tam_leg),
         axis.text.y = element_text(size = tam_val),
         legend.key.height = unit(0.8, "line"))

# Gráfico 5b - Lower limits of confidence intervals
plot_5b <- ggplot(longer_data, aes(x = parameters, y = values)) +
   geom_boxplot(aes(fill = factor(method)), alpha = 0.8) +
   ylab("Lower limits of confidence intervals\n") +
   xlab("\nParameters") +
   scale_fill_manual(values = c('deeppink', 'cornflowerblue'),
                     labels = c("ML ", "REML"),
                     name = "Estimation \n method") +
   scale_x_discrete(labels = c(expression(gamma[1]), expression(gamma[2]), expression(gamma[3]), expression(gamma[4]),
                               expression(beta[0.2]), expression(beta[0.3]), expression(beta[0.4]), expression(beta[0.5]),
                               expression(beta[0.6]), expression(beta[0.7]), expression(beta[0.8]), expression(beta[0.9]),
                               expression(beta[0.1[0]]), expression(beta[0.11]), expression(beta[0.12]), expression(beta[0.13]),
                               expression(beta[0.14]), expression(beta[1]), expression(sigma[O]), expression(sigma[I]),
                               expression(sigma[S]), expression(sigma[IS]))) +
   theme(plot.margin = margin(20, 10, 10, 10),legend.key = element_rect(fill = "white", colour = "white"),
         legend.key.size = unit(1.5, "cm"),
         panel.background = element_rect(fill = "white", colour = "black"),
         plot.title = element_text(hjust = 0.5, size = tam_tit, face = "bold", vjust = 3),
         text = element_text(family = fonte), panel.grid = element_blank(),
         legend.title = element_text(size = tam_leg),
         legend.text = element_text(size = tam_val),
         axis.title.x = element_text(size = tam_leg),
         axis.text.x = element_text(angle = 320, size = tam_leg),
         axis.title.y = element_text(size = tam_leg),
         axis.text.y = element_text(size = tam_val),
         legend.key.height = unit(0.8, "line"))

# Combinar os gráficos usando ggarrange
combined_boxplot <- ggarrange(plot_5a, plot_5b, ncol = 2, common.legend = TRUE, legend = "right",labels = c("(a)", "(b)", "(c)", "(d)"), # Rótulos
                           hjust=(-8.3),
                           vjust=c(1.2,1.2),
                           font.label = list(size = 12, color = "black", face = "plain"))

# Exibir o gráfico combinado
print(combined_boxplot)
