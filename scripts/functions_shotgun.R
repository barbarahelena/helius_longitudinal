## Functions shotgun
plot_feature_importance_shotgun <- function(path_true, top_n){
    theme_Publication <- function(base_size=12, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5,
                                              family = 'Helvetica'),
                    text = element_text(family = 'Helvetica'),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(), 
                    axis.line.x = element_line(colour="black"),
                    axis.ticks.x = element_line(),
                    axis.ticks.y = element_blank(),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "bottom",
                    # legend.direction = "horizontal",
                    legend.key.size= unit(0.2, "cm"),
                    legend.spacing  = unit(0, "cm"),
                    # legend.title = element_text(face="italic"),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = element_text(face="bold")
            ))
    } 
    cols <- list(low = '#ECE7F2',
                 mid = '#0570B0',
                 high = '#034E7B')
    r <- rio::import(file.path(path_true, 'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    shot <- readRDS("data/shotgun_taxtable.RDS")
    r <- r[1:top_n, ]
    r <- r %>% mutate(FeatName = factor(make.unique(FeatName), levels = rev(make.unique(FeatName))))
    mp <- mean(r$RelFeatImp)
    pl <- ggplot(data=r, aes(y=RelFeatImp, x=FeatName, fill=RelFeatImp)) + 
        theme_Publication()+
        scale_fill_gradient2(low=cols$low, mid = cols$mid, high=cols$high, space='Lab', name="",
                             midpoint = 50, guide = "none")+
        geom_bar(stat="identity")+
        coord_flip() +
        ylab("Relative Importance (%)")+
        xlab("") +
        theme(axis.text.x = element_text(size=12)) + 
        theme(axis.text.y = element_text(size=10))+
        theme(legend.key.size= unit(0.5, "cm"))+
        theme(legend.position = 'right')
    # ggsave(path = path_true, filename = 'plot_Feature_Importance.pdf', device = 'pdf', width = 10, height = 7)
}

plot_feature_importance_color_shotgun <- function(path_true, top_n){
    theme_Publication <- function(base_size=12, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5,
                                              family = 'Helvetica'),
                    text = element_text(family = 'Helvetica'),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(), 
                    axis.line.x = element_line(colour="black"),
                    axis.ticks.x = element_line(),
                    axis.ticks.y = element_blank(),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "bottom",
                    # legend.direction = "horizontal",
                    legend.key.size= unit(0.2, "cm"),
                    legend.spacing  = unit(0, "cm"),
                    # legend.title = element_text(face="italic"),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = element_text(face="bold")
            ))
    } 
    r <- rio::import(file.path(path_true, 'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    shot <- readRDS("data/shotgun_taxtable.RDS")
    r$Fam <- shot$Family[match(r$FeatName, shot$Species)]
    r <- r[1:top_n, ]
    r <- r %>% mutate(FeatName = factor(make.unique(FeatName), 
                                        levels = rev(make.unique(FeatName)))) %>% 
        mutate(Fam = as.factor(Fam))
    mp <- mean(r$RelFeatImp)
    colpal <- c("#00468BFF", "#ED0000FF", "#0099B4FF","gold",
                "lightblue1","darkgreen","#AD002AFF","mediumpurple1",
                "#42B540FF", "maroon3", "darkorchid4", "darkorange1",
                "grey40", "brown", "lightgrey", "black")
    colfam <- setNames(colpal[1:length(levels(r$Fam))], levels(r$Fam))
    pl <- ggplot(data=r, aes(y=RelFeatImp, x=FeatName, fill=Fam)) + 
        theme_Publication()+
        scale_fill_manual(values = colfam)+
        geom_bar(stat="identity")+
        coord_flip() +
        ylab("Relative Importance (%)")+
        xlab("") +
        theme(axis.text.x = element_text(size=12)) + 
        theme(axis.text.y = element_text(size=10))+
        theme(legend.key.size= unit(0.5, "cm"))+
        theme(legend.position = 'right')
    pl
    # ggsave(path = path_true, filename = 'plot_Feature_Importance.pdf', device = 'pdf', width = 10, height = 7)
}

plot_features_tests_shotgun <- function(input_path, output_path, top_n=10, labels=c("1", "0")){
    theme_Publication <- function(base_size=11, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5),
                    text = element_text(),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(), 
                    axis.line = element_line(colour="black"),
                    axis.ticks = element_line(),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "bottom",
                    # legend.direction = "horizontal",
                    legend.key.size= unit(0.2, "cm"),
                    legend.spacing  = unit(0, "cm"),
                    # legend.title = element_text(face="italic"),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = element_text(face="bold")
            ))
        
    } 
    plot_path <- file.path(output_path, 'plots')
    dir.create(plot_path)
    r <- rio::import(file.path(output_path,'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    input_data <- rio::import(file.path(input_path, 'X_data.txt'))
    feature_names <- read.csv(file.path(input_path, 'feat_ids.txt'), sep = '\t', header = F)
    shot <- readRDS("data/shotgun_taxtable.RDS")
    names(input_data) <- feature_names$V1
    if(top_n > ncol(input_data)){
        cat('\n\nRequested no. of features is higher than total number of features in model.\n
                 Showing all features in model.\n\n')
        top_n <- ncol(input_data)
    }
    features_tk <- r$FeatName[1:top_n]
    features_tk <- features_tk[! features_tk %in% c('random_variable1', 'random_variable2')]
    ft_tax <- shot$Species[match(features_tk, shot$rowname)]
    dd <- input_data %>% dplyr::select(any_of(features_tk))
    y <- rio::import(file.path(input_path, 'y_binary.txt'))
    dd$y <- y$V1
    dd$y <- factor(ifelse(dd$y==1, labels[1],labels[2]))
    comps <- list(c(labels[1],labels[2]))
    colorguide <- c(pal_simpsons()(2))
    print(colorguide)
    for(j in 1:length(features_tk)){
        asv <- features_tk[j]
        df <- dd %>% dplyr::select(all_of(asv), y)
        names(df)[1] <- 'Feature'
        tax_asv <- shot$Species[match(asv, shot$rowname)]
        pl <- ggplot(df, aes(x=y, y=Feature, fill=y+0.001))+
            geom_violin(trim = TRUE) +
            scale_fill_manual(values = colorguide, guide = FALSE)+
            geom_boxplot(width=0.1, fill="white")+
            theme_Publication()+
            theme(legend.position = 'none')+
            ggpubr::stat_compare_means(comparisons = comps, paired = F)+
            xlab('Group')+
            ylab('Relative abundance (%)')+
            scale_y_log10() + 
            ggtitle(tax_asv)
        cat(j, tax_asv, '\n')
        fname <- str_replace_all(tax_asv, "[*\";,:/\\\\ ]","_")
        ggsave(pl, path = plot_path, filename = paste0(j, '_',fname, '.pdf'), device = 'pdf', width = 5, height = 5)
    }
}

plot_features_top_shotgun <- function(input_path, output_path, top_n=20, nrow=4, labels){
    theme_Publication <- function(base_size=11, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5),
                    text = element_text(),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(), 
                    axis.line = element_line(colour="black"),
                    axis.ticks = element_line(),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "bottom",
                    # legend.direction = "horizontal",
                    legend.key.size= unit(0.2, "cm"),
                    legend.spacing  = unit(0, "cm"),
                    # legend.title = element_text(face="italic"),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = element_text(face="bold", size = rel(0.8))
            ))
        
    } 
    plot_path <- file.path(output_path, 'plots')
    dir.create(plot_path)
    r <- rio::import(file.path(output_path,'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    shot <- readRDS("data/shotgun_taxtable.RDS")
    input_data <- rio::import(file.path(input_path, 'X_data.txt'))
    feature_names <- read.csv(file.path(input_path, 'feat_ids.txt'), sep = '\t', header = F)
    names(input_data) <- feature_names$V1
    if(top_n > ncol(input_data)){
        cat('\n\nRequested no. of features is higher than total number of features in model.\nShowing all features in model.\n\n')
        top_n <- ncol(input_data)
    }
    features_tk <- r$FeatName[1:top_n]
    features_tk <- features_tk[! features_tk %in% c('random_variable1', 'random_variable2')]
    dd <- input_data %>% dplyr::select(any_of(features_tk))
    y <- rio::import(file.path(input_path, 'y_binary.txt'))
    dd$y <- y$V1
    dd$y <- factor(ifelse(dd$y==1, labels[1],labels[2]))
    df <- dd %>% pivot_longer(-y, names_to = 'features', values_to = 'values')
    df <- df %>% mutate(features = as.factor(features),
                        features = fct_inorder(features),
                        values = values
    )
    colorguide <- c(pal_simpsons()(2))
    comps <- list(c(labels[1], labels[2]))
    pl <- ggplot(df, aes(x=y, y=values+0.001))+
        geom_violin(aes(fill=y), trim = TRUE)+
        scale_fill_manual(values = colorguide, guide = FALSE) +
        geom_boxplot(width=0.1, fill="white")+
        theme_Publication()+
        theme(legend.position = 'none')+
        labs(x='Group', y = 'Relative abundance (%)')+
        ggpubr::stat_compare_means(comparisons = comps, paired = F, size = rel(3.0))+
        scale_y_log10() +
        facet_wrap(~ features, nrow=nrow, scales = 'free')
    pl
    # ggsave(pl, path = plot_path, filename = paste0('top_',top_n,'_features.pdf'), device = 'pdf', width=15, height = 14)
    # ggsave(pl, path = plot_path, filename = paste0('top_',top_n,'_features.svg'), device = 'svg', width=15, height = 14)
}


## Functions pathways

plot_feature_importance_pathways <- function(path_true, top_n){
    theme_Publication <- function(base_size=12, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5,
                                              family = 'Helvetica'),
                    text = element_text(family = 'Helvetica'),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(), 
                    axis.line.x = element_line(colour="black"),
                    axis.ticks.x = element_line(),
                    axis.ticks.y = element_blank(),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "bottom",
                    # legend.direction = "horizontal",
                    legend.key.size= unit(0.2, "cm"),
                    legend.spacing  = unit(0, "cm"),
                    # legend.title = element_text(face="italic"),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = element_text(face="bold")
            ))
    } 
    cols <- list(low = '#ECE7F2',
                 mid = '#0570B0',
                 high = '#034E7B')
    r <- rio::import(file.path(path_true, 'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    r <- r[1:top_n, ]
    mp <- mean(r$RelFeatImp)
    pl <- ggplot(data=r, aes(y=RelFeatImp, x=FeatName, fill=RelFeatImp)) + 
        theme_Publication()+
        scale_fill_gradient2(low=cols$low, mid = cols$mid, high=cols$high, space='Lab', name="",
                             midpoint = 50, guide = "none")+
        geom_bar(stat="identity")+
        coord_flip() +
        ylab("Relative Importance (%)")+
        xlab("") +
        ggtitle("Feature importance: best predicting pathways for sex")
    theme(axis.text.x = element_text(size=12)) + 
        theme(axis.text.y = element_text(size=10))+
        theme(legend.key.size= unit(0.5, "cm"))+
        theme(legend.position = 'right')
    pl
    # ggsave(path = path_true, filename = 'plot_Feature_Importance.pdf', device = 'pdf', width = 10, height = 7)
}

plot_features_tests_pathways <- function(input_path, output_path, top_n=10, labels=c("1", "0")){
    theme_Publication <- function(base_size=11, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5),
                    text = element_text(),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(), 
                    axis.line = element_line(colour="black"),
                    axis.ticks = element_line(),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "bottom",
                    # legend.direction = "horizontal",
                    legend.key.size= unit(0.2, "cm"),
                    legend.spacing  = unit(0, "cm"),
                    # legend.title = element_text(face="italic"),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = element_text(face="bold")
            ))
        
    } 
    plot_path <- file.path(output_path, 'plots')
    dir.create(plot_path)
    r <- rio::import(file.path(output_path,'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    input_data <- rio::import(file.path(input_path, 'X_data.txt'))
    feature_names <- read.csv(file.path(input_path, 'feat_ids.txt'), sep = '\t', header = F)
    names(input_data) <- feature_names$V1
    if(top_n > ncol(input_data)){
        cat('\n\nRequested no. of features is higher than total number of features in model.\n
                 Showing all features in model.\n\n')
        top_n <- ncol(input_data)
    }
    features_tk <- r$FeatName[1:top_n]
    features_tk <- features_tk[! features_tk %in% c('random_variable1', 'random_variable2')]
    dd <- input_data %>% dplyr::select(any_of(features_tk))
    y <- rio::import(file.path(input_path, 'y_binary.txt'))
    dd$y <- y$V1
    dd$y <- factor(ifelse(dd$y==1, labels[1],labels[2]))
    comps <- list(c(labels[1],labels[2]))
    colorguide <- c(pal_simpsons()(2))
    #print(colorguide)
    for(j in 1:length(features_tk)){
        asv <- features_tk[j]
        df <- dd %>% dplyr::select(all_of(asv), y)
        names(df)[1] <- 'Feature'
        tax_asv <- paths$label[match(asv, paths$key)]
        pl <- ggplot(df, aes(x=y, y=Feature, fill=y))+
            geom_violin(trim = TRUE) +
            scale_fill_manual(values = colorguide, guide = FALSE)+
            geom_boxplot(width=0.1, fill="white")+
            theme_Publication()+
            theme(legend.position = 'none')+
            ggpubr::stat_compare_means(comparisons = comps, paired = F, method = "wilcox.test",
                                       tip.length = 0)+
            xlab('Group')+
            ylab('Relative abundance (cpm)')+
            ggtitle(tax_asv)
        fname <- tax_asv
        cat(j, fname, '\n')
        fname <- str_replace_all(fname, "[*\";,:/\\\\ ]","_")
        #print(pl)
        ggsave(pl, path = plot_path, filename = paste0(j, '_',fname, '.pdf'), device = 'pdf', width = 5, height = 5)
    }
}

plot_features_top_pathways <- function(input_path, output_path, top_n=20, nrow=4, labels){
    theme_Publication <- function(base_size=11, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5),
                    text = element_text(),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(), 
                    axis.line = element_line(colour="black"),
                    axis.ticks = element_line(),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "bottom",
                    # legend.direction = "horizontal",
                    legend.key.size= unit(0.2, "cm"),
                    legend.spacing  = unit(0, "cm"),
                    # legend.title = element_text(face="italic"),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = element_text(face="bold", size = rel(0.5))
            ))
        
    } 
    plot_path <- file.path(output_path, 'plots')
    dir.create(plot_path)
    r <- rio::import(file.path(output_path,'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    input_data <- rio::import(file.path(input_path, 'X_data.txt'))
    feature_names <- read.csv(file.path(input_path, 'feat_ids.txt'), sep = '\t', header = F)
    names(input_data) <- feature_names$V1
    if(top_n > ncol(input_data)){
        cat('\n\nRequested no. of features is higher than total number of features in model.\nShowing all features in model.\n\n')
        top_n <- ncol(input_data)
    }
    features_tk <- r$FeatName[1:top_n]
    features_tk <- features_tk[! features_tk %in% c('random_variable1', 'random_variable2')]
    dd <- input_data %>% dplyr::select(any_of(features_tk))
    y <- rio::import(file.path(input_path, 'y_binary.txt'))
    dd$y <- y$V1
    dd$y <- factor(ifelse(dd$y==1, labels[1],labels[2]))
    df <- dd %>% pivot_longer(-y, names_to = 'features', values_to = 'values')
    df <- df %>% mutate(features = as.factor(features),
                        features = fct_inorder(features),
                        values = values
    )
    colorguide <- c(pal_simpsons()(2))
    comps <- list(c(labels[1], labels[2]))
    pl <- ggplot(df, aes(x=y, y=values))+
        geom_violin(aes(fill=y), trim = TRUE)+
        scale_fill_manual(values = colorguide, guide = FALSE) +
        geom_boxplot(width=0.1, fill="white")+
        theme_Publication()+
        theme(legend.position = 'none')+
        labs(x='Group', y = 'Relative abundance (cpm)')+
        ggpubr::stat_compare_means(comparisons = comps, paired = F, size = rel(3.0), 
                                   method = "wilcox.test", tip.length = 0)+
        facet_wrap(~ features, nrow=nrow, scales = 'free')
    pl
    # ggsave(pl, path = plot_path, filename = paste0('top_',top_n,'_features.pdf'), device = 'pdf', width=15, height = 14)
    # ggsave(pl, path = plot_path, filename = paste0('top_',top_n,'_features.svg'), device = 'svg', width=15, height = 14)
}


