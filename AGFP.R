########################################################
# Giving Info about the files
sep.char <- ''

var.vec <- c('pBS-Mul1', 
             'p-A-GFP', 
             'p-A-GFP_PXXV', 
             'p-A-GFP_iSD', 
             'p-A-GFP_iSDaa')

# Define which column contains the blank values
blank.var <- 'pBS-Mul1'
wt.var <- 'p-A-GFP'

# Define the GFP signal files and their timepoints
gfp.files <- c('GFP_6h.txt', 'GFP_8h.txt', 'GFP_24h.txt')
gfp.f.desc <- c('6h', '8h', '24h')

# As for GFP but the activity files
act.files <- c('Act_6ha_tab.txt', 'Act_8ha_tab.txt',
               'Act_24h_Calc50.txt')
act.f.desc <- c('6h', '8h', '24h')

# As for GFP but the OD measuments
od.files <- c('OD_6h.txt', 'OD_8h.txt', 'OD_24h.txt')
od.f.desc <- c('6h', '8h', '24h')

########################################################
# Load the needed packages
require(reshape)
require(ggplot2)

# Create a list with the measurement files
# and their descriptions
files <- list(gfp.files, act.files, od.files)
files.desc <- list(gfp.f.desc, act.f.desc, od.f.desc)

names(files) <- c('GFP', 'Act', 'OD')
names(files.desc) <- names(files)


# Reading in the files

df.input <- vector(mode = 'list', length = 3)
names(df.input) <- names(files)


# Go through each group
for(sel.name in names(df.input)){
  
  sel.files <- files[[sel.name]]
  sel.f.desc <- files.desc[[sel.name]]
  
  # Go through each file
  
  tmp.files <- vector(mode = 'list', length = length(sel.files))
  names(tmp.files) <- sel.f.desc
  tmp.out <- NULL
  
  i <- 1
  for(f in sel.files){
    
    tmp <- read.table(sel.files[i], sep=sep.char, dec = ",")  # Read data
    colnames(tmp) <- var.vec  # Get column names
    
    # create a long table
    tmp <- melt(tmp)
    tmp <- data.frame(tmp, 
                      'Cond' = sel.name, 
                      'Time' = sel.f.desc[i],
                      'Bio.Rep' = rep(c('A', 'A', 'B', 'B', 'C', 'C', 'D', 'D'), length(tmp[,1])/8)
                      )
    colnames(tmp)[1] <- 'Plasmid' 
    
    tmp.files[[i]] <- tmp
    tmp.out <- rbind(tmp.out, tmp)
    i = i+1
  }
  
  # save tmp.files in df.input
  df.input[[sel.name]] <- tmp.out
  
}


processed <- NULL
# Go through the activity and GFP files
for(f in c('Act', 'GFP')){
  tmp <- df.input[[f]]
  print(tmp)
  tmp$value[tmp$value < 0 ] <- NA  # Remove values below 0 and set to NA
  print(tmp)
  
  tmp$value <- tmp$value / df.input$OD$value  # OD normalization
  
  blank <- tmp[tmp$Plasmid == blank.var, ]
  
  # Calculate the blank value from the mean of technical replicates
  blank <- condense(blank, c('Plasmid', 'Cond', 'Time', 'Bio.Rep'), mean)
  blank$result <- as.numeric(blank$result)
  
  
  # Substract from each measurement the mean blank
  tmp.out <- NULL
  for(lev in levels(blank$Time)){
    tmp.h <- tmp[tmp$Time == lev, ]
    blank.h <- blank[blank$Time == lev, ]
    blank.h <- mean(blank.h$result)
    
    if(is.na(blank.h)){
      tmp.h$value <- tmp.h$value - 0
    }else{
      tmp.h$value <- tmp.h$value - blank.h
    }
    
    tmp.out <- rbind(tmp.out, tmp.h)
  }
  
  processed <- rbind(processed, tmp.out)
}

# Adjust to activity in 1 ml (here, 100ul was used in the assay)
processed[processed$Cond=='Act',]$value <- processed[processed$Cond=='Act',]$value *10

# Write the resulting long table to a file
write.table(processed,
            file = '_FullBioRepTable.txt',
            sep = '\t',
            quote = F,
            row.names = F)

# Calculate the biological replicates using
# the technical replicates' mean value
means.condense <- condense(processed, c('Plasmid', 'Cond', 'Time', 'Bio.Rep'), mean)
colnames(means.condense)[5] <- 'Mean' 
means.condense$Mean <- as.numeric(means.condense$Mean)

# Calculate the technical
sd.condense <- condense(processed, c('Plasmid', 'Cond', 'Time', 'Bio.Rep'), sd)
colnames(sd.condense)[5] <- 'SD' 
sd.condense$SD <- as.numeric(sd.condense$SD)

df.summary <- data.frame(means.condense,
                         'SD' = sd.condense$SD)

write.table(df.summary,
            file = '_SummaryTable.txt',
            sep = '\t',
            quote = F,
            row.names = F)

cpalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
pdf('_ResultPlot.pdf', width=7, height = 6)
ggplot(df.summary, aes(x=Plasmid, y=Mean, group = Bio.Rep, color=Plasmid)) + 
  facet_grid(Cond ~ Time, scales = "free_y") +
  geom_pointrange(aes(ymin=Mean-SD, ymax=Mean+SD), position=position_dodge(width=0.5)) +
  scale_color_manual(values=cpalette) +
  ylab('Fluorescence [485/520nm]      Converted Substrate [A(620nm)]') +
  xlab('Construct') +
  theme(
    axis.text.x = element_text(angle=90,
                               hjust = 1,
                               vjust = 0.5),
    legend.position="none"
  ) +
  xlab('')
dev.off()