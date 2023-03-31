# ---
# title: "Staphylococcus pseudintermedius AMR ananlysis"
# authors: "Caroline Mevis, Casey Cazer, Ritwik Sadhu; Yuchen Xu"
# ---


require(checkpoint)

checkpoint("2023-03-01")

require(tidyverse)
require(gt)
require(mgsub)
require(writexl)
require(ggpubr)
require(reshape2)
require(directlabels)
require(ggrepel)
require(grDevices)
require(data.table)
require(forcats)
require(plyr)
require(DescTools)
require(sgof)
require(icenReg)
require(readxl)
require(spBayesSurv)
require(viridisLite)
require(splines)
require(gridExtra)
require(gtExtras)
require(webshot2)
require(magick)
require(png)
require(grid)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Icens") #update all [a] if required
require(Icens)
require(interval)
require(survival)
require(sgof)



####Data Import and Cleaning####
#import breakpoints
bp <- read.csv("data/Staph breakpoints.csv", encoding = "UTF-8")
names(bp)
names(bp) <- c("AM.name", "Antimicrobial", "S", "I", "R", "NSbp", "Resource", "Host.species", "Bact.species", "Body.Site", "Comments")

#split signs from values for R (S. value is same as NSbp)
bp$R.value <- unlist(lapply(str_split(bp$R," "), '[', 2))

##original_data <- read.csv(file.choose(), stringsAsFactors=FALSE, header=TRUE, colClasses="character")
staph <- read.csv("data/Staph_fulldata.csv", stringsAsFactors=FALSE, header=TRUE, colClasses= "character", strip.white = TRUE)

#split the column Date and Time.
name1 <- "Date"
name2 <- "Time"
staph <- separate(staph, Date.and.Time, c(name1, name2), sep = " ", extra="merge")

#split the column Date.
name1 <- "Month"
name2 <- "Day"
name3 <- "Year"
staph <- separate(staph, Date, c(name1, name2, name3), sep = "/", extra="merge")
rm(name1,name2,name3)

#filter out body sites to only contain data coded with skin included in the name
staph <- filter(staph, Body.Site %in% c("NASSK", "PUSSK", "SEB", "SKI", "SKIB", 
                                        "SKIF", "SKIR", "SKS", "SKIP", "SKIM", "SK", "SKIS", "ULCSK"))
#body sites
table(staph$Body.Site) #only two body sites given, with majority being labeled SKI

#only include S. intermedius from between 2007 and 2012
unique(staph$Organism.Name) #note spaces after intemedius and the mis-spelling of pseudintermedius
table(staph$Organism.Name, staph$Year) #4 intermedius isolates after 2012
staph$Organism.Name[staph$Organism.Name=="Staphylococcus pseudintermediu    "] <- "Staphylococcus pseudintermedius" #correct mis-spelling
staph <- filter(staph, !(Organism.Name == "Staphylococcus intermedius    " &
                          Year>2011))

#removed antibiotic columns which only consisted of NA values 
staph <- staph[ , colSums(is.na(staph)) < nrow(staph)]

#newer data has MIC columns which only have numbers and no = sign
#need to add the appropriate "=" sign in order to proceed
names(staph)
AM_cols <- 9:49
names(staph)[AM_cols]
staph[,AM_cols] <- sapply(staph[,AM_cols], function(x) str_replace(x, "(^[0-9])", "= \\1"))

#some MIC columns have multiple spaces between sign and value. make just one
staph <- as.data.frame(lapply(staph, str_squish))
#remove any leading or trailing spaces
staph <- as.data.frame(lapply(staph, str_trim))

#some isolates had completely missing susceptibility data. 
#this finds rows with a "= -" or NA in all MIC columns
remove <- staph %>% filter_at(AM_cols, all_vars(. %in% c("= -", NA)))

#remove the rows identified
staph <- anti_join(staph, remove)
rm(remove)

#some isolates have missing MIC values "-" with other MIC having values
#these were replaced with NA
replacedash <- function(df) {
  replace(df, df == "= -", NA)
}

staph <- as.data.frame(sapply(staph, replacedash))

#all analyses: must have at least 30 isolates with MIC values and include only AMs on standard panels
#finding AMs that don't meet the 30 isolate minimum
n_MIC <- sapply(staph[,AM_cols], function(x) sum(!is.na(x)))
AMs_to_exclude <- unique(c(names(n_MIC)[n_MIC<30]))
AMs_to_exclude <- sapply(str_split(AMs_to_exclude, "_"), `[`, 1) #get AM names only


#Are the antibiotics tested on these isolates appropriate and standard for treatment against S. pseudintermedius infection?
#Are there extended panel antibiotics included which may include highly resistant isolates that may bias data?
#utilizing CU AHDC standard companion animal panels for gram positive bacteria
#Standard Panel as of 9/16/2021, COMPGP1F
standard_panel <- c("AMIKAC", "AMOCLA", "AMPICI", "CEFAZO", "CEFOVE", "CEFPOD", "CEPHAL", "CHLORA", "CLINDA", "DOXYCY", "ENROFL", "ERYTH",
                    "GENTAM", "IMIPEN", "MARBOF", "MINOCY", "NITRO", "OXACIL", "PENICI", "PRADOF", "RIFAMP", "TETRA", "TRISUL", "VANCOM")

#on previous standard panel
old_panel <- c("CEFOXI", "CEFTIF", "TICARC", "TICCLA")

#need to remove antibiotics from the dataset that are not used on the standard panels
panels <- unique(c(standard_panel, old_panel)) #create a panel dataset of unique old and current antimicrobials
'%!in%' <- function(x,y)!('%in%'(x,y)) #function allowing us to exclude a factor
names <- colnames(staph[AM_cols]) #take the column names 
nonstandAM <- names[(names %!in% panels)] #match the names of antimicrobials not present in the standard panels
AMs_to_exclude <- unique(c(AMs_to_exclude, nonstandAM))


standard_panel[standard_panel %!in% AMs_to_exclude] #AM on Current Standard and included: AMIKAC, AMOCLA, AMPICI, CEFAZO, CEFOVE, CEFPOD, CEPHAL, CHLORA, CLINDA, DOXYCY, ENROFL, ERYTH,
#GENTAM, IMIPEN, MARBOF, MINOCY, NITRO, OXACIL, PENICI, PRADOF, RIFAMP, TETRA, TRISUL, VANCOM

old_panel[old_panel %!in% AMs_to_exclude] #Old Standard included: CEFOXI, CEFTIF, TICARC, TICCLA

#exclude AMs that don't have at least 30 isolates and are not on the current or prior standard panels
staph <- staph[ , -which(names(staph) %in% AMs_to_exclude)]

#later analysis uses combined sign and MIC, saving dataset to use
susc <- staph

#descriptive results for text
sink("Text Results/descriptive results.txt")
cat("Number of isolates analyzed:  ")
cat(nrow(staph))
sink()


####create survival analysis dataset ####
#split sign from MIC for survival analysis
#I want to start with the metadata columns then I'll add on the AM columns
names(staph)
AM_cols <- 9:36 #new AM_column index after filtering out some AMs
NewData  <- staph %>% select(Bacterial.Species, Organism.Name,
                             Body.Site, Animal.Species, Month, Day, Year, Time)


#split MIC into SIGN and MIC value columns
for(col in AM_cols){
  name <- colnames(staph)[col] #gets the name of the col
  name1 <- paste(name, "SIGN", sep="_") #pastes the name together with "SIGN", separated by "_"
  name2 <- paste(name, "MIC", sep="_") #get "AM_MIC"
  NewColumns <- separate(staph, col, c(name1, name2), sep = " ", extra="merge") #divide sign and number. 
  #sometimes more than one space separated <= and the number. extra="merge" splits at the first space only and puts everything else in the second column
  #paste them together with metadata with cbind
  NewData <- cbind(NewData, NewColumns[,col:(col+1)]) #remember that separate spits out the whole dataframe with the split columns, so I only want to add in the split columns
}
staph_pseud_separated <- NewData
rm(NewColumns, NewData, col, name, name1, name2)

#MIC ranges
#must create start and end column for MIC intervals
#<=4 becomes 0,4
#=4 becomes 2,4
#>4 becomes 4,infinity

#functions need SIGN column (col1) and MIC column (col2)
startcolumn <- function(df, col1, col2) {
  name <- paste(str_split(col1,"_")[[1]][1],"_START", sep="") #create new column name. str_split[[1]][1] returns the AM name
  df %>% dplyr::mutate(!!name:= case_when(.[[get("col1")]] == "<=" ~ 0, 
                                          .[[get("col1")]] == "=" ~ as.numeric(.[[get("col2")]])/ 2, 
                                          .[[get("col1")]] == ">" ~ as.numeric(.[[get("col2")]]), TRUE ~ NA_real_), .keep="none")
}

endcolumn <- function(df, col1, col2) { 
  name <- paste(str_split(col1,"_")[[1]][1],"_END", sep="") 
  df %>% dplyr::mutate(!!name:= case_when(.[[get("col1")]] == "<=" ~ as.numeric(.[[get("col2")]]), 
                                          .[[get("col1")]] == "=" ~ as.numeric(.[[get("col2")]]), 
                                          .[[get("col1")]] == ">" ~ Inf, TRUE ~ NA_real_), .keep="none")
}

names(staph_pseud_separated)
AM_cols_split <- 9:ncol(staph_pseud_separated)
metadata_cols <- 1:8
staphrange <- staph_pseud_separated[,metadata_cols] #metadata (all columns except MICs and SIGNs)

for(i in seq(min(AM_cols_split),ncol(staph_pseud_separated),2)){ #seq creates a sequence of numbers from 9:ncol by 2 (because alternates SIGN, MIC columns)
  col1 <- colnames(staph_pseud_separated)[i] #get the SIGN column
  col2 <- colnames(staph_pseud_separated)[i+1] #get the MIC column
  NewColumnsSt <- startcolumn(df = staph_pseud_separated, col1, col2) #make the start column
  NewColumnsE <- endcolumn(df = staph_pseud_separated, col1, col2) #make the end column
  staphrange <-cbind(staphrange, staph_pseud_separated[,i:(i+1)], NewColumnsSt, NewColumnsE) #save
}

sapply(staphrange, function(x) sum(is.na(x))) #check that number of NA's are the same within each AM SIGN, MIC, start and end col

staph_pseud_separated <- staphrange
rm(staphrange, NewColumnsSt, NewColumnsE, i, col1, col2)
staph <- staph_pseud_separated


####Histograms####
#look at all AMs in the susc dataset (>30 isolates and on standard panel); ok if not meeting MIC variation requirement
#use histograms to visually show above data

#MIC order
MIC_order <- c("<= 0.06", "<= 0.12", "= 0.12", "<= 0.25", "= 0.25", "<= 0.5", "= 0.5", "> 0.5", "<= 1", "= 1", "> 1", "<= 2", "= 2", "> 2", "<= 4", "= 4", "> 4", 
               "<= 8", "= 8", "> 8", "<= 16", "= 16", "> 16", "= 32", "> 32", "= 64", "> 64")
susc[AM_cols] <- map_df(susc[AM_cols], function(x) factor(x, levels=MIC_order))

#create histogram for each column in logstaph
histograms <- imap(susc[,AM_cols], ~{
  ggplot(susc, aes(x = .x))+
    geom_bar()+
    xlab(.y)
})

#make names for saving figs
hist_names <- imap(histograms, ~paste0("Figures and Tables/Histograms/",.y,  ".png")) %>% flatten()

#save figs
walk2(hist_names, histograms, ~ggsave(filename = .x, plot = .y))


####MIC50 and MIC90####
#use staph dataset to examine all AMs with >30 total isolates and on standard panel
#for > MICs, double the MIC value
staph_MICs <- staph %>% select(Year, ends_with(c("MIC","SIGN")))
names(staph_MICs)
AM_cols_staph_MICs <- 2:29 #MIC columns; sign colums + 28
for (i in AM_cols_staph_MICs){
  for(j in 1:nrow(staph_MICs))
    if(isTRUE(staph_MICs[j,i+28]==">")){ #if sign is >
      staph_MICs[j,i] <- as.numeric(staph_MICs[j,i])*2 # double the MIC; otherwise do nothing
    } 
}
#drop SIGN columns
staph_MICs <- staph_MICs %>% select(Year, ends_with(c("MIC")))

#quantiles; have to calculate separately and merge to get a nice dataframe
MIC_q50 <- aggregate(staph_MICs, list(staph_MICs$Year), function(x) quantile(as.numeric(x), probs=c(0.5), na.rm=T))
MIC_q90 <- aggregate(staph_MICs, list(staph_MICs$Year), function(x) quantile(as.numeric(x), probs=c(0.9), na.rm=T))
MIC_qs <- dplyr::full_join(MIC_q50, MIC_q90, by="Year", suffix=c("_50", "_90"))
#clean table
MIC_qs <- MIC_qs %>% select(-"Group.1_50", -"Group.1_90")

#organize columns
MIC_qs <- MIC_qs[,sort(colnames(MIC_qs))]
MIC_qs <- relocate(MIC_qs, "Year")

#rotate; better longer
MIC_qs <- as.data.frame(t(MIC_qs))
colnames(MIC_qs) <- MIC_qs["Year",]
MIC_qs <- MIC_qs[!(row.names(MIC_qs) %in% "Year"),]
MIC_qs$AM <- lapply(str_split(row.names(MIC_qs),"_"), '[',1)
MIC_qs$Q <- paste0("MIC_",lapply(str_split(row.names(MIC_qs),"_"), '[',3))

#later add Log-rank test p-values

rm(staph_MICs, MIC_q50, MIC_q90, AM_cols_staph_MICs)


####tables of MIC values####
#tabulate number of isolates with each MIC value by year and antimicrobial, output list
#use susc (has combo sign-MIC); also staph has had a TETRA 2013 isolate droppe for analysis
MIC_table <- lapply(susc[AM_cols], function(x) as.data.frame(table(susc$Year, x))) 

#pivot wide to have MIC values as columns
MIC_table <- lapply(MIC_table, function(x) pivot_wider(x, names_from=x, values_from=Freq))

#turn 0 into NA (0's come from the individual Year x AM tables if the MIC is present in at least some years but not others)
MIC_table <- lapply(MIC_table, function(x) {x[x==0] <- NA; x})

#drop rows that are all NA in MIC columns (there are 27 MIC columns)
MIC_table <- lapply(MIC_table, function(x) {filter(x, rowSums(is.na(x))<27)})

#then turn NA into blanks
MIC_table <- lapply(MIC_table, function(x) map_df(x,as.character))
MIC_table <- lapply(MIC_table, function(x) {x[is.na(x)]<-""; x})

#change Var1 name to Year
MIC_table <- lapply(MIC_table, function(x) {names(x)[1]<-"Year"; x})

#for each AM, create gt table
#function to create nice tables
MIC_table_gt <- function(x){ #input AM name
  tbl <- MIC_table[[x]] %>% discard(~all(.=="")) %>% tibble() #select AM and drop AM column, drop empty MIC columns, make into tibble
  gt(tbl) %>%
    tab_header(names(MIC_table)[x])
}

#create table for each AM
MIC_gts <- lapply(seq_along(MIC_table), function(x) MIC_table_gt(x))
#list names = AMs
names(MIC_gts) <- names(MIC_table)

#manually add vertical lines for S/I/R and shade MICs that can't be interpreted with current bp (for AM with bps; reference bp)
View(bp)
MIC_gts[["AMIKAC"]] <- gt_add_divider(MIC_gts[["AMIKAC"]], columns="<= 4", side="right", style="dashed") %>% 
  gt_add_divider(columns="= 16", side="left", style="solid") %>%
  tab_style(style=cell_fill(color="#d3d3d3"), locations=cells_body(columns="<= 16"))

MIC_gts[["AMOCLA"]] <- gt_add_divider(MIC_gts[["AMOCLA"]], columns="= 0.25", side="right", style="dashed") %>% 
  gt_add_divider(columns="= 1", side="left", style="solid") %>%
  tab_style(style=cell_fill(color="#d3d3d3"), locations=cells_body(columns="<= 4"))

MIC_gts[["AMPICI"]] <- gt_add_divider(MIC_gts[["AMPICI"]], columns="= 0.25", side="right", style="dashed") #No I

MIC_gts[["CEFAZO"]] <- gt_add_divider(MIC_gts[["CEFAZO"]], columns="<= 2", side="right", style="dashed") %>% 
  gt_add_divider(columns="= 8", side="left", style="solid") %>%
  tab_style(style=cell_fill(color="#d3d3d3"), locations=cells_body(columns=c("<= 8", "<= 4")))

MIC_gts[["CEFOVE"]] <- gt_add_divider(MIC_gts[["CEFOVE"]], columns="= 0.5", side="right", style="dashed") %>% 
  gt_add_divider(columns="= 2", side="left", style="solid")

MIC_gts[["CEFPOD"]] <- gt_add_divider(MIC_gts[["CEFPOD"]], columns="<= 2", side="right", style="dashed") %>% 
  gt_add_divider(columns="= 8", side="left", style="solid") 

MIC_gts[["CEPHAL"]] <- gt_add_divider(MIC_gts[["CEPHAL"]], columns="<= 2", side="right", style="dashed") %>% 
  gt_add_divider(columns="> 4", side="left", style="solid")

MIC_gts[["CHLORA"]] <- gt_add_divider(MIC_gts[["CHLORA"]], columns="= 8", side="right", style="dashed") %>% 
  gt_add_divider(columns="> 16", side="left", style="solid")

MIC_gts[["CLINDA"]] <- gt_add_divider(MIC_gts[["CLINDA"]], columns="= 0.5", side="right", style="dashed") %>% 
  gt_add_divider(columns="= 4", side="left", style="solid") %>%
  tab_style(style=cell_fill(color="#d3d3d3"), locations=cells_body(columns=c("> 1", "> 2"))) #add footnote, these can be interpreted as S/NS

MIC_gts[["DOXYCY"]] <- gt_add_divider(MIC_gts[["DOXYCY"]], columns="<= 0.12", side="right", style="dashed") %>% 
  gt_add_divider(columns="= 0.5", side="left", style="solid") %>%
  tab_style(style=cell_fill(color="#d3d3d3"), locations=cells_body(columns=c("<= 0.25", "<= 2")))

MIC_gts[["ENROFL"]] <- gt_add_divider(MIC_gts[["ENROFL"]], columns="= 0.5", side="right", style="dashed") %>% 
  gt_add_divider(columns="> 2", side="left", style="solid") 

MIC_gts[["ERYTH"]] <- gt_add_divider(MIC_gts[["ERYTH"]], columns="= 0.5", side="right", style="dashed") %>% 
  gt_add_divider(columns="> 4", side="left", style="solid") 

MIC_gts[["GENTAM"]] <- gt_add_divider(MIC_gts[["GENTAM"]], columns="= 4", side="right", style="dashed") %>% 
  gt_add_divider(columns="> 8", side="left", style="solid")

MIC_gts[["MARBOF"]] <- gt_add_divider(MIC_gts[["MARBOF"]], columns="= 1", side="right", style="dashed") %>% 
  gt_add_divider(columns="> 2", side="left", style="solid")

MIC_gts[["MINOCY"]] <- gt_add_divider(MIC_gts[["MINOCY"]], columns="<= 0.5", side="right", style="dashed") %>% 
  gt_add_divider(columns="= 2", side="left", style="solid")

MIC_gts[["NITRO"]] <- gt_add_divider(MIC_gts[["NITRO"]], columns="= 32", side="right", style="dashed") %>% 
  gt_add_divider(columns="> 64", side="left", style="solid") %>%
  tab_footnote(footnote="UTI breakpoint")

MIC_gts[["OXACIL"]] <- gt_add_divider(MIC_gts[["OXACIL"]], columns="<= 0.25", side="right", style="dashed") %>% 
  tab_style(style=cell_fill(color="#d3d3d3"), locations=cells_body(columns=c("<= 2")))

MIC_gts[["PENICI"]] <- gt_add_divider(MIC_gts[["PENICI"]], columns="= 0.12", side="right", style="dashed")

MIC_gts[["PRADOF"]] <- gt_add_divider(MIC_gts[["PRADOF"]], columns="<= 0.25", side="right", style="dashed") %>% 
  gt_add_divider(columns="= 2", side="left", style="solid")

MIC_gts[["RIFAMP"]] <- gt_add_divider(MIC_gts[["RIFAMP"]], columns="<= 1", side="right", style="dashed") %>% 
  gt_add_divider(columns="> 2", side="left", style="solid")

MIC_gts[["TETRA"]] <- gt_add_divider(MIC_gts[["TETRA"]], columns="<= 0.25", side="right", style="dashed") %>% 
  gt_add_divider(columns="= 1", side="left", style="solid") %>%
  tab_style(style=cell_fill(color="#d3d3d3"), locations=cells_body(columns=c("<= 1", "<= 2", "<= 4")))

MIC_gts[["TRISUL"]] <- gt_add_divider(MIC_gts[["TRISUL"]], columns="= 2", side="right", style="dashed")

MIC_gts[["VANCOM"]] <- gt_add_divider(MIC_gts[["VANCOM"]], columns="= 2", side="right", style="dashed")

#make names for saving tables as both html and png
MIC_table_names <- imap(MIC_gts, ~paste0("Figures and Tables/MIC_Tables/",.y,  "_MIC_table.html")) %>% flatten()
walk2(MIC_table_names, MIC_gts, ~gtsave(filename = .x, data = .y)) #save

MIC_table_names <- imap(MIC_gts, ~paste0("Figures and Tables/MIC_Tables//",.y,  "_MIC_table.png")) %>% flatten()
walk2(MIC_table_names, MIC_gts, ~gtsave(filename = .x, data = .y)) #save
#if saving as png gives "Error in s$close() : attempt to apply non-function", re-install webshot2
#NOTE: these MIC tables may not show the same resistance or NS prevalence as the antibiogram for BLs because antibiogram S/NS is aligned for OXA/PEN results and other BL



#### Resistant vs. Susceptible Interpretation####
#for looking at prevalence of resistance, we include all the antimicrobials from the standard panel (even if they had little MIC variation)
#use susc data
#already excluded AMs with < 30 isolates and those not on standard panel
#also exclude AM if there is no breakpoint
AMs_to_exclude_from_Prev <- bp$Antimicrobial[is.na(bp$NSbp)] #missing breakpoints

#keep only breakpoints of included AMs
bp <- filter(bp, Antimicrobial %!in% AMs_to_exclude_from_Prev)
#remove bp for AM not in susc
bp <- filter(bp, Antimicrobial %in% names(susc))

#start table of AM name, class, breakpoints
bp.table <- bp %>% select("AM.name", "Antimicrobial", "S", "I", "R", "NSbp", "Host.species", "Body.Site")

#some AMs in susc do not have a breakpoint; drop them from prevalence analysis
#use AMs_to_include_in_Prev to filter the susc dataset for just the antibiotics used in the prevalence analysis
metadata <- c("Bacterial.Species", "Organism.Name", "Body.Site", "Animal.Species", "Month", "Day", "Year", "Time") #metadata column names
AMs_to_include_in_Prev <- names(susc)[match(bp$Antimicrobial, names(susc))]
susc <- cbind(susc[ , metadata], select(susc, all_of(AMs_to_include_in_Prev)))

#if breakpoints changed over time, older dilution ranges may not cover the new breakpoint
#dilutions that cannot be interpreted should be turned to NA
#if this issue affects many of the isolate MICs for a given AM, it biases the prevalence of resistance
#some year and drug combinations had a low number of isolates with interpretable (non-NA) data which will falsely elevate the resistance percentages of the antibiogram
#a cut off of 30 interpretable isolates per year/drug combination was implemented

#run this before running MIC_to_interpretation() to turn the problematic MIC values for year/drug combos into NA (missing data)
#this checks for dilution errors and changes un-interpretable MIC (due to plate dilutions not covering breakpoints) to NA
#it also replaces year/drug combinations completely with NA if there are <30 isolates that are interpretable or more than 30% of isolates cannot be interpreted
options(warn=1) #print warnings of year/drug combinations replaced with NA
yearexclude <- data.frame()
for (i in match("AMIKAC",colnames(susc)):match("VANCOM",colnames(susc))){ #for each MIC column
  
  bp.index<-match(names(susc)[i],bp$Antimicrobial) #index of which AM bp should be applied to the i MIC column
  
  name <- colnames(susc)[i]
  #Check that the dilutions tested cover the breakpoint appropriately (bp must be > smallest MIC value and < largest MIC value tested). 
  
  #str_match splits the MIC into two columns at "<=". The original string is in the first column returned. If there is a "<=", the MIC is in the second column. if there is not a "<=", it generates an NA in both columns
  min_dilution_too_big = as.numeric(str_match(as.character(susc[,i]), "<=(.*)")[,2])>bp$NSbp[bp.index] #TRUE if minimum dilution is greater than bp; FALSE if it is less than bp; NA if row does not contain a minimum dilution ("<=")
  
  #str_match splits the MIC into two columns at ">". The original string is in the first column. If there is a ">", the MIC is in the second column. if there is not a ">", it generates an NA in both columns
  max_dilution_too_small=as.numeric(str_match(as.character(susc[,i]), ">(.*)")[,2])<bp$NSbp[bp.index] #TRUE if max dilution is less than bp; FALSE if max dilution is greater than bp; NA if does not contain a max dilution (">")
  
  #determine how many isolates in each year have UNinterpretable MICs
  year_dilutions <- as.data.frame(cbind("Year"=susc$Year, min_dilution_too_big, max_dilution_too_small))
  year_dilutions$NotTested <- is.na(susc[i]) #evaluate whether each isolate (row) was tested against the AM
  year_dilutions$min_dilution_too_big <- as.logical(year_dilutions$min_dilution_too_big)
  year_dilutions$max_dilution_too_small <- as.logical(year_dilutions$max_dilution_too_small)
  year_dilutions <- year_dilutions %>% filter(NotTested==FALSE) #remove untested isolates
  year_dilutions$NotOK <- year_dilutions$min_dilution_too_big | year_dilutions$max_dilution_too_small  #identify isolates that cannot be interpreted (either min dilution too big or max dilution too small)
  
  #years that have <30 isolates that can be interpreted (NA in NotOK)
  year_exclude1 <- year_dilutions %>% group_by(Year) %>% filter(sum(is.na(NotOK))<30) %>% select(Year) %>% unique()
  #table(year_dilutions$Year, year_dilutions$NotOK, useNA="always") #check
  
  #years that have more than X% of tested isolates that cannot be interpreted (TRUE in NotOK and FALSE in NotTested)
  threshold <- 0.3
  year_exclude2 <- year_dilutions %>% group_by(Year) %>% filter((sum(NotOK, na.rm=T)/n())>threshold) %>% select(Year) %>% unique()
  #prop.table(table(year_dilutions$Year, year_dilutions$NotOK, useNA="always"),margin=1) #check
  
  year_exclude <- unique(unlist(c(year_exclude1, year_exclude2)))
  yearexclude <- rbind(yearexclude, name, year_exclude) #saving for later analysis
  
  #report affected AM
  if (nrow(year_exclude1)>0){
    warning(colnames(susc)[i], ": ",
            "Years ", year_exclude1, " had fewer than 30 isolates with interpretable MICs. ",
            "All isolate MICs in these years are replaced by NA. ")
  }
  if (nrow(year_exclude2)>0){
    warning(colnames(susc)[i], ": ",
            "Years ", year_exclude2, " had more than ", threshold*100, "% of isolates with UNinterpretable MICs. ",
            "All isolate MICs in these years are replaced by NA. ")
  }
  
  #for excluded years, turn all MIC for the given AM into NA
  susc[susc$Year %in% year_exclude,i] <- as.character(NA)
  
  #report number of remaining interpretation issues
  n_min_dil_too_big <- with(year_dilutions, sum(min_dilution_too_big[!(Year %in% year_exclude)]==TRUE, na.rm=T))
  n_max_dil_too_small <- with(year_dilutions, sum(max_dilution_too_small[!(Year %in% year_exclude)]==TRUE, na.rm=T))
  if (any(year_dilutions$min_dilution_too_big[!(year_dilutions$Year %in% year_exclude)]==TRUE, na.rm=TRUE)){
    warning("minimum dilution greater than breakpoint for ", colnames(susc)[i], ". ", 
            n_min_dil_too_big," values replaced by NA \n")
  }
  
  if (any(year_dilutions$max_dilution_too_small[!(year_dilutions$Year %in% year_exclude)]==TRUE, na.rm=TRUE)){
    warning("maximum dilution less than breakpoint for ", colnames(susc)[i], ". ", 
            n_max_dil_too_small," values replaced by NA \n")
  }
  
  #turn NA into FALSE for indexing
  min_dilution_too_big <- replace_na(min_dilution_too_big, FALSE)
  max_dilution_too_small <- replace_na(max_dilution_too_small, FALSE)
  
  #replace MIC values with NA when it is a min dilution that is too big or a max dilution that is too small
  susc[min_dilution_too_big,i] <- as.character(NA)
  susc[max_dilution_too_small,i] <- as.character(NA)
}

#function to interpret MIC as resistant or susceptible
MIC_to_interpretation <- function (data, index, bp){
  #required packages
  require (stringr)
  require (dplyr)
  require (tidyr)
  require (arules)
  
  data[,index]<-sapply(data[,index], as.character) #first convert MIC to strings
  data[,index] <-sapply(data[,index], function(x) gsub(" ", "", x, fixed = TRUE)) #remove whitespace in the MIC columns
  
  for (i in index){ #for each MIC column
    bp.index<-match(names(data)[i],bp$Antimicrobial) #index of which AM bp should be applied to the i MIC column
    
    #Check that the dilutions tested cover the breakpoint appropriately (bp must be > smallest MIC value and < largest MIC value tested). Print warning if the bp does not fall within the tested dilutions
    dilution_bp_check1=as.numeric(str_match(data[,i], "<=(.*)")[,2])
    #this splits the MIC into two columns at "<=". The original string is in the first column returned. If there is a "<=", the MIC is in the second column. if there is not a "<=", it generates an NA in both columns
    
    dilution_bp_check2=as.numeric(str_match(data[,i], ">(.*)")[,2])
    #this splits the MIC into two columns at ">". The original string is in the first column. If there is a ">", the MIC is in the second column. if there is not a ">", it generates an NA in both columns
    
    #if any of the minimum MIC values are greater than the breakpoint, warn
    if (any(dilution_bp_check1>bp$NSbp[bp.index], na.rm=TRUE)){
      warning("minimum dilution greater than breakpoint for ", colnames(data)[i], "\n")
    }
    
    #if any of the maximum MIC values are greater than the breakpoint, warn
    if (any(dilution_bp_check2<bp$NSbp[bp.index], na.rm=TRUE)){
      warning("maximum dilution less than breakpoint for ", colnames(data)[i], "\n")
    }
    
    
    #next remove >, <=, =; then make MIC numerica
    data[,i]<-as.numeric(
      str_replace_all(data[,i], c(
        "<=" = "", 
        "=" = "",
        ">" = "1000"))) #treat any ">" MIC values as large number to clearly distinguish from "-" MIC values in case different dilution series are tested for the same AM.
    
    #categorize each MIC value as S or NS based on the breakpoint
    if (is.na(bp.index)==TRUE){ #if there bp.index for the i MIC column is NA, then there is no bp for that AM
      data[,i]<-NA #replace all MIC values with NA
      warning("missing breakpoint for ", colnames(data)[i], "\n") #warn that there is no bp for a column
    } else{ #if the bp.index is not NA, then there is a bp for that AM
      data[,i]<-discretize(data[,i], method="fixed", 
                           breaks=c(-Inf, bp$NSbp[bp.index], Inf), 
                           right=TRUE, labels=c("FALSE", "TRUE")) #discretize the MIC values at the bp. Intervals are right closed: (-Inf,bp]; (bp, Inf). Hence, MIC values >bp are given "True" (NS) and MIC values <=bp are given "False" (S)
      data[,i]=as.logical(data[,i]) #make logical
    } 
    
  }
  
  
  
  data #return the discretized data
}

interp_index <- match("AMIKAC", names(susc)):match("VANCOM", names(susc)) ##AM columns in susc
MIC_interp <- MIC_to_interpretation(data = susc, index = interp_index, bp = bp)
#cleaning the environment, but if you want to see the exact year/drug combinations excluded or problematic drug dilutions, don't remove
rm(year_dilutions, year_exclude, year_exclude1, year_exclude2, yearexclude, bp.index, i, interp_index,
   max_dilution_too_small, min_dilution_too_big, n_max_dil_too_small, n_min_dil_too_big, name, threshold)

#since 2007 was fewer than 30 isolates for all antimicrobials, drop from the prevalence dataframe
MIC_interp <- MIC_interp %>% filter(Year!="2007")

#VETO1S standards for beta-lactam resistance
#isolates resistant to oxacillin are considered to be automatically resistant to all b-lactams (except ceftaroline), including penicillins, b-lactam combos, cephems, and carbapenams (Ed 6 Table 2C-1 note 13)
#Penicillin-R isolates are R to other penicillinase-labile penicillins (ampici, ticarc) (Ed 6 Table 2C-1 note 9)
#isolates susceptible to penicillin are susceptible to other beta-lactams, including penicillinase-labile and penicillinase-stable agents (Ed 6 Table 2C-1 note 9)
#Oxacillin-susceptible isolates  can be considered susceptible to: b-lactam combos (amocla, ticcla), cephems (generations I, II, III, IV), carbapenems (IMIPEN) (Ed 6 Table 2C-1 note 13)
#change MIC_interp to reflect this

names(MIC_interp) #available BLs

#first see the problem - OXACIL-R isolates that are S to one or more other b-lactams. (no bp for ticcla or ticarc or CEFTIF)
OXA_R_probs <- MIC_interp %>% filter(OXACIL==TRUE & (PENICI==FALSE | AMPICI==FALSE | AMOCLA==FALSE | CEFAZO==FALSE | CEFOVE==FALSE | CEFPOD==FALSE | CEPHAL==FALSE))
View(OXA_R_probs)
table(MIC_interp$PENICI, MIC_interp$OXACIL, dnn=c("PEN", "OXA"), useNA="always") ##there is 1 OXA-R that is PEN-S. Also 3 OXA-R that are PEN-NA, these could also be changed to PEN-R but we haven't here
table(MIC_interp$AMOCLA, MIC_interp$OXACIL, dnn=c("AMC", "OXA"), useNA="always") #183 isolates OXA-R AMC-S; 50 isolates OXA-R AMC-NA
table(MIC_interp$AMPICI, MIC_interp$OXACIL, dnn=c("AMP", "OXA"), useNA="always") #12 isolates OXA-R AMP-S; 3 isolates OXA-R AMP-NA
table(MIC_interp$CEFAZO, MIC_interp$OXACIL, dnn=c("CEFA", "OXA"), useNA="always") #467 isolates OXA-R CEFA-S; 51 isolates OXA-R CEFA-NA
table(MIC_interp$CEFOVE, MIC_interp$OXACIL, dnn=c("CEFO", "OXA"), useNA="always") #44 isolates OXA-R CEFO-S; 33 isolates OXA-R CEFO-NA
table(MIC_interp$CEFPOD, MIC_interp$OXACIL, dnn=c("CEFP", "OXA"), useNA="always") #298 isolates OXA-R CEFP-S; 9 isolates OXA-R CEFP-NA
table(MIC_interp$CEPHAL, MIC_interp$OXACIL, dnn=c("CEP", "OXA"), useNA="always") #522 isolates OXA-R CEF-S; 17 isolates OXA-R CEP-NA

#OXA-R isolates = resistant to pencillins, B lactam combos, cephems, carbapenems
#changing S to R: as.logical() turns any value > 0 into TRUE and TRUE has value of 1. so if the other b-lactam is S or R, adding TRUE and using as.logical() makes it a TRUE (R)
MIC_interp <- MIC_interp  %>% mutate(PENICI = case_when(OXACIL==TRUE ~ as.logical(PENICI+TRUE),
                                                        TRUE ~ PENICI),
                                     AMOCLA = case_when(OXACIL==TRUE ~ as.logical(AMOCLA+TRUE),
                                                        TRUE ~ AMOCLA),
                                     AMPICI = case_when(OXACIL==TRUE ~ as.logical(AMPICI+TRUE),
                                                        TRUE ~ AMPICI),
                                     CEFAZO = case_when(OXACIL==TRUE ~ as.logical(CEFAZO+TRUE),
                                                        TRUE ~ CEFAZO),
                                     CEFOVE = case_when(OXACIL==TRUE ~ as.logical(CEFOVE+TRUE),
                                                        TRUE ~ CEFOVE),
                                     CEFPOD = case_when(OXACIL==TRUE ~ as.logical(CEFPOD+TRUE),
                                                        TRUE ~ CEFPOD),
                                     CEPHAL = case_when(OXACIL==TRUE ~ as.logical(CEPHAL+TRUE),
                                                        TRUE ~ CEPHAL)) #the last statement, TRUE ~ is like the "else", do this to anything not addressed in the conditions

#check that mutate worked, should be zero isolates that are OXA-R and other beta-lactam-S
table(MIC_interp$PENICI, MIC_interp$OXACIL, dnn=c("PEN", "OXA"), useNA="always") 
table(MIC_interp$AMOCLA, MIC_interp$OXACIL, dnn=c("AMC", "OXA"), useNA="always")
table(MIC_interp$AMPICI, MIC_interp$OXACIL, dnn=c("AMP", "OXA"), useNA="always")
table(MIC_interp$CEFAZO, MIC_interp$OXACIL, dnn=c("CEFA", "OXA"), useNA="always")
table(MIC_interp$CEFOVE, MIC_interp$OXACIL, dnn=c("CEFO", "OXA"), useNA="always") 
table(MIC_interp$CEFPOD, MIC_interp$OXACIL, dnn=c("CEFP", "OXA"), useNA="always")
table(MIC_interp$CEPHAL, MIC_interp$OXACIL, dnn=c("CEP", "OXA"), useNA="always")

#now for PEN-S: these isolates should be susceptible to all other b-lactams
PEN_probs <- MIC_interp %>% filter(PENICI==FALSE & (OXACIL==T | AMPICI==T | AMOCLA==T | CEFAZO==T | CEFOVE==T | CEFPOD==T | CEPHAL==T)) #only AMPICI has inconsistent data for this isolate
View(PEN_probs)
table(MIC_interp$PENICI, MIC_interp$AMPICI, dnn=c("PEN", "AMP"), useNA="always")  #1 isolate, S to PENICI, R to AMPICI. AMPICI and PENICI should be consistent with each other (both are penicillinase-labile penicillins)
#since the PENICI breakpoint is for humans and the AMPICI breakpoint is for dogs, we preferentially took the AMPICI breakpoint interpretation and turned the PENICI to R. 

MIC_interp <- MIC_interp  %>% mutate(PENICI = case_when(AMPICI==TRUE ~ as.logical(PENICI+TRUE),
                                                        TRUE ~ PENICI))
table(MIC_interp$PENICI, MIC_interp$AMPICI, dnn=c("PEN", "AMP"), useNA="always") #changed

#OXA-S: isolates should be S to b-lactam combos, cephems, carbapenems
OXA_S_probs <- MIC_interp %>% filter(OXACIL==F & (AMOCLA==T | CEFAZO==T | CEFOVE==T | CEFPOD==T | CEPHAL==T))
View(OXA_S_probs)
table(MIC_interp$AMOCLA, MIC_interp$OXACIL, dnn=c("AMC", "OXA"), useNA="always") #1 OXA-S that are AMC-R
table(MIC_interp$CEFAZO, MIC_interp$OXACIL, dnn=c("CEFA", "OXA"), useNA="always") #1 isolate
table(MIC_interp$CEFOVE, MIC_interp$OXACIL, dnn=c("CEFO", "OXA"), useNA="always") #8 isolates
table(MIC_interp$CEFPOD, MIC_interp$OXACIL, dnn=c("CEFP", "OXA"), useNA="always") #1 isolate
table(MIC_interp$CEPHAL, MIC_interp$OXACIL, dnn=c("CEP", "OXA"), useNA="always") #3 isolates

#here we need to multiply the logicals, because we want the "T" (R) to be turned to "F" (S)
MIC_interp <- MIC_interp  %>% mutate(AMOCLA = case_when(OXACIL==FALSE ~ as.logical(AMOCLA*F),
                                                        TRUE ~ AMOCLA),
                                     CEFAZO = case_when(OXACIL==FALSE ~ as.logical(CEFAZO*F),
                                                        TRUE ~ CEFAZO),
                                     CEFOVE = case_when(OXACIL==FALSE ~ as.logical(CEFOVE*F),
                                                        TRUE ~ CEFOVE),
                                     CEFPOD = case_when(OXACIL==FALSE ~ as.logical(CEFPOD*F),
                                                        TRUE ~ CEFPOD),
                                     CEPHAL = case_when(OXACIL==FALSE ~ as.logical(CEPHAL*F),
                                                        TRUE ~ CEPHAL)) 
#check that all replaced
table(MIC_interp$AMOCLA, MIC_interp$OXACIL, dnn=c("AMC", "OXA"), useNA="always") 
table(MIC_interp$CEFAZO, MIC_interp$OXACIL, dnn=c("CEFA", "OXA"), useNA="always") 
table(MIC_interp$CEFOVE, MIC_interp$OXACIL, dnn=c("CEFO", "OXA"), useNA="always")
table(MIC_interp$CEFPOD, MIC_interp$OXACIL, dnn=c("CEFP", "OXA"), useNA="always")
table(MIC_interp$CEPHAL, MIC_interp$OXACIL, dnn=c("CEP", "OXA"), useNA="always")

rm(OXA_R_probs, OXA_S_probs, PEN_probs)

##### Antibiogram #####
#generate antibiogram: prevalence of NS for each AM; data grouped by user-defined variable

#data: dataset to analyze. must be in logical form (F=S, T=NS), i.e. from MIC_to_interpretation
#index: indices of AM columns to analyze.
#group: index of the grouping variable column used to split the database (e.g. year)
summary(MIC_interp)

antibiogram <- function (data, index, group){
  #required packages
  require(dplyr)
  require(tidyr)
  
  abgm <- group_by(data[,c(group, index)] , data[,group]) %>%
    summarise_all(list(Prevalence=~mean(.,na.rm=TRUE),Number_Isolates_Tested=~sum(!is.na(.)))) #drop all columns except grouping column and AM columns, group data by the grouping column. Summarise each column with prevalence (# NS / # tested = mean of the logical variable with na.rm=TRUE) and number of isolates tested. This returns the prevalence and number tested within each group value. Columns include: AM_Prevalence and AM_Number_Isolates_Tested, with "AM" being each AM name
  
  #overall prevalences (across all group values)
  overall_prev <- summarise_all(data[,c(group,index)], list(Prevalence=~mean(.,na.rm=TRUE),Number_Isolates_Tested=~sum(!is.na(.)))) #same as above but without grouping
  overall_prev <- cbind("data[, group]"="Overall", overall_prev) #add in a 'group' value of "Overall"
  abgm <- rbind(abgm, overall_prev) #add overall values to the grouped antibiogram
  
  #if an AM is tested but there is no breakpoint, num of isolates tested will show 0 since the data (logical form) will have all NA. doesn't differentiate from AM that were not tested and were in dataset but with no MIC values (logical data also NA and num of isolates tested will be 0).
  #warn if overall num of isolates tested columns (second half of antibiogram) are equal to 0. These AM could be excluded from 'index'
  num_tested_col=c(((ncol(overall_prev)-1)/2+2):ncol(overall_prev)) #column indices for the "Number_Isolates_Tested" columns in the overall antibiogram
  if (any(overall_prev[num_tested_col]==0)){ #if any of these columns equal 0, warn
    warning("All interpretations NA for:", "\n", paste(str_split(names(overall_prev[num_tested_col])[overall_prev[num_tested_col]==0], "_", simplify = TRUE)[,1], collapse=", "), "\n", "May be missing breakpoint or not tested.") #list out AM names of columns with missing bp or no isolates tested
  }
  
  #sort alphabetically to put prevalence next to number tested for each drug.
  abgm <- abgm[,order(names(abgm))]
  
  abgm #return antibiogram
  
}

MIC_interp$Year <- as.numeric(MIC_interp$Year) #make Year numeric
a_index <- match("AMIKAC",colnames(MIC_interp)):match("VANCOM",colnames(MIC_interp)) #AM columns
years <- match("Year", colnames(MIC_interp)) #year column number
antibiogram_staph <- antibiogram(data = MIC_interp, index = a_index, group = years)
#rename the grouping column to year
names(antibiogram_staph)
antibiogram_staph <- antibiogram_staph %>% dplyr::rename(Year="data[, group]")

rm(years)

####Refined Antibiogram Table####
antibiogram_round <- antibiogram_staph
#reorganizing the columns, put year and total isolates first
antibiogram_round <- antibiogram_round %>% relocate(Year, .before = "AMIKAC_Number_Isolates_Tested")
antibiogram_round <- antibiogram_round %>% relocate(Year_Number_Isolates_Tested, .before = "AMIKAC_Number_Isolates_Tested")
antibiogram_round <- antibiogram_round[ , -which(names(antibiogram_round) %in% "Year_Prevalence")] ##remove redundant Year column
##round Prevalence columns to 2 significant digits and make a percentage
antibiogram_round[,seq(4,ncol(antibiogram_round),2)] <- signif(antibiogram_round[,seq(4,ncol(antibiogram_round),2)], digits=2)*100


Table1 <- select(antibiogram_round, Year, Year_Number_Isolates_Tested) #start with Year category and total isolates
#collapse columns of 'number tested' and 'percent resistant' into one column for each AM: %(N)
for (i in seq(3,(ncol(antibiogram_round)),2)){ #for every other column in the overall_abgm (each AM), starting with the Num_Isolates columns
  tab <- apply(antibiogram_round[,c(i+1,i)],1, paste, collapse=" (") #paste the column of Prevalence with column of Number isolates tested; put N tested in parentheses
  tab <- as.data.frame(paste(tab, ")", sep="")) #close parentheses
  colnames(tab) <- str_split(colnames(antibiogram_round)[i], "_")[[1]][1] #make column name the AM
  Table1 <- cbind(Table1, name=tab) #add to Table1
  rm(tab)
}
rm(i)
#rename columns; bp contains AM abbreviations and full names
Table1 <- dplyr::rename_with(Table1, ~bp$AM.name[which(bp$Antimicrobial==.x)], .cols=bp$Antimicrobial)
Table1 <- dplyr::rename(Table1, "Number of Isolates"="Year_Number_Isolates_Tested")

#NaN to '-' for %R when no isolates tested
Table1 <- sapply(Table1, function(x) {gsub("NaN", "-", x)})
Table1 <- as.data.frame(t(Table1)) #swap rows and columns to make longer/not as wide
colnames(Table1) <- Table1["Year",] #rename columns, then drop Year
Table1 <- Table1[!(row.names(Table1) %in% "Year"),]

#later add indicator for Cochran-Armitage test

####Prevalence Plots####
#make a dataframe of only prevalence data for graphing
#for graphs removed the last row of the overall prevalence
#drop the overall data
prev <- antibiogram_staph %>% filter(Year!="Overall")
prev <- prev %>% select(contains("_Prevalence")) #drop columns on number of isolates
#rename the columns for graph
#use str_split() on names(prev) to remove the "_Prevalence"
colnames(prev) <- sapply(str_split(names(prev), "_"), "[")[1,]

#Separate into Beta-Lactams and Other ABX
beta <- prev %>% select("Year", "AMOCLA", "AMPICI","CEFAZO", "CEFOVE", "CEFPOD",
                        "CEPHAL","OXACIL", "PENICI")
other <- prev %>% select("Year", "AMIKAC", "CHLORA","CLINDA", "DOXYCY", "ENROFL",
                         "ERYTH", "GENTAM", "MARBOF","MINOCY", "NITRO", "PRADOF", 
                         "RIFAMP", "TETRA", "TRISUL", "VANCOM")

#long format for plotting
beta_long=reshape2::melt(beta, id.vars=c("Year"))
other_long=reshape2::melt(other, id.vars=c("Year"))

#AM class dataframe
AM <- c("AMIKAC", "AMOCLA", "AMPICI", "CEFAZO", "CEFOVE", "CEFPOD", "CEPHAL", "CHLORA", 
        "CLINDA", "DOXYCY", "ENROFL", "ERYTH", "GENTAM", "MARBOF", "MINOCY", "NITRO", "OXACIL", 
        "PENICI", "PRADOF", "RIFAMP", "TETRA", "TRISUL", "VANCOM")

##Class designations
#A = Aminoglycides
#AN = Ansamycins
#B = Beta-Lactam
#F = Fluoroquinolones
#G = Glycopeptides
#M = Macrolides
#N = Nitrofurans
#P = Phenicols
#S = Sulfonamides
#T = Tetracyclines

Class <- c("A", "B", "B","B", "B", "B", "B", "P", "M","T", "F",
           "M", "A", "F", "T", "N", "B", "B", "F", "AN", "T", "S", "G")

AM_Class <- data.frame(AM, Class)

AM_Class$color <- c("#E69F00FF", "#F0E442FF", "#E69F00FF", "#56B4E9FF", "#009E73FF","#009E73FF", "#56B4E9FF",
                               "#0072B2FF", "#F0E442FF", "#CC79A7FF", "#009E73FF", "#F0E442FF", "#E69F00FF",
                               "#009E73FF", "#CC79A7FF", "#999933FF", "#E69F00FF", "#E69F00FF", "#009E73FF",
                               "#56B4E9FF", "#CC79A7FF", "#D55E00FF", "#882255FF")
                               
#color blind palette c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#882255", "#999933")
#Assigned colors for all classes
##Beta-Lactams
#P = Penicillin alpha("#E69F00", 1)
#1 = 1st gen cephalosporin alpha("#56B4E9", 1)
#3 = 3rd gen cephalosporin alpha("#009E73", 1)
#C = Combinations alpha("#F0E442", 1)

##Other
#A = Aminoglycides alpha("#E69F00", 1)
#AN = Ansamycins alpha("#56B4E9", 1)
#F = Fluoroquinolones alpha("#009E73", 1)
#G = Glycopeptides alpha("#882255", 1)
#M = Macrolides alpha("#F0E442", 1)
#N = Nitrofurans alpha("#999933", 1)
#P = Phenicols alpha("#0072B2", 1)
#S = Sulfonamides alpha("#D55E00", 1)
#T = Tetracyclines alpha("#CC79A7", 1)

##Beta-Lactams
colors <- c(alpha("#F0E442", 1), #AMOCLA C
            alpha("#E69F00", 1), #AMPICI P
            alpha("#56B4E9", 1), #CEFAZO 1
            alpha("#009E73", 1), #CEFOVE 3
            alpha("#009E73", 1), #CEFPOD 3
            alpha("#56B4E9", 1), #CEPHAL 1
            alpha("#E69F00", 1), #OXACIL P
            alpha("#E69F00", 1)) #PENICI P
linetype <- c(1, 1, 1, 1, 2, 2, 2, 3)

prev <- ggplot(beta_long, aes(Year, value)) + 
  geom_line(aes(colour = variable), linewidth=1.5, linetype=rep(linetype,rep(13,length(linetype))))+ 
  geom_point(aes(colour = variable), size=1.5)+
  #geom_label_repel(data = subset(beta_long, Year==2008), aes(label=variable), point.padding=0.25, min.segment.length = 1, nudge_x=0.5, size=5, max.overlaps = Inf)+ ##CC: option to add some missing labels:  max.overlaps = Inf; nudge not working/labels disappearing
  geom_label_repel(data = subset(beta_long, Year==2020), aes(label=variable), point.padding=0.25, min.segment.length =1, nudge_x=1, size=5, max.overlaps = Inf)+  ##CC: option to add some missing labels:  max.overlaps = Inf
  labs(x="Year", y="Proportion of Non-Susceptible")+
  scale_color_manual("Antimicrobial",values=colors)+
  scale_linetype_manual("Antimicrobial", values=linetype)+
  scale_x_continuous(limits=c(2008,2021), breaks=seq(2008, 2020, 1), expand=c(0,0.2))+
  scale_y_continuous(limits=c(0,1.1), breaks=seq(0, 1.0, 0.1))+
  theme_classic(base_size = 20)+
  theme(legend.position="top", plot.margin=unit(c(10,10,10,10), "pt"), plot.title=element_text(size=24, hjust=0.5), axis.text=element_text(size=22))+
  guides(color=guide_legend(nrow=2, keywidth=4, override.aes = list(linetype = linetype)))

png('Figures and Tables/Prevalence/betaprev.png', width=20, height=10, units='in', res=600)
plot.new()
prev
dev.off()

##Others
colors <- c(alpha("#E69F00", 1), #AMIKAC A 
            alpha("#0072B2", 1), #CHLORA P 
            alpha("#F0E442", 1), #CLINDA M 
            alpha("#CC79A7", 1), #DOXYCY T 
            alpha("#009E73", 1), #ENROFL F 
            alpha("#F0E442", 1), #ERYTH  M 
            alpha("#E69F00", 1), #GENTAM A 
            alpha("#009E73", 1), #MARBOF F 
            alpha("#CC79A7", 1), #MINOCY T 
            alpha("#999933", 1), #NITRO N 
            alpha("#009E73", 1), #PRADOF F 
            alpha("#56B4E9", 1), #RIFAMP AN 
            alpha("#CC79A7", 1), #TETRA T 
            alpha("#D55E00", 1), #TRISUL S 
            alpha("#882255", 1)) #VANCOM G 
linetype <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 3, 1, 3, 1, 1)

prev2 <- ggplot(subset(other_long), aes(Year, value)) +
  geom_line(aes(colour = variable), size=1.5, linetype=rep(linetype,rep(13,15)))+ ##also gave linetype length error. changed form rep(13,15) ##CMC: got the same error so I changed it back and it worked?
  geom_point(aes(colour=variable), size=1.5)+
  #geom_label_repel(data=subset(other_long, Year==2008), aes(label=variable), point.padding=0.25, min.segment.length = 1, nudge_x=0, size=5)+
  geom_label_repel(data=subset(other_long, Year==2020), aes(label=variable), point.padding=0.25, min.segment.length =1, nudge_x=1, size=5)+
  labs(x="Year", y="Proportion of Non-Susceptible")+
  scale_color_manual("Antimicrobial",values=colors)+
  scale_linetype_manual("Antimicrobial", values=linetype)+
  scale_x_continuous(limits=c(2008,2021), breaks=seq(2008, 2020, 1), expand=c(0,0.2))+
  scale_y_continuous(limits=c(0,1.1), breaks=seq(0, 1.0, 0.1))+
  theme_classic(base_size = 20)+
  theme(legend.position="top", plot.margin=unit(c(10,10,10,10), "pt"), plot.title=element_text(size=24, hjust=0.5), axis.text=element_text(size=22))+
  guides(color=guide_legend(nrow=2, keywidth=4, override.aes = list(linetype = linetype)))
png('Figures and Tables/Prevalence/otherprev.png', width=20, height=10, units='in', res=600)
plot.new()
prev2
dev.off()

#arrange into Figure2A/B
png('Figures and Tables/Prevalence/Figure 2.png', width=20, height=20, units='in', res=600)
plot.new()
ggarrange(prev, prev2, nrow=2, labels=c("A", "B"), font.label=list(size=30, face="bold"))
dev.off()

rm(prev, prev2, beta, other, beta_long, other_long, colors, linetype)


####MDR profiles####
#For each isolate, identifies the resistance profile (all the drugs to which the isolate is resistant), tabulates the number of drugs to which it is resistant, identifies the class resistance profile (all the AM classes to which the isolate is resistant), the number of drug classes to which it is resistant

#data: dataset of AM interpretations. Note that if there is no breakpoint, then a dataset of AM interpretations will exclude AM that were tested but don't have a breakpoint.
#index: indices of AM columns
#AM_class is a two column dataframe. first column (name = "AM") is the name of the antimicrobial (must match AM names from data). second column (name = "Class") is the assigned class (number or character that can be sorted alphanumerically)

mdr_profile <- function (data, index, AM_class){
  #required packages
  require(stringr)
  require(mgsub)
  
  z <- data[,index] #look at only AM columns
  
  mdr_profile <- data.frame(Profile=rep("", dim(z)[1]), NumRes=numeric(dim(z)[1]), NumClass=numeric(dim(z)[1])) #create dataframe to store mdr profiles
  profile <- rep("", dim(z)[1]) #Initializes a blank character as long as the number of isolates in the data
  
  
  mdr_profile[,2] <- rowSums(z, na.rm=TRUE) #number of drug resistances for each isolate
  
  for(i in 1:ncol(z)) { #for each AM 
    z[which(z[,i]==TRUE), i] <- as.character(colnames(z)[i]) #if resistant, paste AM name into cell instead of True/False (interpreted MIC)
    z[which(z[,i]==FALSE), i] <- as.character("") #if susceptible, paste blank into cell
    z[is.na(z[,i]), i] <- as.character("") #those that weren't tested or have no bp (have NA instead of T/F), paste blank into cell
    profile <- paste(profile, z[,i], sep = " ") #paste all AM columns together to get a character string of all the AM to which each isolate (row) is resistant
  }
  profile <- trimws(profile, which="both") #trim whitespace
  mdr_profile[,1] <- profile #store the resistance profile
  
  
  #warn if the data contains an AM that is not in the AM_class--cannot proceed with class resistance tabulation
  if (any(is.na(match(colnames(z), AM_class$AM)))){ #if there are any AM names in the data that do not have a corresponding entry in AM_class, warn:
    warning("Missing AM class for ", colnames(z)[is.na(match(colnames(z), AM_class$AM))], "\n Inappropriate data from this column is pasted in results.")
  }
  
  #tabulate the class resistance profile and number of classes to which each isolate is resistant
  mdr_profile$Classes <- mgsub(mdr_profile$Profile, as.character(AM_class$AM), as.character(AM_class$Class))#in a new column, replace the name of each AM from the resistance profile with the name of the corresponding class
  
  mdr_profile$Classes.sort <- sapply(sapply(str_extract_all(mdr_profile$Classes, "[a-z0-9A-Z]+"), unique), sort) #list unique classes for each isolate then sort alphanumerically
  
  mdr_profile$NumClass <- sapply(lapply(str_extract_all(mdr_profile$Classes, "[a-z0-9A-Z]+"), unique), length) #number of unique class resistances
  
  mdr_profile #output
}

#generate mdr profiles
MDR <- mdr_profile(data = MIC_interp, index = a_index, AM_class = AM_Class)
sink("Text Results/descriptive results.txt", append=T) ##start saving results to file
cat("\n number of isolates in the prevalence analysis:  ")
cat(nrow(MIC_interp))
cat("\n pan-susceptible prevalence:  ")
cat(sum(MDR$NumClass==0))
cat("\n percent of pan-susceptible isolates: ")
cat(round((sum(MDR$NumClass==0)/nrow(MDR))*100,),"%")
cat("\n number of isolates that are MDR (resistant to 3 or more classes):  ") ##using cat rather than print eliminates printed line numbers
cat(sum(MDR$NumClass>=3))
cat("\n percent of isolates that are MDR:  ") ##\n prints a new line
cat(round((sum(MDR$NumClass>=3)/nrow(MDR))*100,),"%")
cat("\n >=5 class prevalence:  ")
cat(sum(MDR$NumClass>=5))
cat("\n percent of >=5 class isolates: ")
cat(round((sum(MDR$NumClass>=5)/nrow(MDR))*100,),"%")
sink()

#MDR Graphical Representations
#How many of the isolates are resistant to multiple classes?
MDRpct <- MDR %>% dplyr::count(NumClass) %>% mutate(pct = n / sum(n) * 100) #calculate number and percent of isolates resistant to N classes
MDRpct$NumClass <- as_factor(MDRpct$NumClass)
MDRpct$NumClass <- fct_relevel(MDRpct$NumClass, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")

#plot % isolates by N resistance classes
MDRpct_fig <- ggplot(MDRpct, aes(x = as.numeric(NumClass), pct)) + 
  geom_bar(stat = 'identity') + 
  xlab("Number of Non-Susceptible\nAntimicrobial Classes") + 
  ylab("Percent of Isolates") + 
  theme_bw() +
  scale_x_discrete(limits = levels(MDRpct$NumClass))+
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0,25))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("Figures and Tables/MDR/Figure 1.png", MDRpct_fig)


#How many isolates are resistant to each specific antimicrobial class?
##use the Classes.sort column in MDR
View(MDR$Classes.sort) #is a list

#identify presence/absence of each class resistance in each isolate
MDRclass <- data.frame(matrix(NA, nrow=nrow(MIC_interp),
                              ncol=length(unique(AM_Class$Class)))) #create data frame
colnames(MDRclass) <- unique(AM_Class$Class)
MDRclass$A <- unlist(lapply(MDR$Classes.sort, function(x)
  "A" %in% x))
MDRclass$B <- unlist(lapply(MDR$Classes.sort, function(x)
  "B" %in% x))
MDRclass$P <- unlist(lapply(MDR$Classes.sort, function(x)
  "P" %in% x))
MDRclass$M <- unlist(lapply(MDR$Classes.sort, function(x)
  "M" %in% x))
MDRclass$T <- unlist(lapply(MDR$Classes.sort, function(x)
  "T" %in% x))
MDRclass$F <- unlist(lapply(MDR$Classes.sort, function(x)
  "F" %in% x))
MDRclass$N <- unlist(lapply(MDR$Classes.sort, function(x)
  "N" %in% x))
MDRclass$AN <- unlist(lapply(MDR$Classes.sort, function(x)
  "AN" %in% x))
MDRclass$S <- unlist(lapply(MDR$Classes.sort, function(x)
  "S" %in% x))
MDRclass$G <- unlist(lapply(MDR$Classes.sort, function(x)
  "G" %in% x))


#Colors
#A = Aminoglycides alpha("#E69F00", 1)
#AN = Ansamycins alpha("#56B4E9", 1)
#B = Beta-Lactams alpha("#999999", 1)
#F = Fluoroquinolones alpha("#009E73", 1)
#G = Glycopeptides alpha("#882255", 1)
#M = Macrolides alpha("#F0E442", 1)
#N = Nitrofurans alpha("#999933", 1)
#P = Phenicols alpha("#0072B2", 1)
#S = Sulfonamides alpha("#D55E00", 1)
#T = Tetracyclines alpha("#CC79A7", 1)

colors <- c(alpha("#E69F00", 1), #A
            alpha("#56B4E9", 1), #AN
            alpha("#999999", 1), #B
            alpha("#009E73", 1), #F
            alpha("#882255", 1), #G
            alpha("#F0E442", 1), #M
            alpha("#999933", 1), #N
            alpha("#0072B2", 1), #P
            alpha("#D55E00", 1), #S
            alpha("#CC79A7", 1)) #T

labels <- c("Aminoglycides", "Ansamycins", "Beta-Lactams", "Fluoroquinolones", "Glycopeptides", 
            "Macrolides", "Nitrofurans", "Phenicols", "Sulfonamides", "Tetracyclines")

##tabulate number of isolates resistant to each class
MDRclassfreq <- as.data.frame(sapply(MDRclass, function(x) sum(x)))
MDRclassfreq$AMclass <- row.names(MDRclassfreq)
colnames(MDRclassfreq) <- c("value", "AMclass")
MDRclassfreq$pct <- MDRclassfreq$value/nrow(MIC_interp) #note that the denominator here is all isolates, not just isolate tested against that class

MDRclass_fig <- ggplot(MDRclassfreq, aes(factor(AMclass), pct*100, fill = factor(AMclass))) + 
  geom_bar(stat = 'identity') + 
  xlab("Antimicrobial Class") + 
  ylab("Percent Non-susceptible Isolates")+ 
  theme_bw()+
  scale_y_continuous(limits = c(0,100))+
  scale_fill_manual("Class", values=colors, labels = labels)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size=16), legend.title = element_text(size=14), legend.text = element_text(size = 10),
        legend.position = "bottom")
ggsave("Figures and Tables/MDR/MDRclass.png", MDRclass_fig)


rm(AM, Class, colors, labels)

#MDR by year for Prevalence Table1
MDR$Year <- MIC_interp$Year #append Year
nMDR <- MDR %>% group_by(Year) %>% filter(NumClass>=3) %>% tally()
n_yr <- MDR %>% group_by(Year) %>% tally()
nMDR_yr <- full_join(n_yr, nMDR, by="Year", suffix=c("_yr", "_MDR"))
nMDR_yr$prop <- nMDR_yr$n_MDR/nMDR_yr$n_yr

####Cochran-Armitage Tests####

#Utilize antimicrobials used for prevalence analysis that have been filtered based on previous requirements
#excluding those with no bp, < 30 isolates, and not on standard panel (this is MIC_interp data)
#run Cochran Armitage Test for trend over time
CochranArmitage <- sapply(MIC_interp[,a_index], function(x) CochranArmitageTest(table(MIC_interp$Year, x))) %>% t() %>% as.data.frame()

#keep only the antibiotic name P.value, Z.value (statistic)
CochranArmitage$Antibiotic <- row.names(CochranArmitage)
CochranArmitage$statistic <- unlist(CochranArmitage$statistic)
CochranArmitage <- CochranArmitage %>% select(-parameter, -alternative, -method, -data.name)

#reduce significant digits, Benjamini-Hochberg (BH) correction for multiple comparisons
CochranArmitage$p.value <- as.numeric(CochranArmitage$p.value)

#BH tests from the lowest p value to the highest so much arrange the isolates accordingly
CochranArmitage <- CochranArmitage %>% arrange(p.value)
CochranArmitage$BH.p.value <- BH(CochranArmitage$p.value, alpha = 0.05)[["Adjusted.pvalues"]]

#remove the original p-value and Z column, rename columns
CochranArmitage <- CochranArmitage %>% select(-p.value, -statistic) 
CochranArmitage <- dplyr::rename(CochranArmitage, P.Value=BH.p.value)

#changing the number of significant digits
CochranArmitage$P.Value <- signif(CochranArmitage$P.Value, digits=3)

#arrange by alphabetical again
CochranArmitage <- CochranArmitage %>% arrange(Antibiotic)

#add Cochran Armitage test for MDR
CA_mdr <- CochranArmitageTest(matrix(c(nMDR_yr$Year, nMDR_yr$prop),ncol=2))
CochranArmitage <- rbind(CochranArmitage, c("MDR", CA_mdr$p.value))

#For further log-rank and SA testing only want insignificant p values (CA failed to detect trend)
LRABX <- filter(CochranArmitage, P.Value > 0.05)

####Finish Table 1 Prevalence, MDR, Cochran-Armitage####
#first Table1 needs abbreviation
Table1$Abbreviation <- bp$Antimicrobial[match(row.names(Table1), bp$AM.name)]
#Number of Isolates - add to Abbreviation column (currently NA)
Table1$Abbreviation[is.na(Table1$Abbreviation)] <- "Number of Isolates"

#add MDR as row with % (N MDR)
nMDR_yr <- rbind(nMDR_yr, c("Overall", sum(nMDR_yr$n_yr), sum(nMDR_yr$n_MDR), sum(nMDR_yr$n_MDR)/sum(nMDR_yr$n_yr)))
nMDR_yr$prop <- as.numeric(nMDR_yr$prop)
MDR_table1 <- paste0(round(nMDR_yr$prop,2)*100, " (", nMDR_yr$n_MDR,")") #round MDR prevalence, follow with (n MDR isolates)
names(MDR_table1) <- nMDR_yr$Year #name with years
MDR_table1 <- c(MDR_table1, "Abbreviation"="MDR") #add row name
Table1 <- rbind(Table1, MDR_table1)  #add to Table1

#join Table 1 to CA P-value
Table1 <- dplyr::left_join(Table1, CochranArmitage, by=c("Abbreviation"="Antibiotic"))

#move abbreviation to first column
Table1 <- Table1 %>% relocate("Abbreviation")

#round P-value to 3 decimal places
Table1$P.Value <- round(as.numeric(Table1$P.Value),3)
Table1$P.Value <- as.character(Table1$P.Value)
Table1$P.Value[Table1$P.Value=="0"] <- "<0.001"

#nice gt table
Table1_gt <- gt(Table1, rowname_col = "Abbreviation") %>%
  tab_footnote("For each antimicrobial: Non-susceptible isolate prevalence (number of isolates tested)") %>%
  tab_footnote("P-value from Cochran-Armitage test for trend") %>%
  tab_footnote("2007 was excluded because fewer than 30 isolates were available")%>%
  tab_footnote("Multidrug Resistance prevalence (number of MDR isolates)", location=cells_stub("MDR"))

gtsave(Table1_gt, "Figures and Tables/Table1.docx")
gtsave(Table1_gt, "Figures and Tables/Table1.html")


####Breakpoint and Class Table####
names(bp.table)
names(bp.table) <- c("Antimicrobial", "Abbreviation", "S", "I", "R", "Non-susceptible breakpoint", "Species", "Site") #re-name

#add AM classes
AM_Class
bp.table$Class <- AM_Class$Class[match(bp.table$Abbreviation, AM_Class$AM)]

#rename classes (currently one-letter abbreviations)
bp.table$Class <- dplyr::recode(bp.table$Class, "A"="Aminoglycoside", "B"="Beta-Lactam", "P"="Phenicol", "M"="Macrolide/Lincosamide",
                         "T"="Tetracycline", "F"="Fluoroquinolone", "AN"="Ansamycin", "G"="Glycopeptide", "N"="Nitrofuran", "S"="Sulfonamide")

#address mixed case in Body site
bp.table$Species <- str_to_title(bp.table$Species)

#nice gt table
bp_gt <- gt(bp.table) %>%
  cols_align(align="center",
             columns=c(S,I,R))
gtsave(bp_gt, "Figures and Tables/Breakpoint table.docx")
gtsave(bp_gt, "Figures and Tables/Breakpoint table.html")


####Survival analysis data prep####
#Does the MIC data have enough variation for survival analysis?
#Are there an adequate number of isolates with MIC data for analysis?
staph <- staph %>% dplyr::mutate(across(ends_with("MIC"), as.numeric)) #make numeric
summary(staph)

#MIC survival analysis inclusion requirements:
#1) At least 30 isolates (already done)
#2) must have variation in MIC (no more than 90% in any single value)
#3) on standard panel (already done)
AM_cols_MICs <- grep("MIC", colnames(staph)) #get index of MIC columns
print(sapply(staph[,AM_cols_MICs], function(x) prop.table(table(x)))) #tables of MIC by proportion for each antibiotic

#finding AMs that have >90% of isolates with single MIC value
MIC_variation <- sapply(staph[,AM_cols_MICs], function(x) prop.table(table(x)))
MIC_variation <-lapply(MIC_variation, function(x) any(x>0.9))
MICs_to_exclude_from_SA <- unlist(lapply(str_split(names(MIC_variation)[unlist(MIC_variation)],"_"), '[[', 1))
MICs_to_exclude_from_SA

#removing antibiotics from MIC variation analysis
AM_cols_to_exclude <- c(paste(MICs_to_exclude_from_SA, "SIGN", sep="_"), 
                        paste(MICs_to_exclude_from_SA, "MIC", sep="_"),
                        paste(MICs_to_exclude_from_SA, "START", sep="_"),
                        paste(MICs_to_exclude_from_SA, "END", sep="_"))

staph <- staph[ , -which(names(staph) %in% AM_cols_to_exclude)]

#TETRA has only one isolate tested in 2013 and none in 2010-2016; remove the 2013 isolate
staph[(staph$Year==2013 & is.na(staph$TETRA_MIC)==FALSE), ] #one isolate, row 1735
staph[(staph$Year==2013 & is.na(staph$TETRA_SIGN)==FALSE), "TETRA_SIGN"] <- NA
staph[(staph$Year==2013 & is.na(staph$TETRA_START)==FALSE), "TETRA_START"] <- NA
staph[(staph$Year==2013 & is.na(staph$TETRA_END)==FALSE), "TETRA_END"] <- NA
staph[(staph$Year==2013 & is.na(staph$TETRA_MIC)==FALSE), "TETRA_MIC"] <- NA
staph[1735,] #check

#drugs for SA
SA_drugs <- unlist(lapply(str_split(names(staph)[str_detect(names(staph),"_MIC")], "_"), '[[', 1))

write_xlsx(staph,"data/CanineStaphpseud_separated2007-2020.xlsx") #save file ready for survival analysis


####Log-Rank Testing####
#Decided to do log rank tests to test each consecutive year to each other to see if changes from year to year were sig.
#begin with k-way logrank test to see if there is change in the MIC through all years for each antibiotic, then examine years individually

#function to apply LRT to each drug
apply_LRT <- function(drug_name, data) {
  start <- paste0(drug_name, "_START")
  end <- paste0(drug_name, "_END")
  print(drug_name)
  ictest(Surv(get(start),get(end),type="interval2")~Year, scores="logrank2",data=data)
}
LRT <- lapply(as.list(SA_drugs), function(x) apply_LRT(x,staph)) #apply to all drugs in staph
names(LRT) <- SA_drugs

LRT.p <- as.data.frame(cbind(lapply(LRT, "[[", 4))) #get P values as a vector; they are item 4 in each LRT result
LRT.p$Antibiotic <- row.names(LRT.p)
LRT.p <- LRT.p %>% dplyr::rename("P.Value"= "V1")


#BH tests from the lowest p value to the highest so much arrange the isolates accordingly
LRT.p$P.Value <- as.numeric(LRT.p$P.Value)
LRT.p <- LRT.p %>% arrange(P.Value)
LRT.p$BH.p.value <- BH(LRT.p$P.Value, alpha = 0.05)[["Adjusted.pvalues"]]


#consider ranges of years with consistent dilutions (see MIC tables) (2017 - 2020; 2011 - 2016; 2007 - 2010)
#had to set up these LRT tests while referening the MIC tables showing MICs in each year for each drug
#store results
LRT_subsets <- data.frame(AM=character(), Years=character(), P.Value=numeric())
#AMIKAC: still P>0.05
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("AMIKAC", "2007-2016", apply_LRT("AMIKAC", filter(staph, Year<2017))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("AMIKAC", "2017-2020", apply_LRT("AMIKAC", filter(staph, Year>2016))[[4]])

#AMOLCA; P<0.05 for 2011-2016 and 2017-2020
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("AMOCLA", "2007-2010", apply_LRT("AMOCLA", filter(staph, Year<2011))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("AMOCLA", "2011-2016", apply_LRT("AMOCLA", filter(staph, Year>2010 & Year<2017))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("AMOCLA", "2017-2020", apply_LRT("AMOCLA", filter(staph, Year>2016))[[4]])

#AMPICI; P<0.05 for 2011-2016 and 2017-2020
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("AMPICI", "2007-2010", apply_LRT("AMPICI", filter(staph, Year<2011))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("AMPICI", "2011-2016", apply_LRT("AMPICI", filter(staph, Year>2010 & Year<2017))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("AMPICI", "2017-2020", apply_LRT("AMPICI", filter(staph, Year>2016))[[4]])

#CEFAZO; P<0.05 for 2017-2020
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("CEFAZO", "2007-2010", apply_LRT("CEFAZO", filter(staph, Year<2011))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("CEFAZO", "2011-2016", apply_LRT("CEFAZO", filter(staph, Year>2010 & Year<2017))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("CEFAZO", "2017-2020", apply_LRT("CEFAZO", filter(staph, Year>2016))[[4]])

#CEFOVE: p<0.05 for both time periods
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("CEFOVE", "2007-2016", apply_LRT("CEFOVE", filter(staph, Year<2017))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("CEFOVE", "2017-2020", apply_LRT("CEFOVE", filter(staph, Year>2016))[[4]])

#CEFPOD: p<0.05 for both time periods
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("CEFPOD", "2007-2016", apply_LRT("CEFPOD", filter(staph, Year<2017))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("CEFPOD", "2017-2020", apply_LRT("CEFPOD", filter(staph, Year>2016))[[4]])

#CHLORA: P=0.0501 for before 2017; P>0.05 after 2017
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("CHLORA", "2007-2016", apply_LRT("CHLORA", filter(staph, Year<2017))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("CHLORA", "2017-2020", apply_LRT("CHLORA", filter(staph, Year>2016))[[4]])

#CLINDA: P<0.05 for 2010-2020
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("CLINDA", "2007-2009", apply_LRT("CLINDA", filter(staph, Year<2010))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("CLINDA", "2010-2020", apply_LRT("CLINDA", filter(staph, Year>2009))[[4]])

#DOXYCY: P<0.05 2007-2016
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("DOXYCY", "2007-2016", apply_LRT("DOXYCY", filter(staph, Year<2017))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("DOXYCY", "2017-2020", apply_LRT("DOXYCY", filter(staph, Year>2016))[[4]])

#ENROFL: P<0.05 for 2010-2020
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("ENROFL", "2007-2009", apply_LRT("ENROFL", filter(staph, Year<2010))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("ENROFL", "2010-2020", apply_LRT("ENROFL", filter(staph, Year>2009))[[4]])

#ERYTH: P<0.05 2007-2016
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("ERYTH", "2007-2016", apply_LRT("ERYTH", filter(staph, Year<2017))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("ERYTH", "2017-2020", apply_LRT("ERYTH", filter(staph, Year>2016))[[4]])

#GENTAM: P<0.05 for both periods
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("GENTAM", "2007-2016", apply_LRT("GENTAM", filter(staph, Year<2017))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("GENTAM", "2017-2020", apply_LRT("GENTAM", filter(staph, Year>2016))[[4]])

#MARBOF: P<0.05 2007-2016
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("MARBOF", "2007-2016", apply_LRT("MARBOF", filter(staph, Year<2017))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("MARBOF", "2017-2020", apply_LRT("MARBOF", filter(staph, Year>2016))[[4]])

#OXACIL: P<0.05 for 2010-2020
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("OXACIL", "2010-2020", apply_LRT("OXACIL", filter(staph, Year>2009))[[4]])

#TETRA: P>0.05 for all
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("TETRA", "2007-2009", apply_LRT("TETRA", filter(staph, Year<2010))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("TETRA", "2017-2020", apply_LRT("TETRA", filter(staph, Year>2016))[[4]])

#TRISUL: P<0.05 for 2007-2016
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("TRISUL", "2007-2016", apply_LRT("TRISUL", filter(staph, Year<2017))[[4]])
LRT_subsets[nrow(LRT_subsets)+1,] <-  c("TRISUL", "2017-2020", apply_LRT("TRISUL", filter(staph, Year>2016))[[4]])

#BH tests for subsets
LRT_subsets$P.Value <- as.numeric(LRT_subsets$P.Value)
LRT_subsets <- LRT_subsets %>% arrange(P.Value)
LRT_subsets$BH.p.value <- BH(LRT_subsets$P.Value, alpha = 0.05)[["Adjusted.pvalues"]]


#make nice MIC50, MIC90, LRT table
MIC_q_table <- gt(MIC_qs, groupname_col = "AM", rowname_col = "Q") %>%
  fmt_number(columns=contains("20"), drop_trailing_zeros=T) %>% 
  tab_options(row_group.as_column=TRUE) %>%
  tab_footnote("P < 0.05, asymptotic logrank k-sample test (permutation form) across all years, Finkelstein's scores.",
               cells_row_groups(LRT.p$Antibiotic[LRT.p$BH.p.value<0.05])) %>%
  tab_footnote("Logrank tests were not performed on antimicrobials with more than 90% of isolates having a single MIC value",
               cells_row_groups(unique(unlist(MIC_qs$AM[MIC_qs$AM %!in% LRT.p$Antibiotic])))) %>%
  tab_footnote("P < 0.05 Logrank test for year range 2007-2016",
               cells_row_groups(LRT_subsets$AM[LRT_subsets$BH.p.value<0.05 & LRT_subsets$Years=="2007-2016"])) %>%
  tab_footnote("P < 0.05 Logrank test for year range 2011-2016",
               cells_row_groups(LRT_subsets$AM[LRT_subsets$BH.p.value<0.05 & LRT_subsets$Years=="2011-2016"])) %>%
  tab_footnote("P < 0.05 Logrank test for year range 2010-2020",
               cells_row_groups(LRT_subsets$AM[LRT_subsets$BH.p.value<0.05 & LRT_subsets$Years=="2010-2020"])) %>%
  tab_footnote("P < 0.05 Logrank test for year range 2017-2020",
               cells_row_groups(LRT_subsets$AM[LRT_subsets$BH.p.value<0.05 & LRT_subsets$Years=="2017-2020"])) %>%
  tab_footnote("NA: not applicable, no isolates tested against this antimicrobial in this year") %>%
  text_transform(
    locations = cells_stub(rows=everything()),
    fn = function(x){
      sub <- str_split(x, '_', simplify = TRUE)[,2]
      text <- str_split(x, '_', simplify = TRUE)[,1]
      glue::glue("{text}<sub>{sub}</sub>")
    })

gtsave(MIC_q_table, "Figures and Tables/MIC_quantiles.docx")
gtsave(MIC_q_table, "Figures and Tables/MIC_quantiles.html")

####Ritwik and Yuchen MIC analysis####
#select only Year, full date, and MIC columns
MIC_data_req_columns <-
  staph %>% 
  mutate(YMD = as.POSIXlt(paste(paste(Year,Month,Day, sep='-'), Time))) %>%
  select(c('Year', 'YMD', ends_with('_START'), ends_with('_END')))

#create a row/sample id
MIC_data_req_columns$sample_id <- 1:nrow(MIC_data_req_columns)

## Ensure MIC value columns are numeric (some are all NA) 
MIC_data_req_columns <-
  MIC_data_req_columns %>% mutate_at(vars(c(ends_with('START'), ends_with('END'))), 
                                     as.numeric)

## Pivot the data to long form
MIC_data_req_columns <-
  MIC_data_req_columns %>% pivot_longer(
    cols = c(ends_with('START'), ends_with('END')),
    names_to = c('drug', 'start_end'),
    names_pattern = "(.*)_(.*)",
    values_to = "MIC"
  ) %>% pivot_wider(names_from = start_end,
                    values_from = MIC)

#keep only drugs on the standard panel - already done
MIC_data_std_panel_complete <- na.omit(MIC_data_req_columns) #remove NAs

## Run drug specific categorical regression models on year
## ic_np fits Turnbull's estimator, while ic_sp fits a semiparametric
## hazard function estimator based on proportional hazards ('ph') or
## proportional odds ('po') assumption.

fit_drug_specific_model <- function(drug_name, model_type) {
  print(drug_name)
  drug_data <- MIC_data_std_panel_complete %>%
    dplyr::filter(drug == drug_name) %>%
    mutate(Year = as.factor(Year))
  ## Non-parametric fit
  return(ic_np(cbind(START, END) ~ Year,
               data = drug_data))
}

## Fit drug specific models for all drugs
SA_drugs <- unique(MIC_data_std_panel_complete$drug)
drug_specific_models <-
  sapply(as.list(SA_drugs),
         fit_drug_specific_model)
names(drug_specific_models) <- SA_drugs


## Plot drug specific categorical regression model fits
plot_drug_specific_fit <- function(drug_name) {
  print(drug_name)
  drug_data <- MIC_data_std_panel_complete %>%
    dplyr::filter(drug == drug_name)
  png(paste0('Figures and Tables/Nonparametric Fits/',
             drug_name, 
             '.png'))
  plot(
    drug_specific_models[[drug_name]],
    newdata = data.frame(Year = levels(as.factor(drug_data$Year)), 
                         row.names = levels(as.factor(drug_data$Year))),
    xlab = paste(drug_name, 'concentration in ug/ml (C)'),
    ylab = 'P[MIC > C]',
    col = viridisLite::turbo(n = nlevels(as.factor(drug_data$Year)))
  )
  dev.off()
}

lapply(SA_drugs, plot_drug_specific_fit)



# Create additional columns (centered at year 2014; numeric continuous year) for modelling
MIC_data_std_panel_complete$Yeargrp = as.factor(
  as.numeric(as.numeric(MIC_data_std_panel_complete$Year) >= 2014))
MIC_data_std_panel_complete$year_num <- 
  decimal_date(MIC_data_std_panel_complete$YMD)

# Fitting smooth effect of Year by drug. 
##### Create spline regressor columns for full year range
years <- seq(2007, 2021, by = 0.01) #continuous numeric year
year_scaled <- scale(years) #centers years at 2014
year_splines <- bs(year_scaled, degree = 2, df = 5, intercept = T)
spline_df <- data.frame(
  year_num = years,
  year_spl1 = year_splines[, 1],
  year_spl2 = year_splines[, 2],
  year_spl3 = year_splines[, 3],
  year_spl4 = year_splines[, 4],
  year_spl5 = year_splines[, 5]
)
rownames(spline_df) <- spline_df$year_num


## Plot splines
plot(years, year_splines[,1], type = 'l', 
     main = 'Spline regressors as function of time',
     xlab = 'Year',
     ylab = 'Spline regressor',
     col = 1)
lines(years, year_splines[,2], type = 'l', col = 2)
lines(years, year_splines[,3], type = 'l', col = 3)
lines(years, year_splines[,4], type = 'l', col = 4)
lines(years, year_splines[,5], type = 'l', col = 5)

##### Create spline regressor columns for partial year range #####
years2 <- seq(2017, 2021, by = 0.01)
year2_scaled <- scale(years2)
year2_splines <- bs(year2_scaled, degree = 2, df = 5, intercept = T)
spline2_df <- data.frame(
  year_num = years2,
  year_spl1 = year2_splines[, 1],
  year_spl2 = year2_splines[, 2],
  year_spl3 = year2_splines[, 3],
  year_spl4 = year2_splines[, 4],
  year_spl5 = year2_splines[, 5]
)
rownames(spline2_df) <- spline2_df$year_num

## Plot splines
plot(years2, year2_splines[,1], type = 'l', 
     main = 'Spline regressors as function of time',
     xlab = 'Year',
     ylab = 'Spline regressor',
     col = 1)
lines(years2, year2_splines[,2], type = 'l', col = 2)
lines(years2, year2_splines[,3], type = 'l', col = 3)
lines(years2, year2_splines[,4], type = 'l', col = 4)
lines(years2, year2_splines[,5], type = 'l', col = 5)

#CC: examine baselines
baselines_fit <- function(drug_name) {
  print(drug_name)
  ## Subset data to drug
  drug_data = MIC_data_std_panel_complete %>% filter(drug == drug_name)
  
  if(min(drug_data$Year) <= 2010){ #select year range
    drug_data <-
      cbind(drug_data, spline_df[match(round(drug_data$year_num, 2), 
                                       round(spline_df$year_num, 2)),])
  }else{
    drug_data <-
      cbind(drug_data, spline2_df[match(round(drug_data$year_num, 2), 
                                        round(spline2_df$year_num, 2)),])
  }
  png(paste0("Figures and Tables/Baseline Fits/",drug_name,".png"))
  diag_baseline(cbind(START, END) ~ year_num,
                model="ph",
                data=drug_data,
                lgdLocation = "topright")
  title(drug_name)
  dev.off()
  # assign(paste0(drug_name,"_baseline_plot"), plot)
}

#graphics::layout(matrix(seq_len(6),3,1, byrow=T), respect=F)
baseline_plots <- lapply(as.list(SA_drugs), baselines_fit)

#weibull, lnorm, and loglogistic all seem appropriate/equivalent based on visual examination of plots

## Set baseline distribution
baseline_dist <- 'weibull'

#### Function for drug specific fits  ####
get_smooth_drug_specific_fit <- function(drug_name) {
  print(drug_name)
  ## Subset data to drug
  drug_data = MIC_data_std_panel_complete %>% filter(drug == drug_name)
  
  if(min(drug_data$Year) <= 2010){ #select year range
    drug_data <-
      cbind(drug_data, spline_df[match(round(drug_data$year_num, 2), 
                                       round(spline_df$year_num, 2)),])
  }else{
    drug_data <-
      cbind(drug_data, spline2_df[match(round(drug_data$year_num, 2), 
                                        round(spline2_df$year_num, 2)),])
  }
  
  
  ## Prior function: current version is flat priors on scale and shape parameter
  ## of weibull followed by a laplace prior on spline coefficients.
  prior_log <- function(x) { #CC: x is ?
    dexp(abs(x[-c(1,2)]), rate = 0.1) ##CC: this is an exponential prior? not laplace prior?; -c(1,2) ignores the first two parameters (keeps flat)?
  }
  bayesControls(useMLE_start = F)
  bayes_formula <- as.formula(cbind(START, END) ~
                                year_spl1 + year_spl2 + year_spl3 + year_spl4 + year_spl5)
  
  bayes_model <- ic_bayes(
    bayes_formula,
    model = 'aft',
    dist = baseline_dist,
    controls = bayesControls(useMLE_start = F),
    data = drug_data,
    logPriorFxn = prior_log
  )
  
  return(bayes_model)
}

#get model fits
set.seed(1234)
smooth_fits <- lapply(as.list(SA_drugs), get_smooth_drug_specific_fit)
names(smooth_fits) <- SA_drugs

#regression summaries, but note that these are not easily interpretable. Splines must be visualized
lapply(smooth_fits, summary)

##CC: MCMC chains
get_MCMC <- function(drug_name){
  pdf(paste0("Figures and Tables/Smooth_Bayesian_fits/MCMC/",drug_name,".pdf"))
  plot(smooth_fits[[drug_name]]$mcmcList)
  title(drug_name)
  dev.off()
}
lapply(as.list(SA_drugs), get_MCMC)


##CC: models with numeric year
baseline_dist <- 'lnorm'
get_smooth_drug_specific_fit_year_num <- function(drug_name) {
  print(drug_name)
  ## Subset data to drug
  drug_data = MIC_data_std_panel_complete %>% filter(drug == drug_name)
  
  if(min(drug_data$Year) <= 2010){ #select year range
    drug_data <-
      cbind(drug_data, spline_df[match(round(drug_data$year_num, 2), 
                                       round(spline_df$year_num, 2)),])
  }else{
    drug_data <-
      cbind(drug_data, spline2_df[match(round(drug_data$year_num, 2), 
                                        round(spline2_df$year_num, 2)),])
  }
  
  
  ## Prior function: current version is flat priors on scale and shape parameter
  ## of weibull followed by a laplace prior on spline coefficients.
  # prior_log <- function(x) { #CC: x is ?
  #  dexp(abs(x[-c(1,2)]), rate = 0.1) ##CC: this is an exponential prior? not laplace prior?; -c(1,2) ignores the first two parameters (keeps flat)?
  #}
  bayesControls(useMLE_start = F)
  bayes_formula <- as.formula(cbind(START, END) ~
                                year_num)
  
  bayes_model <- ic_bayes(
    bayes_formula,
    model = 'aft',
    dist = baseline_dist,
    controls = bayesControls(useMLE_start = F),
    data = drug_data,
    # logPriorFxn = prior_log
  )
  
  return(bayes_model)
}

#get model fits with year as numeric predictor rather than splines
set.seed(1234)
year_num_fits_ln <- lapply(as.list(SA_drugs), get_smooth_drug_specific_fit_year_num)
names(year_num_fits_ln) <- SA_drugs
#regression model summaries with numeric year predictor
lapply(year_num_fits_ln, summary)


plot_drug_specific_smooth_fit <- function(drug_name){
  print(drug_name)
  bayes_model <- smooth_fits[[drug_name]]
  pred_vars <- c('year_spl1',
                 'year_spl2',
                 'year_spl3',
                 'year_spl4',
                 'year_spl5')
  
  min_year <- min((MIC_data_std_panel_complete %>%
                     filter(drug == drug_name))$year_num)
  
  max_conc <- max((MIC_data_std_panel_complete %>%
                     filter(drug == drug_name))$START) *2
  
  drug_years = floor(min_year):2020
  
  png(paste0('Figures and Tables/Smooth_Bayesian_fits/',
             drug_name,
             '.png'))
  if(min_year <= 2010){
    plot(
      bayes_model,
      newdata = spline_df[spline_df$year_num %in% (drug_years + 0.5),], #plot each line at mid of the year
      xlab = paste0('C = ', drug_name, ' concentration (ug/ml)'),
      ylab = 'P[MIC > C]',
      xlim = c(0,max_conc),
      ylim = c(0,1),
      cis = F,
      col = viridisLite::turbo(n = length(drug_years))
    )
  } else {
    plot(
      bayes_model,
      newdata = spline2_df[spline2_df$year_num %in% (drug_years + 0.5),],
      xlab = paste0('C = ', drug_name, ' concentration (ug/ml)'),
      ylab = 'P[MIC > C]',
      xlim = c(0,max_conc),
      ylim = c(0,1),
      cis = F,
      col = viridisLite::turbo(n = length(drug_years))
    )
  }
  
  dev.off()
  
  
  
  ## Interpretable time effect: plot estimated survival probability 
  ## at different concentration values over year
  drug_years <- seq(floor(min_year),2020, by = 0.25)
  bp_drug <- bp[match(drug_name, bp$Antimicrobial), "NSbp"] 
  min_conc <- min((MIC_data_std_panel_complete %>%
                     filter(drug == drug_name))$END)
  max_conc <- max((MIC_data_std_panel_complete %>%
                     filter(drug == drug_name))$START)
  if (is.na(bp_drug)==TRUE){bp_drug <- floor((min_conc+max_conc)/2)} #for drugs without bp -- CEFTIF
  c1 <- if(min_conc>round(0.25*bp_drug,2)){min_conc} else{c(min_conc, round(0.25*bp_drug,2))}
  c2 <- if(max_conc>round(4*bp_drug,1)){c(round(4*bp_drug,1), max_conc)} else{max_conc}
  concentrations <- sort(unique(c(c1, bp_drug, c2))) #plot bp concentration, 1/4 x bp (if >min), 4xbp (if <max), minimum conc tested, max conc tested
  
  if(min_year <= 2010){
    surv_probs <- survCIs(bayes_model,
                          newdata = spline_df[spline_df$year_num %in% drug_years, 
                                              pred_vars],
                          q = concentrations,
                          ci_level = 0.95)
  } else {
    surv_probs <- survCIs(bayes_model,
                          newdata = spline2_df[spline2_df$year_num %in% drug_years, 
                                               pred_vars],
                          q = concentrations,
                          ci_level = 0.95)
  }
  
  png(paste0('Figures and Tables/Smooth_Bayesian_fits/',
             drug_name,
             '_year_survival.png'))
  
  plot_colors <- viridisLite::viridis(n = length(concentrations))
  
  i = length(concentrations)
  surv_med = 1 - sapply(surv_probs$cis, function(x) as.data.frame(x)[["estimate (median)"]][i])
  surv_low = 1 - sapply(surv_probs$cis, function(x) as.data.frame(x)[["lower"]][i])
  surv_upp = 1 - sapply(surv_probs$cis, function(x) as.data.frame(x)[["upper"]][i])
  
  plot(drug_years, surv_med, 
       type = 'l', lty = 1, 
       ylim = c(0, 1),
       xlab = paste0('Year'),
       ylab = 'P[MIC > C]',
       xaxt = 'n',
       main = paste0('Effect of Year on ', drug_name),
       col = plot_colors[i])
  lines(drug_years, surv_low, lty = 2,
        col = plot_colors[i])
  lines(drug_years, surv_upp, lty = 2,
        col = plot_colors[i])
  for (i in 1:length(concentrations)) {
    surv_med = 1 - sapply(surv_probs$cis, function(x) as.data.frame(x)[["estimate (median)"]][i])
    surv_low = 1 - sapply(surv_probs$cis, function(x) as.data.frame(x)[["lower"]][i])
    surv_upp = 1 - sapply(surv_probs$cis, function(x) as.data.frame(x)[["upper"]][i])
    
    lines(drug_years, surv_med, 
          #type = 'l',
          lty = 1, 
          #ylim = c(0,1),
          #xlab = paste0('Year'),
          #ylab = 'P[MIC > C]',
          #main = paste0('Effect of Year on ', drug_name),
          col = plot_colors[i])
    lines(drug_years, surv_low, lty = 2,
          col = plot_colors[i])
    lines(drug_years, surv_upp, lty = 2,
          col = plot_colors[i])
  }
  
  legend('top',
         lty = rep(c(1,2), each = length(concentrations)),
         col = rep(plot_colors, 2), 
         legend = c(paste0(concentrations, ', median'),
                    paste0(concentrations, ', CI 95%')),
         ncol = 2,
         title = paste0('C=', drug_name, ' conc.(ug/ml)')
  )
  
  x_ticks = seq(floor(min_year), 2020, by = 1)
  axis(side = 1, at = x_ticks)
  dev.off()
}

lapply(as.list(SA_drugs), plot_drug_specific_smooth_fit)
