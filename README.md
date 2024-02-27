---
editor_options: 
  markdown: 
    wrap: 72
---

# RHRV_SDU

------------------------------------------------------------------------

The **RHRV_SDU** extension to the **RHRV** R package allows for the easy
time- and frequency-analysis of heart rate variability (HRV) within
episodes defined by start time and duration. Briefly, by removing the
simultaneous computation of HRV outside the episodes of interest, we
increase processing speed. Moreover, we include the
`CreateTimeAnalysisByEpisodes` function which was so far only mentioned
in an RHRV admin forum.

#### This is package is a modification of the RHRV package --- Heart Rate Variability Analysis of ECG Data. Homepage: <http://rhrv.r-forge.r-project.org/>

We modified RHRV while working on the publication: "*Facilitating
ambulatory heart rate variability analysis using accelerometry-based
classifications of body position and self-reported sleep*" (INSERT DOI
WHEN PUBLISHED)

------------------------------------------------------------------------

## Installation

The extension can be easily installed from GitHub.

```{r setup}
remotes::install_github("marleriee/RHRV_SDU")
library(RHRV)
```

```{r theme, include = FALSE}
theme_set(theme_bw(base_size = 10))
```

------------------------------------------------------------------------

The **RHRV_SDU** extension to the **RHRV** R package allows for the easy
time- and frequency-analysis of heart rate variability (HRV) within
episodes defined by start time and duration. Briefly, by removing the
simultaneous computation of HRV outside the episodes of interest, we
increase processing speed. Moreover, we include the
`CreateTimeAnalysisByEpisodes` function which was so far only mentioned
in an RHRV admin forum.

#### This is package is a modification of the RHRV package --- Heart Rate Variability Analysis of ECG Data. Homepage: <http://rhrv.r-forge.r-project.org/>

We modified RHRV while working on the publication: "*Facilitating
ambulatory heart rate variability analysis using accelerometry-based
classifications of body position and self-reported sleep*" (INSERT DOI
WHEN PUBLISHED)

------------------------------------------------------------------------

## About

The **RHRV_SDU** extension to the **RHRV** package can be employed to
analyze HRV in defined episodes, such as episodes created using a
classification algorithm for accelerometry data in free-living settings.
As long as start time and duration of episodes were defined, e.g., even
via self-report, the **RHRV_SDU** modifications allow for the easy
analysis of within episode HRV. For information regarding the general
use of the RHRV-functions, we recommend the RHRV documentation:
<https://github.com/cran/RHRV>

In the tutorial we present here, we highlight how **RHRV_SDU** may be
used to compute and present ambulatory HRV, specifically RMSSD [ms], in
episodes of accelerometry-classified behavior in three anonymized
individuals participating in the SCREENS randomized controlled trial,
step-by-step.

------------------------------------------------------------------------

### Load required packages

```{r "Load packages", warning=FALSE, message=FALSE}
library(ggplot2)
library(reshape2)
library(dplyr)
library(limma)
library(readxl)
library(ggpubr)
library(tidyr)
library(viridis)
library(hrbrthemes)
library(ggthemes)
library(lubridate)
```

------------------------------------------------------------------------

### Load Required Datasets and Files

#### Datasets

-   This example requires two datasets integrated in the RHRV_SDU
    package:

1.  **Ex_Eps** data containing of the variables:

    | Variable         | Explanation                                     |
    |------------------|-------------------------------------------------|
    | **HRV_time**     | POSIXct variable for the start of the episode   |
    | **HRV_duration** | Episode duration (seconds)                      |
    | **HRV_init**     | Initiation time of episode from .sdf (seconds)  |
    | **HRV_endtime**  | POSIXct variable for the end of the episode     |
    | **HRV_tag**      | Unique identifier for Episode                   |
    | **Subject_ID**   | Unique identifier for Subject                   |
    | **Acc_AID**      | Activity Type                                   |
    | **X**            | Vector containing numerical values for episodes |

2.  **Overview** data containing of the variables:

    | Variable       |                               |
    |----------------|-------------------------------|
    | **Subject_ID** | Unique identifier for Subject |
    | **age**        | Age of Subject                |
    | **sex**        | Biological Sex of Subject     |

The data sets can be easily loaded into R:

```{r "Load Tutorial Data"}
data(INSERT)
head(INSERT)
data(INSERT)
head(INSERT)
```

------------------------------------------------------------------------

#### Files

-   There are sample single lead ECG files added in the package Ex_M1.sdf, Ex_M2.sdf and Ex_F1.sdf
-   The data is stored in the extdata folder of the package
-   The list of files is avaliable using:

system.file("extdata", package = "RHRV") |> list.files()

-   The filename and folde of a specific file can be genereated using:

system.file("extdata", "Ex_F1.sdf", package = "RHRV")

------------------------------------------------------------------------

### Description of Accelerometry Data / Episodes

When analyzing HRV, different HRV outcomes are relevant for different
lengths of recorded episodes. For our analysis, we focus on the main
outcome RMSSD as well as several secondary HRV parameters. In the
Ex_Eps data frame, all episodes in the positions standing, sitting, and
lying were extracted for three individuals. Moreover, we even include a
24-h episode and total sleep.

In the manuscript, we describe the pre-processing of the
accelerometry-determined episodes: "*To qualify for time and frequency
analysis, a minimum episode duration of 360 sec for each bout of
physical behavior was required. The first and last 30 seconds of each
episode were removed to adjust for immediate changes in HRV associated
with hemodynamic changes rather than the behavior itself, and this
ensures a minimum epoch duration of five minutes (300 sec).*" All
episodes in the Ex_Eps dataset have been pre-processed.

------------------------------------------------------------------------

### Using RHRV to compute RHRV within episodes

**To speed up processing and analysis of several HRV episodes for
multiple participants, we employ a for-loop.**

1.  Initially, the location of the ECG or PPG files (e.g. .sdf file) on
    your computer needs to be defined under *datafolder*. Moreover, the
    files describing *episodes-of-interest* as well as *participant
    characteristics* need to be loaded into R as explained above.

To run the following code on your own datasets, you need to define the
location of your HR files, similar to the tutorial. Moreover, the
episode and sample characteristics datasets need to be loaded into R
containing all columns included in the example data.

```{r "Define Datafolder"}
#datafolder <- "/Users/marlenerietz/Library/CloudStorage/OneDrive-KarolinskaInstitutet/HRV/Revision Manuscript/Files for README/"
#"FilePath to Folder to ECG/PPG FILES"
```

```{r "Load Example Data"}
overview <- RHRV::overview
Ex_Eps <- RHRV::Ex_Eps
```

------------------------------------------------------------------------

**Automatic Analysis using For-Loops**

2.  We use individual for-loops for time-analysis and frequency
    analysis, respectively. The for-loop cycles through each participant
    included in the overview file. A file path to the sdf.file with ECG
    data of each participant is created, and it is confirmed that a file
    exists. Next, the *load_HRV* and the *HRV_proc* command combine
    several commands of the original RHVR package into one.

-   **load_HRV** combines *CreateHRVData()*, *SetVerbose()*,
    *LoadBeatSuunto()*
-   **HRV_proc** combines *BuildNIHR()*,

2x *filter_HRV()*, *InterpolateNIHR()*, *CreateTimeAnalysis()*,
*CreateFreqAnalysis()*, and *CalculatePowerBand()* - **filter_HRV** in
the RHRV_SDU extensions filters the ECG data using the parameters
long=50, minbpm=25, maxbpm=(220 - age), last=13.

Then, all episodes for a participant are stored in the dataframe
*allepisodes*. *List1* is defined using the RHRV *AddEpisodes* command.
Then, another for-loop cycles through all episodes, runs time-analysis,
and stores the results for time-analyses in the *results* dataframe. The
results dataframe can be merged with the Ex_Eps dataframe using the
episode Tags.

##### Time Analysis

```{r "Time Analysis", warning=FALSE}

load_HRV = function(title, input) {
  HRVdata <- CreateHRVData()
  HRVdata <- SetVerbose(HRVdata, TRUE)
  HRVdata <- LoadBeatSuunto(HRVdata, input, RecordPath =".")
  
  return(HRVdata)
}

filter_HRV = function (title, age) {
  title = FilterNIHR(title, long=50, minbpm =25, maxbpm=(220 - age), last=13)
  
  return(title)
}

#For Powerband, windows are 300sec (5 min) wide, and frame shift (displacement of window for calculations) is 60 sec. (Updated after Meeting)
HRV_proc = function (title, age) {
  title <- BuildNIHR(title)
  title <- filter_HRV(title, age)
  title <- filter_HRV(title, age)
  title <- InterpolateNIHR (title, freqhr = 4)
  title <- CreateTimeAnalysis(title)
  title <- CreateFreqAnalysis(title)
  title <- CalculatePowerBand(title, indexFreqAnalysis = length(title$FreqAnalysis), size=300, type="fourier", shift =60)
  
}

add_episodes = function (title) {
  title = AddEpisodes(title, InitTimes = c(list_episodes$HRV_init), Tags = c(list_episodes$HRV_tag), Durations = c(list_episodes$HRV_duration), Values = c(list_episodes$X))
  
  return (title)
}


time_analysis_eps = function(title) {
  results_time_sitting <- CreateTimeAnalysisByEpisodes(list1, Tag=c(list_episodes$HRV_tag), size=300, interval=7.8125)
}

overview <- RHRV::overview
Ex_Eps <- RHRV::Ex_Eps

results = data.frame();

datafolder <- system.file("extdata", package = "RHRV") |> list.dirs()

for (i in 1:length(overview))  {
  
  #Creating path for file
  fpath = paste0(datafolder,'/',overview$Subject_ID[i],".sdf")
  
  haveFile = file.exists(fpath)
  print(paste("File path:", fpath))
  print(paste("File exists:", haveFile))
  
  if (haveFile == TRUE) {
    
    #Loading the data
    list1 <- load_HRV(overview$Subject_ID[i], fpath)
    
    list1 <- HRV_proc (list1, overview$age)
    
    print(overview$Subject_ID[i])
    
    
    #Lets get all episodes for the subject
    
    allepisodes <- Ex_Eps[which(Ex_Eps$Subject_ID == overview$Subject_ID[i]),] 
    
    if (nrow(allepisodes)>0) {
      
      print("Adding episodes")
      
      list1 = AddEpisodes(list1, InitTimes = c(allepisodes$HRV_init), Tags = 
                            c(allepisodes$HRV_tag), Durations = c(allepisodes$HRV_duration), Values = c(allepisodes$X),verbose = NULL)
      
      print(allepisodes$HRV_tag)
      
      for (p in 1:length(allepisodes$HRV_tag)) {
        
        results_time <- CreateTimeAnalysisByEpisodes(list1, Tag=allepisodes$HRV_tag[p], size=60, interval=7.8125)
        
        results_time$resultIn$HRV_tag = allepisodes$HRV_tag[p]
        results_time$resultIn$Subject_ID = overview$Subject_ID[i]
        results_time$resultIn$HRV_duration <- allepisodes$HRV_duration[p]
        results_time$resultIn$HRV_time <- allepisodes$HRV_time[p]
        
        results = rbind(results,as.data.frame(results_time$resultIn))
      }
    }
  }
  
}

#merge Ex_Eps file with time-analysis results 
TA_results <- merge (Ex_Eps, results, by=c("HRV_tag", "HRV_duration", "HRV_time", "Subject_ID"))

```

------------------------------------------------------------------------

##### Frequency Analysis

```{r "Frequency Analysis", warning=FALSE}
results = data.frame();

for (i in 1:length(overview$Subject_ID))  {
  
  #Creating path for file
  fpath = paste0(datafolder,overview$Subject_ID[i],".sdf")
  
  haveFile = file.exists(fpath)
  
  if (haveFile == TRUE) {

    #Loading the data
    list1 <- load_HRV(overview$Subject_ID[i], fpath)
    list1 <- HRV_proc (list1, overview$Alder)
    print(overview$Subject_ID[i])
    
    #Lets get all episodes for the subject
  
    allepisodes <- Ex_Eps[which(Ex_Eps$Subject_ID == overview$Subject_ID[i]),] 
    
    if (nrow(allepisodes)>0) {
    
      print("Adding episodes")
      
      list1 = AddEpisodes(list1, InitTimes = c(allepisodes$HRV_init), Tags = 
      c(allepisodes$HRV_tag), Durations = c(allepisodes$HRV_duration), Values = 
      c(allepisodes$X),verbose = NULL)
      
      print(allepisodes$HRV_tag)
       
      for (p in 1:length(allepisodes$HRV_tag)) {
        
        HR <- AnalyzeHRbyEpisodes(list1, 
        Tag=allepisodes$HRV_tag[p],doOutOfEpisodeAnalysis = FALSE, func=mean)
        
        results_freq <-  SplitPowerBandByEpisodes(list1, Tag =allepisodes$HRV_tag[p])
        #doOutOfEpisodeAnalysis = FALSE) @Jan - had to drop this for it to run. Is that         #in RHRV_SDU??
        
        results_freq$resultIn$HRV_tag = allepisodes$HRV_tag[p]
        results_freq$resultIn$Subject_ID = overview$Subject_ID[i]
        results_freq$resultIn$HRV_duration <- allepisodes$HRV_duration[p]
        results_freq$resultIn$HRV_time <- allepisodes$HRV_time[p]
        
        results_freq$resultIn$meanHR <- HR$resultIn
        results_freq$resultIn$mean_ULF <- mean(results_freq$InEpisodes$ULF)
        results_freq$resultIn$sd_ULF <- sd(results_freq$InEpisodes$ULF) 
        results_freq$resultIn$mean_VLF <- mean(results_freq$InEpisodes$VLF)
        results_freq$resultIn$sd_VLF <- sd(results_freq$InEpisodes$VLF) 
        results_freq$resultIn$mean_LF <- mean(results_freq$InEpisodes$LF)
        results_freq$resultIn$sd_LF <- sd(results_freq$InEpisodes$LF) 
        results_freq$resultIn$mean_HF <- mean(results_freq$InEpisodes$HF)
        results_freq$resultIn$sd_HF <- sd(results_freq$InEpisodes$HF) 
        results_freq$resultIn$frames_freq <- length(results_freq$InEpisodes$ULF)
        
        
        results = rbind(results,as.data.frame(results_freq$resultIn))
      }
    }
  }
          
}

#bind results of frequency analysis with time-analysis results
all_results <- merge(TA_results, results, by=c("HRV_tag", "HRV_time", "HRV_duration","Subject_ID"))
```

------------------------------------------------------------------------

**HRV Results have been obtained.**

4.  In the *all_results* dataframe, we have now extracted HRV parameters
    (frequency- and time-analysis) for each subject. We extract all
    outcomes available in *RHRV*, which are: **SDNN, SDANN, SDNNIDX,
    pNN50, SDSD, rMSSD, IRRR, MADRR, TINN, HRVi** for time-analysis &
    **mean (±SD) for ULF, VLF, LF, and HF power** for frequency analysis
    as well as the number of frames analyzed for each episode.

5.  After all results have been obtained, the results may now be cleaned
    based on pre-defined criteria. However, since this example is based
    on selected episodes after cleaning, this is not required in this
    case. For detailed information on our post-processing methodology,
    we refer to our publication.

------------------------------------------------------------------------

### Presenting continuous HRV

**To present continuous HRV, we may graph RMSSD by daytime for each
individual.**

Since the plotting may take a lot of time when done manually, we have
even compiled a for-loop to plot the continuous HRV across episodes for
each individual.

```{r Plotting, warning=FALSE}
subject_ids <- unique(overview$Subject_ID)
plots <- list()
subsets <- list()

for (subject_id in subject_ids) {
  
  # Subset for Subject 
  subset_data <- subset(all_results, Subject_ID %in% subject_id)
  subset_data$Acc_AID <- factor(subset_data$Acc_AID, levels=c("1", "3", "8", "TS", "24"), labels=c("Sitting", "Standing", "Lying", "Total Sleep", "24h"))
  subsets[[subject_id]] <- subset_data
  
  # Extract variables for plotting (Wake-up, Bedtime, 24h and Total Sleep RMSSD)
  wake1 <- as.POSIXct(subset_data$HRV_time[subset_data$Acc_AID == "24h"])
  endday <- as.POSIXct(wake1 + hours(24))
  bedtime <- as.POSIXct(subset_data$HRV_time[subset_data$Acc_AID == "Total Sleep"])
  wake2 <- as.POSIXct(bedtime + seconds(subset_data$HRV_duration[subset_data$Acc_AID == "Total Sleep"]))
  TS_RMSSD <- as.numeric(subset_data$rMSSD[subset_data$Acc_AID == "Total Sleep"])
  Day_RMSSD <- as.numeric(subset_data$rMSSD[subset_data$Acc_AID == "24h"])

  #Extract variables for short episodes
  activity_data <- subset(subset_data, subset_data$Acc_AID %in% c("Sitting", "Standing", "Lying"))
  
  #Plot all 
  plot <- ggplot(activity_data, aes(x = HRV_time, y = rMSSD, col = Acc_AID), size = 2) +
    geom_vline(xintercept = as.POSIXct(wake1), col = "black") +
    annotate("rect", xmin = as.POSIXct(bedtime), xmax = as.POSIXct(wake2), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "lightgrey") +
    annotate("segment", x = bedtime, xend = as.POSIXct(wake2), y = TS_RMSSD, yend = TS_RMSSD, colour = "grey", linetype = "dotted") +
    annotate("segment", x = as.POSIXct(wake1), xend = as.POSIXct(wake1 + hours(24)), y = Day_RMSSD, yend = Day_RMSSD, col = "black", linetype = "dashed") +
    geom_vline(xintercept = as.POSIXct(wake2), linetype = "dotted", col = "grey") +
    geom_vline(xintercept = as.POSIXct(bedtime), linetype = "dotted", col = "grey") +
    geom_point(size = 3.2) +
    labs(x = "Daytime", y = "RMSSD (ms)", title = paste(overview$Subject_ID[overview$Subject_ID == subject_id], ",", overview$age[overview$Subject_ID == subject_id])) +
    geom_linerange(aes(xmin = HRV_time, xmax = (HRV_time + HRV_duration)), linewidth = 0.8) +
    theme_few() +
    scale_x_datetime(date_labels = '%T', 
                     limits = c(as.POSIXct(wake1 - hours(1)), as.POSIXct(wake2 + hours(1))), 
                     breaks = c(as.POSIXct(wake1 - seconds(30)), as.POSIXct(bedtime - seconds(30)), as.POSIXct(wake2 - seconds(30)))) +
    scale_y_continuous(limits = c(20, 70), breaks = c(20, 30, 40, 50, 60, 70)) +
    theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title = element_text(size = 15),  
          legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 10), legend.position = "bottom", 
          legend.background = element_rect(linetype = 1, linewidth = 0.5, colour = 1), plot.title = element_text(size = 15), plot.subtitle = element_text(size = 10)) +
    scale_color_viridis(discrete = T)
  
  plots[[subject_id]] <- plot
}

```

After the for-loop was run, the results can be easily compiled using
`ggarrange`.

```{r "All Plots", warnings=F}
comparison <- ggarrange(plots$Ex_M1, plots$Ex_M2, plots$Ex_F1, nrow = 1, common.legend = TRUE, legend = "bottom")
print(comparison)
```

![Combined Plots - Continuous
HRV](https://github.com/marleriee/RHRV_SDU/blob/master/README_PLOT.png)

------------------------------------------------------------------------

### Conclusion

##### We hope that this tutorial makes it easy to understand how **RHRV** may be used to analyze HRV within episodes, such as episodes classified using self-report or accelerometry. In case of any questions, please do not hesitate to contact us.

\`\`\`
