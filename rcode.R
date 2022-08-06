#### Mate-analysis Acc and RT ####
# Overall effect sizes and subgroup-analysis
# r1=0.5, r2=0
# Date: 18.06.2022
# Studies: 53

#Load packages
library(metafor)
library(meta)
library(ggplot2)
library(dmetar)
library(dplyr)
library(forestplot)

setwd("~/OneDrive - UNSW/Brain_Stimulation/Meta-analysis/4 Data analyses/2. Data Analyses/2022.06.13_53")

#### Acc_Effect sizes calculation ####
acc <- read.csv ("~/OneDrive - UNSW/Brain_Stimulation/Meta-analysis/4 Data analyses/2. Data Analyses/2022.0.09/Datasets/20220510_Raw_Dataset_Acc.csv")

acc <- acc %>%
  mutate(Cognitive.domain=as.factor(Cognitive.domain), Study.design=as.factor(Study.design),
         Targeting.method=as.factor(Targeting.method),
         Frequency=as.factor(Frequency), Ftype=as.factor(Ftype),
         Fsubgroup=as.factor(Fsubgroup), Type.of.sham=as.factor(Type.of.sham))

#parallel trials
acc.parallel <- filter (acc, Study.design == "Parallel")
acc.parallel.es <- escalc(n1i = nACT,n2i = nCTRL, measure = "SMD", vtype = "LS",m1i = mACT_Ch, m2i = mCTRL_Ch,
                          sd1i = sdACT_Ch,sd2i = sdCTRL_Ch,data = acc.parallel)

#cross-over trials
acc.cross <- filter (acc, Study.design == "Cross-over") 
acc.cross.es <- escalc(ni = nACT, measure = "SMCC", vtype = "LS", m1i = mACT_Ch, m2i = mCTRL_Ch,
                       sd1i = sdACT_Ch,sd2i = sdCTRL_Ch, ri = r, data = acc.cross)

# Merge two dataframes together and save the acc_es data #
acc.es <- rbind(acc.parallel.es, acc.cross.es)
# Export the merged dataframe
write.csv(acc.es, file = "Accuracy_changes_es_r")

#### RT_Effect sizes calculation ####
rt <- read.csv("~/OneDrive - UNSW/Brain_Stimulation/Meta-analysis/4 Data analyses/2. Data Analyses/2022.05.09/Datasets/20220510_Raw_Dataset_RT.csv")
rt <- rt %>%
  mutate(Cognitive.domain=as.factor(Cognitive.domain), Study.design=as.factor(Study.design),
         Targeting.method=as.factor(Targeting.method),
         Frequency=as.factor(Frequency), Ftype=as.factor(Ftype),
         Fsubgroup=as.factor(Fsubgroup), Type.of.Sham=as.factor(Type.of.Sham))

#parallel trials
rt.parallel <- filter (rt, Study.design == "Parallel")
rt.parallel.es <- escalc (n1i = nACT,n2i = nCTRL, measure = "SMD", vtype = "UB",m1i = mACT_Ch, m2i = mCTRL_Ch,
                          sd1i = sdACT_Ch,sd2i = sdCTRL_Ch,data = rt.parallel)

#cross-over trials
rt.cross <- filter (rt, Study.design == "Cross-over") 
rt.cross.es <- escalc (ni = nACT, measure = "SMCC", vtype = "LS", m1i = mACT_Ch, m2i = mCTRL_Ch,
                       sd1i = sdACT_Ch,sd2i = sdCTRL_Ch, ri = r, data = rt.cross)

# Merge two dataframes together and save the acc_es data #
rt.es <- rbind(rt.parallel.es, rt.cross.es)
# Export the merged dataframe
write.csv(rt.es, file = "RT_changes_es_r")

#### Acc_Ex ####
acc.ces <- read.csv("~/OneDrive - UNSW/Brain_Stimulation/Meta-analysis/4 Data analyses/2. Data Analyses/2022.06.13_53/Datasets/20220613_CES_Acc.csv")
acc.ces <- acc.ces %>%
  mutate(Cognitive.domain=as.factor(Cognitive.domain), Study.design=as.factor(Study.design),
         Targeting.method=as.factor(Targeting.method),Ftype=as.factor(Ftype),
         Frequency=as.factor(Frequency), Type.of.sham=as.factor(Type.of.sham))

#Accuracy_Excitatory_Overall
acc.ces.ex <- filter (acc.ces, Ftype == "excitatory")
acc_ex <- metagen (TE = acc.ces.ex$cyi, seTE = acc.ces.ex$cvi, data = acc.ces.ex, 
                   studlab = acc.ces.ex$Authors, 
                   complab = acc.ces.ex$Site, 
                   title = acc.ces.ex$Frequency,
                   outclab = acc.ces.ex$Total.pulses.per.session,
                   fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD", 
                   method.tau = "PM", n.e = acc.ces.ex$nACT, n.c = acc.ces.ex$nCTRL)
forest(acc_ex, col.diamond = "red", study.results = TRUE, rightlabs = c("SMD", "95% CI", "Weight"), 
       leftcols = c("studlab","complab", "title", "outclab", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Dosing","Active", "Control"),
       sortvar = TE, just = "center", label.left = "Favors Control", label.right = "Favors Active",
       fontsize = 10, spacing = 0.8, ff.lr = "bold")   

#identify outliers: sig
acc_ex_fo <- find.outliers(acc_ex)

funnel(acc_ex,xlab = "Hedges' g",
       contour = c(0.95, 0.975, 0.99),
       col.contour=c("darkblue","blue","lightblue"))+
  legend(-0.75, 0, c("p < 0.05", "p < 0.025", "p < 0.01"), bty = "n",
         fill=c("darkblue","blue","lightblue"))

acc_ex_egger <- eggers.test(x = acc_ex)

#Cognitive domain
acc_ex_sg_cd <- update.meta (acc_ex, subgroup = acc.ces.ex$Cognitive.domain, random = TRUE, fixed = FALSE)
forest(acc_ex_sg_cd,col.diamond = "red", study.results = TRUE,
       leftcols = c("studlab","complab", "title", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Active", "Control"),
       rightlabs = c("SMD", "95% CI", "Weight"), sortvar = TE, just = "center", col.by = "black",
       label.left = "Favors Control", label.right = "Favors Active", ff.lr = "bold",
       test.subgroup = FALSE, overall = FALSE, subgroup = TRUE, sep.subgroup = ": ", hetstat = TRUE, overall.hetstat = FALSE,
       fontsize = 9, spacing = 0.6)

# Attention
acc.ces.ex.att <- filter (acc.ces.ex, Cognitive.domain == "Attention")
acc_ex_att <- metagen (TE = acc.ces.ex.att$cyi, seTE = acc.ces.ex.att$cvi, data = acc.ces.ex.att, studlab = acc.ces.ex.att$Authors,
                       fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                       method.tau = "PM", n.e = acc.ces.ex.att$nACT, n.c = acc.ces.ex.att$nCTRL)

acc_ex_att_fo <- find.outliers(acc_ex_att)

# Executive function
acc.ces.ex.ef <- filter (acc.ces.ex, Cognitive.domain == "Executive function")
acc_ex_ef <- metagen (TE = acc.ces.ex.ef$cyi, seTE = acc.ces.ex.ef$cvi, data = acc.ces.ex.ef, 
                      studlab = acc.ces.ex.ef$Authors, 
                      complab = acc.ces.ex.ef$Site, 
                      title = acc.ces.ex.ef$Frequency,
                      outclab = acc.ces.ex.ef$Total.pulses.per.session,
                      fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD", 
                      method.tau = "PM", n.e = acc.ces.ex.ef$nACT, n.c = acc.ces.ex.ef$nCTRL)
forest(acc_ex_ef, col.diamond = "red", study.results = TRUE, rightlabs = c("SMD", "95% CI", "Weight"), 
       leftcols = c("studlab","complab", "title", "outclab", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Dosing","Active", "Control"),
       sortvar = TE, just = "center", label.left = "Favors Control", label.right = "Favors Active",
       fontsize = 10, spacing = 1.0, ff.lr = "bold")  

acc_ex_ef_fo <- find.outliers(acc_ex_ef)

# Working memory
acc.ces.ex.wm <- filter (acc.ces.ex, WM == "Working memory")
acc_ex_wm <- metagen (TE = acc.ces.ex.wm$cyi, seTE = acc.ces.ex.wm$cvi, data = acc.ces.ex.wm, studlab = acc.ces.ex.wm$Authors,
                      comb.fixed = FALSE, comb.random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                      method.tau = "PM", n.e = acc.ces.ex.wm$nACT, n.c = acc.ces.ex.wm$nCTRL)

acc_ex_wm_fo <- find.outliers(acc_ex_wm)

# Perception
acc.ces.ex.per <- filter (acc.ces.ex, Cognitive.domain == "Perception")
acc_ex_per <- metagen (TE = acc.ces.ex.per$cyi, seTE = acc.ces.ex.per$cvi, data = acc.ces.ex.per, studlab = acc.ces.ex.per$Authors,
                      fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                      method.tau = "PM", n.e = acc.ces.ex.per$nACT, n.c = acc.ces.ex.per$nCTRL)

acc_ex_per_fo <- find.outliers(acc_ex_per)

# Memory
acc.ces.ex.mem <- filter (acc.ces.ex, Cognitive.domain == "Memory")
acc_ex_mem <- metagen (TE = acc.ces.ex.mem$cyi, seTE = acc.ces.ex.mem$cvi, data = acc.ces.ex.mem, studlab = acc.ces.ex.mem$Authors,
                       fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                       method.tau = "PM", n.e = acc.ces.ex.mem$nACT, n.c = acc.ces.ex.mem$nCTRL)

acc_ex_mem_fo <- find.outliers(acc_ex_mem)

# Motor
acc.ces.ex.mot <- filter (acc.ces.ex, Cognitive.domain == "Motor")
acc_ex_mot <- metagen (TE = acc.ces.ex.mot$cyi, seTE = acc.ces.ex.mot$cvi, data = acc.ces.ex.mot, studlab = acc.ces.ex.mot$Authors,
                       fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                       method.tau = "PM", n.e = acc.ces.ex.mot$nACT, n.c = acc.ces.ex.mot$nCTRL)

acc_ex_mot_fo <- find.outliers(acc_ex_mot)

#### Acc_Subgroup analysis ####
# Frequency_Ftype
acc_ex_sg_fre <- update.meta (acc_ex, subgroup = acc.ces.ex$Frequency, random = TRUE, fixed = FALSE)
forest(acc_ex_sg_fre,col.diamond = "red", study.results = TRUE,
       leftcols = c("studlab","complab", "title", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Active", "Control"),
       rightlabs = c("SMD", "95% CI", "Weight"), sortvar = TE, just = "center", col.by = "black",
       label.left = "Favors Control", label.right = "Favors Active", ff.lr = "bold",
       test.subgroup = FALSE, overall = FALSE, subgroup = TRUE, sep.subgroup = ": ", hetstat = TRUE, overall.hetstat = FALSE,
       fontsize = 9, spacing = 0.6)

# Frequency_Fsubgroup
acc_ex_sg_fs <- update.meta (acc_ex, subgroup = acc.ces.ex$Fsubgroup, random = TRUE, fixed = FALSE)
forest(acc_ex_sg_fs,col.diamond = "red", study.results = TRUE,
       leftcols = c("studlab","complab", "title", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Active", "Control"),
       rightlabs = c("SMD", "95% CI", "Weight"), sortvar = TE, just = "center", col.by = "black",
       label.left = "Favors Control", label.right = "Favors Active", ff.lr = "bold",
       test.subgroup = FALSE, overall = FALSE, subgroup = TRUE, sep.subgroup = ": ", hetstat = TRUE, overall.hetstat = FALSE,
       fontsize = 9, spacing = 0.6)

# Targeting methods
acc_ex_sg_tm <- update.meta (acc_ex, subgroup = acc.ces.ex$Targeting.method, random = TRUE, fixed = FALSE)

# Sham
acc_ex_sg_sham <- update.meta (acc_ex, subgroup = acc.ces.ex$Type.of.sham, random = TRUE, fixed = FALSE)
forest(acc_ex_sg_sham,col.diamond = "red", study.results = TRUE,
       leftcols = c("studlab","complab", "title", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Active", "Control"),
       rightlabs = c("SMD", "95% CI", "Weight"), sortvar = TE, just = "center", col.by = "black",
       label.left = "Favors Control", label.right = "Favors Active", ff.lr = "bold",
       test.subgroup = FALSE, overall = FALSE, subgroup = TRUE, sep.subgroup = ": ", hetstat = TRUE, overall.hetstat = FALSE,
       fontsize = 9, spacing = 0.7)

acc.ces.ex %>%
  group_by(Type.of.sham) %>%
  forestplot(clip = c(-.1, 0.075),
             shapes_gp = fpShapesGp(box = c("blue", "darkred") %>% lapply(function(x) gpar(fill = x, col = "#555555")),
                                    default = gpar(vertices = TRUE)),
             ci.vertices = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             xlab = "Favor Control              Favor Active")

tabletex <- cbind(c("Subgroup","\n",data$Variable), 
                  c("No. of Patients (%)","\n",np), 
                  c("4-Yr Cum. Event Rate\n PCI","\n",data$PCI.Group), 
                  c("4-Yr Cum. Event Rate\n Medical Therapy","\n",data$Medical.Therapy.Group), 
                  c("P Value","\n",data$P.Value))

forestplot(labeltext=tabletext, graph.pos=3, 
           mean=c(NA,NA,data$Point.Estimate), 
           lower=c(NA,NA,data$Low), upper=c(NA,NA,data$High),
           title="Hazard Ratio",
           xlab="     <---PCI Better---    ---Medical Therapy Better--->",
           hrzl_lines=list("3" = gpar(lwd=1, col="#99999922"), 
                           "7" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
                           "15" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
                           "23" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
                           "31" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922")),
           txt_gp=fpTxtGp(label=gpar(cex=1.25),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           col=fpColors(box="black", lines="black", zero = "gray50"),
           zero=1, cex=0.9, lineheight = "auto", boxsize=0.5, colgap=unit(6,"mm"),
           lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.4)

#Stimulation sites
acc_ex_sg_ss <- update.meta (acc_ex,subgroup = acc.ces.ex$SS, random = TRUE, fixed = FALSE)

#Number of sessions
acc_ex_sg_session <- update.meta (acc_ex, subgroup = acc.ces.ex$Type.of.session, random = TRUE, fixed = FALSE)

#### Acc_In #####
acc.ces.in <- filter (acc.ces, Ftype == "inhibitory")
acc_in <- metagen (TE = acc.ces.in$cyi, seTE = acc.ces.in$cvi, data = acc.ces.in, 
                   studlab = acc.ces.in$Authors, 
                   complab = acc.ces.in$Site, 
                   title = acc.ces.in$Frequency,
                   outclab = acc.ces.in$Total.pulses.per.session,
                   fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD", 
                   method.tau = "PM", n.e = acc.ces.in$nACT, n.c = acc.ces.in$nCTRL)

forest(acc_in, col.diamond = "red", study.results = TRUE, rightlabs = c("SMD", "95% CI", "Weight"), 
       leftcols = c("studlab","complab", "title", "outclab", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Dosing","Active", "Control"),
       sortvar = TE, just = "center", label.left = "Favors Control", label.right = "Favors Active",
       fontsize = 10, spacing = 0.8, ff.lr = "bold")

#identify outliers
acc_in_fo <- find.outliers(acc_in)

#Cognitive domain
acc_in_sg_cd <- update.meta (acc_in, subgroup = acc.ces.in$Cognitive.domain, random = TRUE, fixed = FALSE)

forest(acc_in_sg_cd,col.diamond = "red", study.results = TRUE,
       leftcols = c("studlab","complab", "title", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Active", "Control"),
       rightlabs = c("SMD", "95% CI", "Weight"), sortvar = TE, just = "center", col.by = "black",
       label.left = "Favors Control", label.right = "Favors Active", ff.lr = "bold",
       test.subgroup = FALSE, overall = FALSE, subgroup = TRUE, sep.subgroup = ": ", hetstat = TRUE, overall.hetstat = FALSE,
       fontsize = 10, spacing = 0.8)

# Perception
acc.ces.in.per <- filter (acc.ces.in, Cognitive.domain == "Perception")
acc_in_per <- metagen (TE = acc.ces.in.per$cyi, seTE = acc.ces.in.per$cvi, data = acc.ces.in.per, studlab = acc.ces.in.per$Authors,
                       fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                       method.tau = "PM", n.e = acc.ces.in.per$nACT, n.c = acc.ces.in.per$nCTRL)

acc_in_per_fo <- find.outliers(acc_in_per)

# Executive function
acc.ces.in.ef <- filter (acc.ces.in, Cognitive.domain == "Executive function")
acc_in_ef <- metagen (TE = acc.ces.in.ef$cyi, seTE = acc.ces.in.ef$cvi, data = acc.ces.in.ef, studlab = acc.ces.in.ef$Authors,
                      fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                      method.tau = "PM", n.e = acc.ces.in.ef$nACT, n.c = acc.ces.in.ef$nCTRL)

acc_in_ef_fo <- find.outliers(acc_in_ef)

# Attention
acc.ces.in.att <- filter (acc.ces.in, Cognitive.domain == "Attention")
acc_in_att <- metagen (TE = acc.ces.in.att$cyi, seTE = acc.ces.in.att$cvi, data = acc.ces.in.att, studlab = acc.ces.in.att$Authors,
                      fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                      method.tau = "PM", n.e = acc.ces.in.att$nACT, n.c = acc.ces.in.att$nCTRL)

acc_in_att_fo <- find.outliers(acc_in_att)

# Memory
acc.ces.in.mem <- filter (acc.ces.in, Cognitive.domain == "Memory")
acc_in_mem <- metagen (TE = acc.ces.in.mem$cyi, seTE = acc.ces.in.mem$cvi, data = acc.ces.in.mem, studlab = acc.ces.in.mem$Authors,
                       fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                       method.tau = "PM", n.e = acc.ces.in.mem$nACT, n.c = acc.ces.in.mem$nCTRL)

acc_in_mem_fo <- find.outliers(acc_in_mem)

# Motor
acc.ces.in.mot<- filter (acc.ces.in, Cognitive.domain == "Motor")
acc_in_mot <- metagen (TE = acc.ces.in.mot$cyi, seTE = acc.ces.in.mot$cvi, data = acc.ces.in.mot, studlab = acc.ces.in.mot$Authors,
                       fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                       method.tau = "PM", n.e = acc.ces.in.mot$nACT, n.c = acc.ces.in.mot$nCTRL)

acc_in_mot_fo <- find.outliers(acc_in_mot)

#### RT_Ex ####
# Overall
rt.ces <- read.csv ("~/OneDrive - UNSW/Brain_Stimulation/Meta-analysis/4 Data analyses/2. Data Analyses/2022.06.13_53/Datasets/20220613_CES_RT.csv")
rt.ces <- rt.ces %>%
  mutate(Cognitive.domain=as.factor(Cognitive.domain), Study.design=as.factor(Study.design),
         Targeting.method=as.factor(Targeting.method),Ftype=as.factor(Ftype),
         Frequency=as.factor(Frequency), Fsubgroup=as.factor(Fsubgroup), Type.of.sham=as.factor(Type.of.sham))

rt.ces.ex <- filter (rt.ces, Ftype == "excitatory")
rt_ex <- metagen (TE = rt.ces.ex$cyi, seTE = rt.ces.ex$cvi, data = rt.ces.ex, 
                   studlab = rt.ces.ex$Authors, 
                   complab = rt.ces.ex$Site, 
                   title = rt.ces.ex$Frequency,
                   outclab = rt.ces.ex$Total.pulses.per.session,
                   fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD", 
                   method.tau = "PM", n.e = rt.ces.ex$nACT, n.c = rt.ces.ex$nCTRL)

rt_ex_fo <- find.outliers(rt_ex)

forest(rt_ex, col.diamond = "red", study.results = TRUE, rightlabs = c("SMD", "95% CI", "Weight"), 
       leftcols = c("studlab","complab", "title", "outclab", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Dosing","Active", "Control"),
       sortvar = TE, just = "center", label.left = "Favors Active", label.right = "Favors Control",
       fontsize = 10, spacing = 0.8, ff.lr = "bold")   

funnel(rt_ex,xlab = "Hedges' g",
       contour = c(0.95, 0.975, 0.99),
       col.contour=c("darkblue","blue","lightblue"))+
  legend(-1.0, 0, c("p < 0.05", "p < 0.025", "p < 0.01"), bty = "n",
         fill=c("darkblue","blue","lightblue"))

rt_ex_egger <- eggers.test(x = rt_ex)

#### RT_Subgroup analysis ####
# Frequency
rt_ex_sg_fre <- update.meta (rt_ex, subgroup = rt.ces.ex$Frequency, random = TRUE, fixed = FALSE)
forest(rt_ex_sg_fre,col.diamond = "red", study.results = TRUE,
       leftcols = c("studlab","complab", "title", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Active", "Control"),
       rightlabs = c("SMD", "95% CI", "Weight"), sortvar = TE, just = "center", col.by = "black",
       label.left = "Favors Active", label.right = "Favors Control", ff.lr = "bold",
       test.subgroup = FALSE, overall = FALSE, subgroup = TRUE, sep.subgroup = ": ", hetstat = TRUE, overall.hetstat = FALSE,
       fontsize = 9, spacing = 0.6)

# Frequency_Fsubgroup
rt_ex_sg_fs <- update.meta (rt_ex, subgroup= rt.ces.ex$Fsubgroup, random = TRUE, fixed = FALSE)
forest(rt_ex_sg_fs,col.diamond = "red", study.results = TRUE,
       leftcols = c("studlab","complab", "title", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Active", "Control"),
       rightlabs = c("SMD", "95% CI", "Weight"), sortvar = TE, just = "center", col.by = "black",
       label.left = "Favors Active", label.right = "Favors Control", ff.lr = "bold",
       test.subgroup = FALSE, overall = FALSE, subgroup = TRUE, sep.subgroup = ": ", hetstat = TRUE, overall.hetstat = FALSE,
       fontsize = 9, spacing = 0.6)

# Targeting methods
rt_ex_sg_tm <- update.meta (rt_ex, subgroup = rt.ces.ex$Targeting.method, random = TRUE, fixed = FALSE)
forest(rt_ex_sg_tm,col.diamond = "red", study.results = TRUE,
       leftcols = c("studlab","complab", "title", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Active", "Control"),
       rightlabs = c("SMD", "95% CI", "Weight"), sortvar = TE, just = "center", col.by = "black",
       label.left = "Favors Active", label.right = "Favors Control", ff.lr = "bold",
       test.subgroup = FALSE, overall = FALSE, subgroup = TRUE, sep.subgroup = ": ", hetstat = TRUE, overall.hetstat = FALSE,
       fontsize = 9, spacing = 0.6)

# Sham
rt_ex_sg_sham <- update.meta (rt_ex, subgroup = rt.ces.ex$Type.of.sham, random = TRUE, fixed = FALSE)
forest(rt_ex_sg_sham,col.diamond = "red", study.results = TRUE,
       leftcols = c("studlab","complab", "title", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Active", "Control"),
       rightlabs = c("SMD", "95% CI", "Weight"), sortvar = TE, just = "center", col.by = "black",
       label.left = "Favors Active", label.right = "Favors Control", ff.lr = "bold",
       test.subgroup = FALSE, overall = FALSE, subgroup = TRUE, sep.subgroup = ": ", hetstat = TRUE, overall.hetstat = FALSE,
       fontsize = 9, spacing = 0.7)

#Stimulation sites
rt_ex_sg_ss <- update.meta (rt_ex, subgroup = rt.ces.ex$SS, random = TRUE, fixed = FALSE)

#Number of sessions
rt_ex_sg_session <- update.meta (rt_ex, subgroup = rt.ces.ex$Type.of.session, random = TRUE, fixed = FALSE)

# Cognitive domain
rt_ex_sg_cd <- update.meta (rt_ex, subgroup = rt.ces.ex$Cognitive.domain, comb.random = TRUE, comb.fixed = FALSE)
forest(rt_ex_sg_cd,col.diamond = "red", study.results = TRUE,
       leftcols = c("studlab","complab", "title", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Active", "Control"),
       rightlabs = c("SMD", "95% CI", "Weight"), sortvar = TE, just = "center", col.by = "black",
       label.left = "Favors Active", label.right = "Favors Control", ff.lr = "bold",
       test.subgroup = FALSE, overall = FALSE, subgroup = TRUE, sep.subgroup = ": ", hetstat = TRUE, overall.hetstat = FALSE,
       fontsize = 9, spacing = 0.6)

# Attention
rt.ces.ex.att <- filter (rt.ces.ex, Cognitive.domain == "Attention")
rt_ex_att <- metagen (TE = rt.ces.ex.att$cyi, seTE = rt.ces.ex.att$cvi, data = rt.ces.ex.att, studlab = rt.ces.ex.att$Authors,
                      fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                      method.tau = "PM", n.e = rt.ces.ex.att$nACT, n.c = rt.ces.ex.att$nCTRL)

rt_ex_att_fo <- find.outliers(rt_ex_att)

# Executive function
rt.ces.ex.ef <- filter (rt.ces.ex, Cognitive.domain == "Executive function")
rt_ex_ef <- metagen (TE = rt.ces.ex.ef$cyi, seTE = rt.ces.ex.ef$cvi, data = rt.ces.ex.ef, 
                       studlab = rt.ces.ex.ef$Authors, 
                       complab = rt.ces.ex.ef$Site, 
                       title = rt.ces.ex.ef$Frequency,
                       outclab = rt.ces.ex.ef$Total.pulses.per.session,
                       fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD", 
                       method.tau = "PM", n.e = rt.ces.ex.ef$nACT, n.c = rt.ces.ex.ef$nCTRL)

forest(rt_ex_ef, col.diamond = "red", study.results = TRUE, rightlabs = c("SMD", "95% CI", "Weight"), 
       leftcols = c("studlab","complab", "title", "outclab", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Dosing","Active", "Control"),
       sortvar = TE, just = "center", label.left = "Favors Active", label.right = "Favors Control",
       fontsize = 10, spacing = 1.0, ff.lr = "bold")  

rt_ex_ef_fo <- find.outliers(rt_ex_ef)

# Motor
rt.ces.ex.mot <- filter (rt.ces.ex, Cognitive.domain == "Motor")
rt_ex_mot <- metagen (TE = rt.ces.ex.mot$cyi, seTE = rt.ces.ex.mot$cvi, data = rt.ces.ex.mot, 
                     studlab = rt.ces.ex.mot$Authors, 
                     complab = rt.ces.ex.mot$Site, 
                     title = rt.ces.ex.mot$Frequency,
                     outclab = rt.ces.ex.mot$Total.pulses.per.session,
                     fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD", 
                     method.tau = "PM", n.e = rt.ces.ex.mot$nACT, n.c = rt.ces.ex.mot$nCTRL)

forest(rt_ex_mot, col.diamond = "red", study.results = TRUE, rightlabs = c("SMD", "95% CI", "Weight"), 
       leftcols = c("studlab","complab", "title", "outclab", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Dosing","Active", "Control"),
       sortvar = TE, just = "center", label.left = "Favors Active", label.right = "Favors Control",
       fontsize = 10, spacing = 1.0, ff.lr = "bold") 

rt_ex_mot_fo <- find.outliers(rt_ex_mot)

#### RT_In ####
rt.ces.in <- filter (rt.ces, Ftype == "inhibitory")
rt_in <- metagen (TE = rt.ces.in$cyi, seTE = rt.ces.in$cvi, data = rt.ces.in, 
                   studlab = rt.ces.in$Authors, 
                   complab = rt.ces.in$Site, 
                   title = rt.ces.in$Frequency,
                   outclab = rt.ces.in$Total.pulses.per.session,
                   fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD", 
                   method.tau = "PM", n.e = rt.ces.in$nACT, n.c = rt.ces.in$nCTRL)

forest(rt_in, col.diamond = "red", study.results = TRUE, rightlabs = c("SMD", "95% CI", "Weight"), 
       leftcols = c("studlab","complab", "title", "outclab", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Dosing","Active", "Control"),
       sortvar = TE, just = "center", label.left = "Favors Active", label.right = "Favors Control",
       fontsize = 10, spacing = 0.8, ff.lr = "bold")

#identify outliers: non-sig
rt_in_fo <- find.outliers(rt_in)

#Cognitive domain
rt_in_sg_cd <- update.meta (rt_in, subgroup = rt.ces.in$Cognitive.domain, random = TRUE, fixed = FALSE)
forest(rt_in_sg_cd,col.diamond = "red", study.results = TRUE,
       leftcols = c("studlab","complab", "title", "n.e","n.c"), 
       leftlabs = c("Authors", "Site", "Frequency", "Active", "Control"),
       rightlabs = c("SMD", "95% CI", "Weight"), sortvar = TE, just = "center", col.by = "black",
       label.left = "Favors Active", label.right = "Favors Control", ff.lr = "bold",
       test.subgroup = FALSE, overall = FALSE, subgroup = TRUE, sep.subgroup = ": ", hetstat = TRUE, overall.hetstat = FALSE,
       fontsize = 10, spacing = 0.8)

# Attention
rt.ces.in.att <- filter (rt.ces.in, Cognitive.domain == "Attention")
rt_in_att <- metagen (TE = rt.ces.in.att$cyi, seTE = rt.ces.in.att$cvi, data = rt.ces.in.att, studlab = rt.ces.in.att$Authors,
                      fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                      method.tau = "PM", n.e = rt.ces.in.att$nACT, n.c = rt.ces.in.att$nCTRL)

rt_in_att_fo <- find.outliers(rt_in_att)

# Executive function
rt.ces.in.ef <- filter (rt.ces.in, Cognitive.domain == "Executive function")
rt_in_ef <- metagen (TE = rt.ces.in.ef$cyi, seTE = rt.ces.in.ef$cvi, data = rt.ces.in.ef, studlab = rt.ces.in.ef$Authors,
                     fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                     method.tau = "PM", n.e = rt.ces.in.ef$nACT, n.c = rt.ces.in.ef$nCTRL)

rt_in_ef_fo <- find.outliers(rt_in_ef)

# Language
rt.ces.in.lan <- filter (rt.ces.in, Cognitive.domain == "Language")
rt_in_lan <- metagen (TE = rt.ces.in.lan$cyi, seTE = rt.ces.in.lan$cvi, data = rt.ces.in.lan, studlab = rt.ces.in.lan$Authors,
                      fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                      method.tau = "PM", n.e = rt.ces.in.lan$nACT, n.c = rt.ces.in.lan$nCTRL)

rt_in_lan_fo <- find.outliers(rt_in_lan)

# Perception
rt.ces.in.per <- filter (rt.ces.in, Cognitive.domain == "Perception")
rt_in_per <- metagen (TE = rt.ces.in.per$cyi, seTE = rt.ces.in.per$cvi, data = rt.ces.in.per, studlab = rt.ces.in.per$Authors,
                      fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                      method.tau = "PM", n.e = rt.ces.in.per$nACT, n.c = rt.ces.in.per$nCTRL)

rt_in_per_fo <- find.outliers(rt_in_per)

# Working memory
rt.ces.ex.wm <- filter (rt.ces.ex, WM == "Working memory")
rt_ex_wm <- metagen (TE = rt.ces.ex.wm$cyi, seTE = rt.ces.ex.wm$cvi, data = rt.ces.ex.wm, studlab = rt.ces.ex.wm$Authors,
                      fixed = FALSE, random = TRUE, hakn = FALSE, prediction = TRUE,sm ="SMD",
                      method.tau = "PM", n.e = rt.ces.ex.wm$nACT, n.c = rt.ces.ex.wm$nCTRL)

rt_ex_wm_fo <- find.outliers(rt_ex_wm)

#### Summary plots ####
#accuracy
data_acc <- read.csv ("~/OneDrive - UNSW/Brain_Stimulation/Meta-analysis/4 Data analyses/2. Data Analyses/2022.06.13_53/Datasets/summary_acc.csv", stringsAsFactors=FALSE)

## Labels defining subgroups are a little indented!
subgps <- c(2,3,4,5,8,9,12,13,14,15,16,19,20,21)
data_acc$Variable[subgps] <- paste("  ",data_acc$Variable[subgps]) 

## The rest of the columns in the table. 
tabletext <- cbind(c("Subgroup","\n",data_acc$Variable), 
                   c("k","\n",data_acc$K), 
                   c("SMD","\n",data_acc$SMD), 
                   c("τ^2","\n",data_acc$τ.2),
                   c("Q","\n",data_acc$Q),
                   c("I^2","\n",data_acc$I.2),
                   c("p","\n",data_acc$P.Value))

forestplot(labeltext=tabletext,  
           mean=c(NA,NA,data_acc$SMD), graph.pos= 4,
           lower=c(NA,NA,data_acc$Low), upper=c(NA,NA,data_acc$High),
           grid =structure(c(-0.5, 0, 1),
                           gp = gpar(col = "black", lty = 2)), 
           graphwidth = unit(4, "cm"),
           xlab="Favors Control          Favors Active",
           txt_gp=fpTxtGp(label=gpar(cex = 1.2),
                          ticks=gpar(cex = 1.2),
                          xlab=gpar(cex = 1.0),
                          title=gpar(cex = 1.2)),
           col=fpColors(box="royalblue", lines="darkblue", zero = "gray50"),
           zero= 1.000000000001, lwd.zero = 0 , cex= 1, lineheight = "auto", boxsize= 0.3, colgap=unit(3,"mm"),
           lwd.ci= 2, ci.vertices=TRUE, ci.vertices.height = 0.2)

# reaction time
data_rt <- read.csv ("~/OneDrive - UNSW/Brain_Stimulation/Meta-analysis/4 Data analyses/2. Data Analyses/2022.06.13_53/Datasets/summary_rt.csv", stringsAsFactors=FALSE)

## Labels defining subgroups are a little indented!
subgps <- c(2,3,4,5,8,9,12,13,14,15,16,19,20,21)
data_rt$Variable[subgps] <- paste("  ",data_rt$Variable[subgps]) 

## The rest of the columns in the table. 
tabletext <- cbind(c("Subgroup","\n",data_rt$Variable), 
                   c("k","\n",data_rt$K), 
                   c("SMD","\n",data_rt$SMD), 
                   c("τ^2","\n",data_rt$τ.2),
                   c("Q","\n",data_rt$Q),
                   c("I^2","\n",data_rt$I.2),
                   c("p","\n",data_rt$P.Value))


forestplot(labeltext=tabletext,  
           mean=c(NA,NA,data_rt$SMD), graph.pos= 4,
           lower=c(NA,NA,data_rt$Low), upper=c(NA,NA,data_rt$High),
           graphwidth = unit(5, "cm"),
           xlab="                 Favors Active      Favors Control",
           grid =structure(c(-1, 0, 0.5),
                                  gp = gpar(col = "black", lty = 2)), 
           txt_gp=fpTxtGp(label=gpar(cex = 1.2),
                          ticks=gpar(cex = 1.2),
                          xlab=gpar(cex = 1.0),
                          title=gpar(cex = 1.2)),
           col=fpColors(box="royalblue", lines="darkblue", zero = "gray50"), 
           zero= 0.500000000001, lwd.zero = 0, cex= 1, lineheight = "auto", boxsize= 0.3, colgap=unit(6,"mm"),
           lwd.ci= 2, ci.vertices=TRUE, ci.vertices.height = 0.2)

