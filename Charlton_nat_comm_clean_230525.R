
# All necessary packages to load ------------------------------------------
if(!require(readxl)) install.packages("readxl") ;library(readxl)
if(!require(ggplot2)) install.packages("ggplot2") ;library(ggplot2)
if(!require(forcats)) install.packages("forcats"); library(forcats)
if(!require(dplyr)) install.packages("dplyr") ;library(dplyr)
if(!require(arsenal)) install.packages("arsenal") ;library(arsenal)
if(!require(RColorBrewer)) install.packages("RColorBrewer") ;library(RColorBrewer)
if(!require(reshape2)) install.packages("reshape2") ;library(reshape2)
if(!require(ggpubr)) install.packages("ggpubr") ;library(ggpubr)
if(!require(car)) install.packages("car"); data(car, package = "car")
if(!require(tidyverse)) install.packages("tidyverse") ;library(tidyverse)
if(!require(ggprism)) install.packages("ggprism") ;library(ggprism)
if(!require(scales)) install.packages("scales") ;library(scales)
if(!require(cowplot)) install.packages("cowplot"); library(cowplot)
if(!require(magick)) install.packages("magick"); library(magick)
if(!require(lme4)) install.packages("lme4") ;library(lme4)
if(!require(MASS)) install.packages("MASS") ;library(MASS)
if(!require(SciViews)) install.packages("SciViews"); library(SciViews)
if(!require(stats)) install.packages("stats") ;library(stats)
if(!require(rstatix)) install.packages("rstatix"); library(rstatix)
if(!require(emmeans)) install.packages("emmeans"); library(emmeans)
if(!require(nlme)) install.packages("nlme"); library(nlme)
if(!require(glmm)) install.packages("glmm"); library(glmm)
if(!require(cocor)) install.packages("cocor"); library(cocor)
if(!require(ggrepel)) install.packages("ggrepel"); library(ggrepel)
if(!require(plotly)) install.packages("plotly"); library(plotly)
if(!require(ggbreak)) install.packages("ggbreak"); library(ggbreak)
if(!require(gsubfn)) install.packages("gsubfn"); library(gsubfn)  # need 0.7-0 or later
if(!require(DescTools)) install.packages("DescTools"); library(DescTools)
if(!require(car)) install.packages("car"); library(car)
if(!require(lsmeans)) install.packages("lsmeans"); library(lsmeans)
if(!require(multcomp)) install.packages("multcomp"); library(multcomp)
if(!require(plotrix)) install.packages("plotrix") ;library(plotrix)
if(!require(gg.gap)) install.packages("gg.gap"); library(gg.gap)
if(!require(gapminder)) install.packages("gapminder"); data(gapminder, package = "gapminder")
devtools::install_github("sysilviakim/Kmisc")
windowsFonts("Arial" = windowsFont("Arial"))
windowsFonts("Times New Roman" = windowsFont("TT Times New Roman"))
line_size <- 3
base_size <- 70
axis_text_rel_size = -1 
title_text_rel_size = +2
colours3<-c("#FFFFFF","grey42","#d16014","#00798C")
date<-"100625"
setwd("D:/PhD/AGBRESA_LC_ME_CFS/Physical_inactivity_paper/")
output_folder<-(paste("Physical_inactivity_graphs",date,sep="/"))
dir.create(file.path(output_folder), showWarnings = FALSE)


# box_cox_transform -------------------------------------------------------


box_cox_transform<-function (input){
  # input<-dat
  
  #data_lambda<-input %>% select(-Subject, -Session, -Group)
  
  
  
  data_lambda<-input
  lambda<-data.frame(matrix(nrow = length(colnames(data_lambda))))
  rownames(lambda)<-colnames(data_lambda)
  for (i in 1:length(colnames(data_lambda))){
    
    try({ 
      set.seed(123)
      b<-boxcox(lm(data_lambda[,i]~1,na.rm=TRUE))
      set.seed(123)
      lambda[i,] <- b$x[which.max(b$y)]}, silent = TRUE)
  }
  
  
  new_data<-data.frame(matrix(ncol=length(colnames(data_lambda)), nrow = length(rownames(data_lambda))))
  colnames(new_data)<-colnames(data_lambda)
  rownames(new_data)<-rownames(data_lambda)
  ## transformation
  for (i in 1:length(colnames(data_lambda))){
    
    # i=12
    
    try({
      
      if(is.na(lambda[i,]==TRUE)){
        
        lambda_i<-0
      }
      
      else{ lambda_i<-lambda[i,]}
      
      for (j in 1:length(rownames(data_lambda))) {
        
        
        if(lambda_i==0) {
          
          if(is.na(data_lambda[j,i])==TRUE){
            new_data[j,i]<-NA
          }
          
          else if(data_lambda[j,i]==0){
            
            
            new_data[j,i]<-0
          }
          
          else{
            new_data[j,i]<-ln(data_lambda[j,i])
          }
          
          
          
        }
        
        else{
          try({
            new_data[j,i]<-(data_lambda[j,i]^ lambda_i - 1) / lambda_i
          },silent=TRUE)
        }
        
      }
    }, silent=TRUE)
    
    
    
    
  }
  
  
  pvals<-as.data.frame(matrix(nrow=5, ncol=length(colnames(input))))
  colnames(pvals)<-colnames(data_lambda)
  # rownames(pvals)<-c("CON","LC","ME","BDC","HDT55")
  
  
  for (i in 1:length(colnames(input))){
    
    try({
      a<-shapiro.test(input[input$Session=="CON",i])
      pvals[1,i]<-a$p.value
      b<-shapiro.test(input[input$Session=="LC",i])
      pvals[2,i]<-b$p.value
      c<-shapiro.test(input[input$Session=="ME",i])
      pvals[3,i]<-c$p.value
      d<-shapiro.test(input[input$Session=="BDC",i])
      pvals[4,i]<-d$p.value
      e<-shapiro.test(input[input$Session=="HDT55",i])
      pvals[5,i]<-e$p.value
      for (j in 1:length(rownames(input))){
        if(is.na(input[j,i])) {
          new_data[j,i]<-NA
        }
        else if(a$p.value>0.05 && b$p.value>0.05 && c$p.value>0.05 && d$p.value>0.05 && e$p.value>0.05){
          new_data[j,i]<-input[j,i]
        }
        else{
          new_data[j,i]<-new_data[j,i]
        }
      }
      
      
      
    })
    
    
    
    
  }
  
  
  
  
  # new_data[,c("Subject","Session","Group","Sex")]<-input[,c("Subject","Session","Group","Sex")]
  output<-new_data
  return(c(output))
  
}


# normality_test ----------------------------------------------------------


normality_test<-function (input){
  
  # input=data
  
  pvals<-as.data.frame(matrix(nrow=5, ncol=length(colnames(input))))
  colnames(pvals)<-colnames(input)
  rownames(pvals)<-c("CON","LC","ME","BDC","HDT55")
  
  for (i in 4:length(colnames(input))){
    
    try({
      
      a<-shapiro.test(input[input$Session=="CON",i])
      pvals[1,i]<-a$p.value
      b<-shapiro.test(input[input$Session=="LC",i])
      pvals[2,i]<-b$p.value
      c<-shapiro.test(input[input$Session=="ME",i])
      pvals[3,i]<-c$p.value
      d<-shapiro.test(input[input$Session=="BDC",i])
      pvals[4,i]<-d$p.value
      e<-shapiro.test(input[input$Session=="HDT55",i])
      pvals[5,i]<-e$p.value
      
      
      
    },silent=TRUE)
    
    
    
  }
  return(pvals)
}



# Normality test residuals ------------------------------------------------


normality_test_resid<-function (input){
  
  # input=data
  
  pvals<-as.data.frame(matrix(nrow=5, ncol=length(colnames(input))))
  colnames(pvals)<-colnames(input)
  rownames(pvals)<-c("CON","LC","ME","BDC","HDT55")
  
  for (i in 4:length(colnames(input))){
    # i=5
    try({parameter<-colnames(input[i])
    model1<-lm(as.formula(paste(parameter,"~Session",sep="")), data=input[input$Group=="POST-VIRAL",])
    res1<-as.data.frame(resid(model1))
    colnames(res1)<-"Residual"
    res1$row_no<-rownames(res1)
    model2<-lm(as.formula(paste(parameter,"~Session",sep="")), data=input[input$Group=="BED REST",])
    res2<-as.data.frame(resid(model2))
    colnames(res2)<-"Residual"
    res2$row_no<-rownames(res2)
    res_comb<-rbind(res2,res1)
    
    input2<-as.data.frame(matrix(ncol=4, nrow=length(rownames(input))))
    colnames(input2)<-c("Subject","Group","Session","Residual")
    input2$Subject<-input$Subject
    input2$Group<-input$Group
    input2$Session<-input$Session
    input2$Residual<-res_comb[match(rownames(input),res_comb$row_no),"Residual"]},silent=TRUE)
    
    try({
      
      a<-shapiro.test(input2[input2$Session=="CON","Residual"])
      pvals[1,i]<-a$p.value
      b<-shapiro.test(input2[input2$Session=="LC","Residual"])
      pvals[2,i]<-b$p.value
      c<-shapiro.test(input2[input2$Session=="ME","Residual"])
      pvals[3,i]<-c$p.value
      d<-shapiro.test(input2[input2$Session=="BDC","Residual"])
      pvals[4,i]<-d$p.value
      e<-shapiro.test(input2[input2$Session=="HDT55","Residual"])
      pvals[5,i]<-e$p.value
      
      
    },silent=TRUE)
    
    
    
  }
  return(pvals)
}


data<-as.data.frame(read_xlsx("Manuscript_data_clean_230525.xlsx", sheet="Sheet2"))

data<-data %>%
  mutate(Group=dplyr::recode(Group, "AGBRESA"="BED REST","MUSCLE-ME"="POST-VIRAL"))

data$Oxphos_norm<-data$Oxphos/data$SDH

data$E_L_coup<-(data$Uncoupled-data$Leak)/(data$Uncoupled)

data$PN_PNS<-(apply(data[,c("N_linked","ADP")],1, max)/data$Oxphos)*100

data$PS_PNS<-(data$S_linked/data$Uncoupled)*100


data$Leak_norm<-data$Leak/data$SDH
data$membrane_intact<-data$N_linked/data$ADP
data$Percent_IIa_IIx_IIx<-data$Percent_IIa_IIx+data$Percent_IIx

all_pvals<-as.data.frame(matrix(ncol=length(colnames(data)), nrow=7))

colnames(all_pvals)<-colnames(data)
rownames(all_pvals)<-c("BDC-HDT55","BED_REST_test","CON-ME","CON-LC","LC-ME","ANOVA","POST_VIRAL_test")

cor_pvals<-as.data.frame(matrix(ncol=6))
colnames(cor_pvals)<-c("p_value","r","Session","Group","p_value.signif","comparison")


# Check normality ---------------------------------------------------------

op<-par(mfrow=c(1,5), mar=c(2,2,2,2))
sessions<-unique(data$Session)

for (i in 5:length(colnames(data)) ){
  
  
  parameter<-colnames(data[i])
  for (j in 1:length(sessions)){
    try( {session<-sessions[j]
         tmp<- data[data$Session==session,parameter]
         qqnorm(tmp,xlab=parameter, main=paste(session,parameter,sep="_"))
         qqline(tmp)
         recordPlot()
         
         }, silent=T)}
   

  
  
}


norm_pvals<-normality_test(data)
resid_norm_pvals<-normality_test_resid(data)
write.csv(norm_pvals, file=paste(output_folder,"/norm_pvals_",date,".csv", sep = ""))
write.csv(resid_norm_pvals, file=paste(output_folder,"/resid_norm_pvals_",date,".csv", sep = ""))
test<-as.data.frame(box_cox_transform(data))
test[,c("Subject","Session","Group", "Sex")]<-data[,c("Subject","Session","Group","Sex")]
test_norm<-normality_test(test)


# Medians for table 1 -----------------------------------------------------

median(data[data$Session=="BDC","Age"],na.rm=TRUE); iqr(data[data$Session=="BDC","Age"], na.rm=TRUE); quantile(data[data$Session=="BDC","Age"], na.rm=TRUE)
median(data[data$Session=="HDT55","Age"],na.rm=TRUE); IQR(data[data$Session=="HDT55","Age"], na.rm=TRUE);quantile(data[data$Session=="HDT55","Age"], na.rm=TRUE)
median(data[data$Session=="CON","Age"],na.rm=TRUE); IQR(data[data$Session=="CON","Age"], na.rm=TRUE); quantile(data[data$Session=="CON","Age"], na.rm=TRUE)
median(data[data$Session=="LC","Age"],na.rm=TRUE); IQR(data[data$Session=="LC","Age"], na.rm=TRUE); quantile(data[data$Session=="LC","Age"], na.rm=TRUE)
median(data[data$Session=="ME","Age"],na.rm=TRUE); IQR(data[data$Session=="ME","Age"], na.rm=TRUE); quantile(data[data$Session=="ME","Age"], na.rm=TRUE)

median(data[data$Session=="BDC","Height"],na.rm=TRUE); IQR(data[data$Session=="BDC","Height"], na.rm=TRUE); quantile(data[data$Session=="BDC","Height"], na.rm=TRUE)
median(data[data$Session=="HDT55","Height"],na.rm=TRUE); IQR(data[data$Session=="HDT55","Height"], na.rm=TRUE); quantile(data[data$Session=="HDT55","Height"], na.rm=TRUE)
median(data[data$Session=="CON","Height"],na.rm=TRUE); IQR(data[data$Session=="CON","Height"], na.rm=TRUE); quantile(data[data$Session=="CON","Height"], na.rm=TRUE)
median(data[data$Session=="LC","Height"],na.rm=TRUE); IQR(data[data$Session=="LC","Height"], na.rm=TRUE); quantile(data[data$Session=="LC","Height"], na.rm=TRUE)
median(data[data$Session=="ME","Height"],na.rm=TRUE); IQR(data[data$Session=="ME","Height"], na.rm=TRUE); quantile(data[data$Session=="ME","Height"], na.rm=TRUE)

median(data[data$Session=="BDC","Weight"],na.rm=TRUE); IQR(data[data$Session=="BDC","Weight"], na.rm=TRUE); quantile(data[data$Session=="BDC","Weight"], na.rm=TRUE)
median(data[data$Session=="HDT55","Weight"],na.rm=TRUE); IQR(data[data$Session=="HDT55","Weight"], na.rm=TRUE); quantile(data[data$Session=="HDT55","Weight"], na.rm=TRUE)
median(data[data$Session=="CON","Weight"],na.rm=TRUE); IQR(data[data$Session=="CON","Weight"], na.rm=TRUE); quantile(data[data$Session=="CON","Weight"], na.rm=TRUE)
median(data[data$Session=="LC","Weight"],na.rm=TRUE); IQR(data[data$Session=="LC","Weight"], na.rm=TRUE); quantile(data[data$Session=="LC","Weight"], na.rm=TRUE)
median(data[data$Session=="ME","Weight"],na.rm=TRUE); IQR(data[data$Session=="ME","Weight"], na.rm=TRUE); quantile(data[data$Session=="ME","Weight"], na.rm=TRUE)

median(data[data$Session=="BDC","Steps"],na.rm=TRUE); IQR(data[data$Session=="BDC","Steps"], na.rm=TRUE); quantile(data[data$Session=="BDC","Steps"], na.rm=TRUE)
median(data[data$Session=="HDT55","Steps"],na.rm=TRUE); IQR(data[data$Session=="HDT55","Steps"], na.rm=TRUE); quantile(data[data$Session=="HDT55","Steps"], na.rm=TRUE)
median(data[data$Session=="CON","Steps"],na.rm=TRUE); IQR(data[data$Session=="CON","Steps"], na.rm=TRUE); quantile(data[data$Session=="CON","Steps"], na.rm=TRUE)
median(data[data$Session=="LC","Steps"],na.rm=TRUE); IQR(data[data$Session=="LC","Steps"], na.rm=TRUE); quantile(data[data$Session=="LC","Steps"], na.rm=TRUE)
median(data[data$Session=="ME","Steps"],na.rm=TRUE); IQR(data[data$Session=="ME","Steps"], na.rm=TRUE); quantile(data[data$Session=="ME","Steps"], na.rm=TRUE)


median(data[data$Session=="LC","Sx_duration"],na.rm=TRUE); IQR(data[data$Session=="LC","Sx_duration"], na.rm=TRUE); quantile(data[data$Session=="LC","Sx_duration"], na.rm=TRUE)
median(data[data$Session=="ME","Sx_duration"],na.rm=TRUE); IQR(data[data$Session=="ME","Sx_duration"], na.rm=TRUE); quantile(data[data$Session=="ME","Sx_duration"], na.rm=TRUE)


# Age and Sex between cohorts ---------------------------------------------
remove<-test[is.na(test$Age)==TRUE,"Subject"]
test2<-test[!test$Subject %in% remove, ]
test2<-test2[test2$Session!="HDT55",]

kruskal.test(Age~Session, test2)

a<-pairwise.wilcox.test(test2$Age, test2$Session, p.adjust.method="BH")


remove<-test[is.na(test$Sex)==TRUE,"Subject"]
test2<-test[!test$Subject %in% remove, ]
test2<-test2[test2$Session!="HDT55",]

test2$Sex<-ifelse(test2$Sex=="Female",test2$Sex<-1,test2$Sex<-0)

kruskal.test(Sex~Session, test2)

a<-pairwise.wilcox.test(test2$Sex, test2$Session, p.adjust.method="BH")

test2 %>% tukey_hsd(Sex~Session)


remove<-test[is.na(test$Steps)==TRUE,"Subject"]
test2<-test[!test$Subject %in% remove, ]
kruskal.test(Steps~Session, test2)
pairwise.wilcox.test(test2$Steps, test2$Session, p.adjust.method="BH")


remove<-test[is.na(test$Sx_duration)==TRUE,"Subject"]
test2<-test[!test$Subject %in% remove, ]

qqnorm(data[data$Session=="LC","Sx_duration"]);qqline(data[data$Session=="LC","Sx_duration"]) #non normal
qqnorm(data[data$Session=="ME","Sx_duration"]);qqline(data[data$Session=="ME","Sx_duration"]) #non normal

qqnorm(test2[test2$Session=="LC","Sx_duration"]);qqline(test2[test2$Session=="LC","Sx_duration"]) #non normal
qqnorm(test2[test2$Session=="ME","Sx_duration"]);qqline(test2[test2$Session=="ME","Sx_duration"]) #normal

wilcox.test(test2[test2$Session=="LC","Sx_duration"],test2[test2$Session=="ME","Sx_duration"],paired=FALSE )




# Figure 1 ----------------------------------------------------------------
# VO2 absolute ------------------------------------------------------------

# Stats -------------------------------------------------------------------

mean(data[data$Session=="BDC","VO2_abs"],na.rm=TRUE); sd(data[data$Session=="BDC","VO2_abs"],na.rm=TRUE)
mean(data[data$Session=="HDT55","VO2_abs"],na.rm=TRUE);sd(data[data$Session=="HDT55","VO2_abs"],na.rm=TRUE)
mean(data[data$Session=="CON","VO2_abs"],na.rm=TRUE);sd(data[data$Session=="CON","VO2_abs"],na.rm=TRUE)
mean(data[data$Session=="LC","VO2_abs"],na.rm=TRUE);sd(data[data$Session=="LC","VO2_abs"],na.rm=TRUE)
mean(data[data$Session=="ME","VO2_abs"],na.rm=TRUE);sd(data[data$Session=="ME","VO2_abs"],na.rm=TRUE)


remove<-data[is.na(data$VO2_abs)==TRUE,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","VO2_abs"]);qqline(test2[test2$Session=="BDC","VO2_abs"]) 
# qqnorm(test2[test2$Session=="HDT55","VO2_abs"]);qqline(test2[test2$Session=="HDT55","VO2_abs"]) 
# qqnorm(test2[test2$Session=="CON","VO2_abs"]);qqline(test2[test2$Session=="CON","VO2_abs"]) 
# qqnorm(test2[test2$Session=="LC","VO2_abs"]);qqline(test2[test2$Session=="LC","VO2_abs"]) 
# qqnorm(test2[test2$Session=="ME","VO2_abs"]);qqline(test2[test2$Session=="ME","VO2_abs"]) 


a<-kruskal.test(VO2_abs~Session, test2[test2$Group=="POST-VIRAL",])
all_pvals["ANOVA","VO2_abs"]<-a$p.value

a<-pairwise.wilcox.test(test2[test2$Group=="POST-VIRAL","VO2_abs"], test2[test2$Group=="POST-VIRAL","Session"], p.adjust.method="BH")

b<-
  wilcox.test(test2[test2$Session=="BDC","VO2_abs"], test2[test2$Session=="HDT55","VO2_abs"],paired=TRUE)

all_pvals["BDC-HDT55","VO2_abs"]<-b$p.value
all_pvals["CON-ME","VO2_abs"]<-a$p.value[2]
all_pvals["CON-LC","VO2_abs"]<-a$p.value[1]
all_pvals["LC-ME","VO2_abs"]<-a$p.value[4]
all_pvals["BED_REST_test","VO2_abs"]<-"wilcoxon"
all_pvals["POST_VIRAL_test","VO2_abs"]<-"kruskal_wilcoxon_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55","p<0.001","BED REST",
  "CON","ME","p<0.001","POST-VIRAL",
  "CON","LC",paste(format(round(a$p.value[1],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.value[4],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

VO2_abs_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=VO2_abs, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("V\U0307" ~O[2]*"max (L"~min^-1*")"))+
  theme(legend.position = "none",
        # aspect.ratio = 3/1,
        axis.line=element_line(colour="black", size = line_size/2),
        axis.ticks = element_line(colour="black", size = line_size),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size+4, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size+4) / base_size)), 
        axis.title.y = element_text(angle = 90, vjust = 1.5), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( size =rel((base_size+axis_text_rel_size-4)/base_size)),
        strip.background = element_blank(),
        strip.text = element_blank())+
  scale_y_continuous(limits=c(0,5), breaks=seq(0,4,1),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(4), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        # axis.ticks = element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=VO2_abs_a,
       filename = paste(output_folder,"/VO2_abs_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

VO2_abs_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=VO2_abs, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("V\U0307" ~O[2]*"max (L"~min^-1*")"))+
  scale_y_continuous(limits=c(0,5), breaks=seq(0,4,1),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(4.75), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(4.25), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[4,], label = "p.adj",y.position=c(3.75), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        # axis.ticks = element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 15, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=VO2_abs_b,
       filename = paste(output_folder,"/VO2_abs_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

# gas exchange threshold --------------------------------------------------
# Stats -------------------------------------------------------------------

mean(data[data$Session=="BDC","GET"],na.rm=TRUE); sd(data[data$Session=="BDC","GET"],na.rm=TRUE)
mean(data[data$Session=="HDT55","GET"],na.rm=TRUE);sd(data[data$Session=="HDT55","GET"],na.rm=TRUE)
mean(data[data$Session=="CON","GET"],na.rm=TRUE);sd(data[data$Session=="CON","GET"],na.rm=TRUE)
mean(data[data$Session=="LC","GET"],na.rm=TRUE);sd(data[data$Session=="LC","GET"],na.rm=TRUE)
mean(data[data$Session=="ME","GET"],na.rm=TRUE);sd(data[data$Session=="ME","GET"],na.rm=TRUE)

remove<-data[is.na(data$GET)==TRUE,"Subject"]
remove2<-c("P110117","P110118","P110127")
test2<-test[!test$Subject %in% remove & !test$Subject %in% remove2,]

# qqnorm(test2[test2$Session=="BDC","GET"]);qqline(test2[test2$Session=="BDC","GET"])
# qqnorm(test2[test2$Session=="HDT55","GET"]);qqline(test2[test2$Session=="HDT55","GET"])
# qqnorm(test2[test2$Session=="CON","GET"]);qqline(test2[test2$Session=="CON","GET"])
# qqnorm(test2[test2$Session=="LC","GET"]);qqline(test2[test2$Session=="LC","GET"])
# qqnorm(test2[test2$Session=="ME","GET"]);qqline(test2[test2$Session=="ME","GET"])

model <-lme(GET~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
all_pvals["ANOVA","GET"]<-a$`Pr(>Chisq)`
# post-hoc testing 
a<-
  test2[test2$Group=="POST-VIRAL",] %>% tukey_hsd(GET~Session)
        

b<-
  t_test(test2[test2$Group=="BED REST",], GET~Session, paired=TRUE)
b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","GET"]<-b$p
all_pvals["CON-ME","GET"]<-a$p.adj[2]
all_pvals["CON-LC","GET"]<-a$p.adj[1]
all_pvals["LC-ME","GET"]<-a$p.adj[3]
all_pvals["BED_REST_test","GET"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","GET"]<-"anova_tukey_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55","p<0.001","BED REST",
  "CON","ME","p<0.001","POST-VIRAL",
  "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

GET_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=GET, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("GET (L"~min^-1*")"))+
  scale_y_continuous(limits=c(0,3.5), breaks=seq(0,3,1),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(3.1), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        # axis.ticks = element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 40, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=GET_a,
       filename = paste(output_folder,"/GET_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

GET_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=GET, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("GET (L"~min^-1*")"))+
  scale_y_continuous(limits=c(0,3.5), breaks=seq(0,3,1),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(3.3), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(3.1), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[4,], label = "p.adj",y.position=c(2.8), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        # axis.ticks = element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 40, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=GET_b,
       filename = paste(output_folder,"/GET_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")


# minute ventilation ------------------------------------------------------
# Stats -------------------------------------------------------------------

mean(data[data$Session=="BDC","VE_max"],na.rm=TRUE); sd(data[data$Session=="BDC","VE_max"],na.rm=TRUE)
mean(data[data$Session=="HDT55","VE_max"],na.rm=TRUE);sd(data[data$Session=="HDT55","VE_max"],na.rm=TRUE)
mean(data[data$Session=="CON","VE_max"],na.rm=TRUE);sd(data[data$Session=="CON","VE_max"],na.rm=TRUE)
mean(data[data$Session=="LC","VE_max"],na.rm=TRUE);sd(data[data$Session=="LC","VE_max"],na.rm=TRUE)
mean(data[data$Session=="ME","VE_max"],na.rm=TRUE);sd(data[data$Session=="ME","VE_max"],na.rm=TRUE)


remove<-data[is.na(data$VE_max)==TRUE ,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","VE_max"]);qqline(test2[test2$Session=="BDC","VE_max"]) 
# qqnorm(test2[test2$Session=="HDT55","VE_max"]);qqline(test2[test2$Session=="HDT55","VE_max"])
# qqnorm(test2[test2$Session=="CON","VE_max"]);qqline(test2[test2$Session=="CON","VE_max"])
# qqnorm(test2[test2$Session=="LC","VE_max"]);qqline(test2[test2$Session=="LC","VE_max"])
# qqnorm(test2[test2$Session=="ME","VE_max"]);qqline(test2[test2$Session=="ME","VE_max"])

model <-lme(VE_max~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
all_pvals["ANOVA","VE_max"]<-a$`Pr(>Chisq)`


b<-
  t_test(test2[test2$Group=="BED REST",], VE_max~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","VE_max"]<-b$p
# all_pvals["CON-ME","VE_max"]<-a$p.adj[2]
# all_pvals["CON-LC","VE_max"]<-a$p.adj[1]
# all_pvals["LC-ME","VE_max"]<-a$p.adj[3]
all_pvals["BED_REST_test","VE_max"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","VE_max"]<-"anova_tukey_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55","p<0.001","BED REST",
  "CON","ME",paste(round(a$`Pr(>Chisq)`,3)),"POST-VIRAL")
# "CON","LC","p<0.001","POST-VIRAL",
# "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")

plot_data<-data[!data$Subject %in% remove,]

VE_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=VE_max, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("V\U0307" ~E[max]* "(L "*~min^-1*")"))+
  scale_y_continuous(limits=c(0,250), breaks=seq(0,200,50),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(175), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 40, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=VE_a,
       filename = paste(output_folder,"/VE_a",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

VE_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=VE_max, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("V\U0307" ~E[max]* "(L "*~min^-1*")"))+
  scale_y_continuous(limits=c(0,250), breaks=seq(0,200,50),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(230), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 40, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=VE_b,
       filename = paste(output_folder,"/VE_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")


# VE/VCO2 max -------------------------------------------------------------
# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","EqCO2_max"],na.rm=TRUE); sd(data[data$Session=="BDC","EqCO2_max"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","EqCO2_max"],na.rm=TRUE);sd(data[data$Session=="HDT55","EqCO2_max"],na.rm=TRUE)
# mean(data[data$Session=="CON","EqCO2_max"],na.rm=TRUE);sd(data[data$Session=="CON","EqCO2_max"],na.rm=TRUE)
# mean(data[data$Session=="LC","EqCO2_max"],na.rm=TRUE);sd(data[data$Session=="LC","EqCO2_max"],na.rm=TRUE)
# mean(data[data$Session=="ME","EqCO2_max"],na.rm=TRUE);sd(data[data$Session=="ME","EqCO2_max"],na.rm=TRUE)


remove<-data[is.na(data$EqCO2_max)==TRUE,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","EqCO2_max"]);qqline(test2[test2$Session=="BDC","EqCO2_max"])
# qqnorm(test2[test2$Session=="HDT55","EqCO2_max"]);qqline(test2[test2$Session=="HDT55","EqCO2_max"])
# qqnorm(test2[test2$Session=="CON","EqCO2_max"]);qqline(test2[test2$Session=="CON","EqCO2_max"])
# qqnorm(test2[test2$Session=="LC","EqCO2_max"]);qqline(test2[test2$Session=="LC","EqCO2_max"])
# qqnorm(test2[test2$Session=="ME","EqCO2_max"]);qqline(test2[test2$Session=="ME","EqCO2_max"])

model <-lme(EqCO2_max~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
all_pvals["ANOVA","EqCO2_max"]<-a$`Pr(>Chisq)`
      
b<-
  t_test(test2[test2$Group=="BED REST",], EqCO2_max~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","EqCO2_max"]<-b$p
# all_pvals["CON-ME","EqCO2_max"]<-a$p.adj[2]
# all_pvals["CON-LC","EqCO2_max"]<-a$p.adj[1]
# all_pvals["LC-ME","EqCO2_max"]<-a$p.adj[3]
all_pvals["BED_REST_test","EqCO2_max"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","EqCO2_max"]<-"anova_no_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55","p<0.001","BED REST",
  "CON","ME",paste(round(a$`Pr(>Chisq)`,3)),"POST-VIRAL")
# "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
# "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

EqCO2_max_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=EqCO2_max, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("V\U0307"~E[max]*"/V\U0307"~CO[2][max]*""))+
  scale_y_continuous(limits=c(0,70), breaks=seq(0,60,20),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(50), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(15),
              space=c(0.5),
              scales=c(20),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 40, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=EqCO2_max_a,
       filename = paste(output_folder,"/EqCO2_max_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

EqCO2_max_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=EqCO2_max, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("V\U0307"~E[max]*"/V\U0307"~CO[2][max]*""))+
  scale_y_continuous(limits=c(0,70), breaks=seq(0,60,20),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(60), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(15),
              space=c(0.5),
              scales=c(20),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=EqCO2_max_b,
       filename = paste(output_folder,"/EqCO2_max_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")



# Figure 2B ----------------------------------------------------------------

# FCSA --------------------------------------------------------------------

# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","FCSA"],na.rm=TRUE); sd(data[data$Session=="BDC","FCSA"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","FCSA"],na.rm=TRUE);sd(data[data$Session=="HDT55","FCSA"],na.rm=TRUE)
# mean(data[data$Session=="CON","FCSA"],na.rm=TRUE);sd(data[data$Session=="CON","FCSA"],na.rm=TRUE)
# mean(data[data$Session=="LC","FCSA"],na.rm=TRUE);sd(data[data$Session=="LC","FCSA"],na.rm=TRUE)
# mean(data[data$Session=="ME","FCSA"],na.rm=TRUE);sd(data[data$Session=="ME","FCSA"],na.rm=TRUE)


remove<-data[is.na(data$FCSA)==TRUE ,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","FCSA"]);qqline(test2[test2$Session=="BDC","FCSA"])
# qqnorm(test2[test2$Session=="HDT55","FCSA"]);qqline(test2[test2$Session=="HDT55","FCSA"])
# qqnorm(test2[test2$Session=="CON","FCSA"]);qqline(test2[test2$Session=="CON","FCSA"])
# qqnorm(test2[test2$Session=="LC","FCSA"]);qqline(test2[test2$Session=="LC","FCSA"])
# qqnorm(test2[test2$Session=="ME","FCSA"]);qqline(test2[test2$Session=="ME","FCSA"])

model <-lme(FCSA~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
all_pvals["ANOVA","FCSA"]<-a$`Pr(>Chisq)`

b<-
  t_test(test2[test2$Group=="BED REST",], FCSA~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))


all_pvals["BDC-HDT55","FCSA"]<-b$p
# all_pvals["CON-ME","FCSA"]<-a$p.adj[2]
# all_pvals["CON-LC","FCSA"]<-a$p.adj[1]
# all_pvals["LC-ME","FCSA"]<-a$p.adj[3]
all_pvals["BED_REST_test","FCSA"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","FCSA"]<-"anova_tukey_no_hoc"



# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55","p<0.001","BED REST",
  "CON","ME",format(paste(round(a$`Pr(>Chisq)`,3)),nsmall=3),"POST-VIRAL")
# "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
# "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

FCSA_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=FCSA, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("FCSA ("*mu*m^2*")"))+
  scale_y_continuous(limits=c(0,9000), breaks=seq(0,8000,2000),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(7500), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))



ggsave(plot=FCSA_a,
       filename = paste(output_folder,"/FCSA_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

FCSA_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=FCSA, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("FCSA ("*mu*m^2*")"))+
  scale_y_continuous(limits=c(0,9000), breaks=seq(0,8000,2000),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(8300), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=FCSA_b,
       filename = paste(output_folder,"/FCSA_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")


# Figure 3 ----------------------------------------------------------------
# SDH ---------------------------------------------------------------------

# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","SDH"],na.rm=TRUE); sd(data[data$Session=="BDC","SDH"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","SDH"],na.rm=TRUE);sd(data[data$Session=="HDT55","SDH"],na.rm=TRUE)
# mean(data[data$Session=="CON","SDH"],na.rm=TRUE);sd(data[data$Session=="CON","SDH"],na.rm=TRUE)
# mean(data[data$Session=="LC","SDH"],na.rm=TRUE);sd(data[data$Session=="LC","SDH"],na.rm=TRUE)
# mean(data[data$Session=="ME","SDH"],na.rm=TRUE);sd(data[data$Session=="ME","SDH"],na.rm=TRUE)


remove<-data[is.na(data$SDH)==TRUE ,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","SDH"]);qqline(test2[test2$Session=="BDC","SDH"])
# qqnorm(test2[test2$Session=="HDT55","SDH"]);qqline(test2[test2$Session=="HDT55","SDH"])
# qqnorm(test2[test2$Session=="CON","SDH"]);qqline(test2[test2$Session=="CON","SDH"])
# qqnorm(test2[test2$Session=="LC","SDH"]);qqline(test2[test2$Session=="LC","SDH"])
# qqnorm(test2[test2$Session=="ME","SDH"]);qqline(test2[test2$Session=="ME","SDH"])

model <-lme(SDH~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")

a<-car::Anova(model)

all_pvals["ANOVA","SDH"]<-a$`Pr(>Chisq)`
# post-hoc testing
a<-
  test2[test2$Group=="POST-VIRAL",] %>% tukey_hsd(SDH~Session)

b<-
  t_test(test2[test2$Group=="BED REST",], SDH~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))


all_pvals["BDC-HDT55","SDH"]<-b$p
all_pvals["CON-ME","SDH"]<-a$p.adj[2]
all_pvals["CON-LC","SDH"]<-a$p.adj[1]
all_pvals["LC-ME","SDH"]<-a$p.adj[3]
all_pvals["BED_REST_test","SDH"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","SDH"]<-"anova_tukey_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55",paste(format(round(b$p,3), drop0trailing=F)),"BED REST",
  "CON","ME",paste(format(round(a$p.adj[2],3), drop0trailing=F)),"POST-VIRAL",
  "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

SDH_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=SDH*10^5, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(atop("SDH Activity", "(\U0394"~A[660]*~mu*m^-1*""*~s^-1*~10^-5*")")))+
  scale_y_continuous(limits=c(0,1.8), breaks=seq(0,1.5,0.5),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(1.2), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        # axis.ticks = element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=SDH_a,
       filename = paste(output_folder,"/SDH_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

SDH_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=SDH*10^5, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(atop("SDH Activity", "(\U0394"~A[660]*~mu*m^-1*""*~s^-1*~10^-5*")")))+
  scale_y_continuous(limits=c(0,1.8), breaks=seq(0,1.5,0.5),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(1.7), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(1.6), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[4,], label = "p.adj",y.position=c(1.5), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=SDH_b,
       filename = paste(output_folder,"/SDH_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")


# oxphos ------------------------------------------------------------------

# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","Oxphos"],na.rm=TRUE); sd(data[data$Session=="BDC","Oxphos"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","Oxphos"],na.rm=TRUE);sd(data[data$Session=="HDT55","Oxphos"],na.rm=TRUE)
# mean(data[data$Session=="CON","Oxphos"],na.rm=TRUE);sd(data[data$Session=="CON","Oxphos"],na.rm=TRUE)
# mean(data[data$Session=="LC","Oxphos"],na.rm=TRUE);sd(data[data$Session=="LC","Oxphos"],na.rm=TRUE)
# mean(data[data$Session=="ME","Oxphos"],na.rm=TRUE);sd(data[data$Session=="ME","Oxphos"],na.rm=TRUE)


remove<-data[is.na(data$Oxphos)==TRUE | data$membrane_intact>1.1 ,"Subject"]
resp_data<-data[!data$Subject %in% remove, ]
resp_test<-as.data.frame(box_cox_transform(resp_data))
resp_test[,c("Subject","Session","Group", "Sex")]<-resp_data[,c("Subject","Session","Group","Sex")]
resp_test_norm<-normality_test(resp_test)
test2<-resp_test[!resp_test$Subject %in% remove ,]


# qqnorm(test2[test2$Session=="BDC","Oxphos"]);qqline(test2[test2$Session=="BDC","Oxphos"])
# qqnorm(test2[test2$Session=="HDT55","Oxphos"]);qqline(test2[test2$Session=="HDT55","Oxphos"])
# qqnorm(test2[test2$Session=="CON","Oxphos"]);qqline(test2[test2$Session=="CON","Oxphos"])
# qqnorm(test2[test2$Session=="LC","Oxphos"]);qqline(test2[test2$Session=="LC","Oxphos"])
# qqnorm(test2[test2$Session=="ME","Oxphos"]);qqline(test2[test2$Session=="ME","Oxphos"])


a<-kruskal.test(Oxphos~Session, test2[test2$Group=="POST-VIRAL",])
all_pvals["ANOVA","Oxphos"]<-a$p.value

a<-pairwise.wilcox.test(test2[test2$Group=="POST-VIRAL","Oxphos"], test2[test2$Group=="POST-VIRAL","Session"], p.adjust.method="BH")


b<-
  t_test(test2[test2$Group=="BED REST",], Oxphos~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))


all_pvals["BDC-HDT55","Oxphos"]<-b$p
all_pvals["CON-ME","Oxphos"]<-a$p.value[2]
all_pvals["CON-LC","Oxphos"]<-a$p.value[1]
all_pvals["LC-ME","Oxphos"]<-a$p.value[4]
all_pvals["BED_REST_test","Oxphos"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","Oxphos"]<-"kruskal_wilcoxon_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55","p<0.001","BED REST",
  "CON","ME","p<0.001","POST-VIRAL",
  "CON","LC",paste(format(round(a$p.value[1],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME","p<0.001","POST-VIRAL")



plot_data<-data[!data$Subject %in% remove,]

Oxphos_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=Oxphos, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(bquote(atop("Oxidative Phosphorylation","(pmol"*~s^-1*~mg^-1*")")))+
  scale_y_continuous(limits=c(0,180), breaks=seq(0,150,50),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(115), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=Oxphos_a,
       filename = paste(output_folder,"/Oxphos_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

Oxphos_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=Oxphos, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(bquote(atop("Oxidative Phosphorylation","(pmol"*~s^-1*~mg^-1*")")))+
  scale_y_continuous(limits=c(0,180), breaks=seq(0,150,50),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(160), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(150), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[4,], label = "p.adj",y.position=c(130), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=Oxphos_b,
       filename = paste(output_folder,"/Oxphos_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")


# SDH vs VO2max -----------------------------------------------------------

remove<-data[is.na(data$SDH)==TRUE | is.na(data$VO2_rel)==TRUE  ,"Subject"]
plot_data<-data[!data$Subject %in% remove, ]


# df<-list(plot_data[plot_data$Session=="BDC",],plot_data[plot_data$Session=="HDT55",])
# cocor(~VO2_rel+SDH| VO2_rel  + SDH, df )
# 
# df<-list(plot_data[plot_data$Session=="CON",],plot_data[plot_data$Session=="LC",])
# cocor(~VO2_rel+SDH| VO2_rel  + SDH, df )
# 
# df<-list(plot_data[plot_data$Session=="CON",],plot_data[plot_data$Session=="ME",])
# cocor(~VO2_rel+SDH| VO2_rel  + SDH, df )
# 
# df<-list(plot_data[plot_data$Session=="LC",],plot_data[plot_data$Session=="ME",])
# cocor(~VO2_rel+SDH| VO2_rel  + SDH, df )

model<-lm(VO2_rel~SDH*Session, data=plot_data[plot_data$Group=="BED REST",])
car::Anova(model)
model<-lm(VO2_rel~SDH*Session, data=plot_data[plot_data$Group=="POST-VIRAL",])
car::Anova(model)


y<-plot_data$VO2_rel
x<-plot_data$SDH
Session <- plot_data$Session

BDC_cor<-cor.test(x[Session == "BDC"], y[Session == "BDC"], method="pearson")
HDT55_cor<-cor.test(x[Session == "HDT55"], y[Session == "HDT55"], method="pearson")
CON_cor<- cor.test(x[Session == "CON"], y[Session == "CON"], method="pearson")
LC_cor<- cor.test(x[Session == "LC"], y[Session == "LC"], method="pearson")
ME_cor<- cor.test(x[Session == "ME"], y[Session == "ME"], method="pearson")
cor_data_1<-rbind(BDC_cor$p.value,HDT55_cor$p.value,CON_cor$p.value,LC_cor$p.value,ME_cor$p.value)
cor_data_2<-rbind(BDC_cor$estimate,HDT55_cor$estimate,CON_cor$estimate,LC_cor$estimate,ME_cor$estimate)
cor_data_3<-c("BDC","HDT55","CON","LC","ME")
cor_data_4<-c("BED REST","BED REST","POST-VIRAL","POST-VIRAL","POST-VIRAL")
cor_data<-as.data.frame(cbind(cor_data_1,cor_data_2, cor_data_3, cor_data_4))
colnames(cor_data)<-c("p_value","r","Session","Group")
cor_data$p_value<-as.numeric(cor_data$p_value)
cor_data$p_value<-signif(cor_data$p_value,2)
cor_data$r<-as.numeric(cor_data$r)
cor_data$r<-round(cor_data$r,2)
cor_data<-add_significance(cor_data, p.col = "p_value" ,cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
cor_data$comparison<-"SDH_v_VO2"

cor_pvals<-rbind(cor_pvals,cor_data)


# Make graph --------------------------------------------------------------
lab_set_a<-data.frame(x1=c(0.9,0.9),y1=c(65,60),Group=c("BED REST","BED REST"),
                      labs=c(paste("PRE: r=",format(round(cor_data[cor_data$Session=="BDC","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="BDC","p_value"],3),nsmall=3)),
                             paste("POST: r=",format(round(cor_data[cor_data$Session=="HDT55","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="HDT55","p_value"],3),nsmall=3))))

lab_set_b<-data.frame(x1=c(0.9,0.9,0.9),y1=c(65,60,55),Group=c("POST-VIRAL","POST-VIRAL","POST-VIRAL"),
                      labs=c(paste("CON: r=",format(round(cor_data[cor_data$Session=="CON","r"],3),nsmall=2),", p<0.001"),
                             paste("LC: r=",format(round(cor_data[cor_data$Session=="LC","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="LC","p_value"],3),nsmall=3)),
                             paste("ME: r=",format(round(cor_data[cor_data$Session=="ME","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="ME","p_value"],3),nsmall=3))))

VO2_v_SDH_a<-
  ggplot(data=plot_data[plot_data$Group=="BED REST",], aes(y=VO2_rel, x=SDH*10^5))+
  geom_jitter(size=7.5, shape=21,aes(fill=fct_relevel(Session, "BDC","HDT55","CON","LC","ME")),width = 0, height=0, stroke=2)+
  scale_fill_manual(values = c("CON"="white", "LC"=colours3[3], "ME"=colours3[4],"BDC"="white","HDT55"=colours3[2]),
                    labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON","LC"="LC","ME"="ME"))+
  scale_colour_manual(values = c("CON"="white", "LC"=colours3[3], "ME"=colours3[4],"BDC"="white","HDT55"=colours3[2]),
                      labels=c("BDC"="Bedrest-PRE", "HDT55"="POST","CON"="CON","LC"="LC","ME"="ME"))+
  geom_smooth(data=data[data$Session=="HDT55",],method="lm",formula = y~x, se=FALSE, linetype= "solid",fill="white",colour=colours3[2], size=2, fullrange=F)+
  xlab(expression(atop("SDH Activity", "(\U0394"~A[660]*~mu*m^-1*""*~s^-1*~10^-5*")")))+
  ylab(expression("V\U0307" ~O[2][max]*"(mL"~min^-1*~kg^-1*")"))+
  scale_x_continuous(limits=c(0,2.2), breaks=seq(0.5,2.0,0.5), expand=c(0,0))+
  scale_y_continuous(limits = c(0,80), breaks=seq(0,60,20), expand=c(0,0))+
  theme(
    aspect.ratio = 1/1.2,
    legend.position = "none",
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.line=element_line(colour="black", size=line_size/2),
    axis.ticks = element_line(colour="black", size = line_size/2),
    text=element_text(size=base_size, family = "Arial"),
    panel.background  = element_rect(fill="white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      size = rel((title_text_rel_size + base_size) / base_size),
      hjust = 0.5
    ),
    axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
    axis.title.y = element_text(angle = 90, vjust = 1,size = rel((title_text_rel_size + base_size) / base_size),),
    axis.text.y = element_text( size = rel((title_text_rel_size + base_size) / base_size)
    ),
    axis.text.x=element_text(angle=0, size=rel((title_text_rel_size + base_size) / base_size)),
    strip.background = element_blank(),
    strip.text = element_blank())+
  guides(y=guide_axis(cap="upper"))+
  geom_text(data=lab_set_a,aes(y=y1,x=x1, group=Group,label=labs), family="Arial", fontface="bold", size=base_size/4, inherit.aes = TRUE )


ggsave(plot=VO2_v_SDH_a,
       filename = paste(output_folder,"/VO2_v_SDH_a_",date,".png", sep = ""),
       device="png",  width = 20, height = 16, units = "in")



VO2_v_SDH_b<-
  ggplot(data=plot_data[plot_data$Group=="POST-VIRAL",], aes(y=VO2_rel, x=SDH*10^5))+
  geom_jitter(size=7.5, shape=21,aes(fill=fct_relevel(Session, "BDC","HDT55","CON","LC","ME")),width = 0, height=0, stroke=2)+  
  scale_fill_manual(values = c("CON"="white", "LC"=colours3[3], "ME"=colours3[4],"BDC"="white","HDT55"=colours3[2]),
                    labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON","LC"="LC","ME"="ME"))+
  scale_colour_manual(values = c("CON"="white", "LC"=colours3[3], "ME"=colours3[4],"BDC"="white","HDT55"=colours3[2]),
                      labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON","LC"="LC","ME"="ME"))+
  geom_smooth(data=data[data$Session=="CON",],method="lm",formula = y~x, se=FALSE, linetype= "solid",fill="white",colour="black", size=2, fullrange=F)+
  xlab(expression(atop("SDH Activity", "(\U0394"~A[660]*~mu*m^-1*""*~s^-1*~10^-5*")")))+
  ylab(expression("V\U0307" ~O[2][max]*"(mL"~min^-1*~kg^-1*")"))+
  scale_x_continuous(limits=c(0,2.2), breaks=seq(0.5,2.0,0.5), expand=c(0,0))+
  scale_y_continuous(limits = c(0,80), breaks=seq(0,60,20), expand=c(0,0))+
  theme(
    aspect.ratio = 1/1.2,
    legend.position = "none",
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.line=element_line(colour="black", size=line_size/2),
    axis.ticks = element_line(colour="black", size = line_size/2),
    text=element_text(size=base_size, family = "Arial"),
    panel.background  = element_rect(fill="white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      size = rel((title_text_rel_size + base_size) / base_size),
      hjust = 0.5
    ),
    axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
    axis.title.y = element_text(angle = 90, vjust = 1,size = rel((title_text_rel_size + base_size) / base_size),),
    axis.text.y = element_text( size = rel((title_text_rel_size + base_size) / base_size)
    ),
    axis.text.x=element_text(angle=0, size=rel((title_text_rel_size + base_size) / base_size)),
    strip.background = element_blank(),
    strip.text = element_blank())+
  guides(y=guide_axis(cap="upper"))+
  geom_text(data=lab_set_b,aes(y=y1,x=x1, group=Group,label=labs), family="Arial", fontface="bold", size=base_size/4, inherit.aes = TRUE )


ggsave(plot=VO2_v_SDH_b,
       filename = paste(output_folder,"/VO2_v_SDH_b_",date,".png", sep = ""),
       device="png",  width = 20, height = 16, units = "in")


# Oxphos vs VO2max --------------------------------------------------------

remove<-data[is.na(data$Oxphos)==TRUE | is.na(data$VO2_rel)==TRUE | data$membrane_intact>1.1 ,"Subject"]
plot_data<-data[!data$Subject %in% remove, ]

model<-lm(VO2_rel~Oxphos*Session, data=plot_data[plot_data$Group=="BED REST",])
car::Anova(model)
model<-lm(VO2_rel~Oxphos*Session, data=plot_data[plot_data$Group=="POST-VIRAL",])
car::Anova(model)


# df<-list(plot_data[plot_data$Session=="BDC",],plot_data[plot_data$Session=="HDT55",])
# cocor(~VO2_rel+Oxphos| VO2_rel  + Oxphos, df )
# 
# df<-list(plot_data[plot_data$Session=="CON",],plot_data[plot_data$Session=="LC",])
# cocor(~VO2_rel+Oxphos| VO2_rel  + Oxphos, df )
# 
# df<-list(plot_data[plot_data$Session=="CON",],plot_data[plot_data$Session=="ME",])
# cocor(~VO2_rel+Oxphos| VO2_rel  + Oxphos, df )
# 
# df<-list(plot_data[plot_data$Session=="LC",],plot_data[plot_data$Session=="ME",])
# cocor(~VO2_rel+Oxphos| VO2_rel  + Oxphos, df )

y<-plot_data$VO2_rel
x<-plot_data$Oxphos
Session <- plot_data$Session

BDC_cor<-cor.test(x[Session == "BDC"], y[Session == "BDC"], method="pearson")
HDT55_cor<-cor.test(x[Session == "HDT55"], y[Session == "HDT55"], method="pearson")
CON_cor<- cor.test(x[Session == "CON"], y[Session == "CON"], method="pearson")
LC_cor<- cor.test(x[Session == "LC"], y[Session == "LC"], method="pearson")
ME_cor<- cor.test(x[Session == "ME"], y[Session == "ME"], method="pearson")
cor_data_1<-rbind(BDC_cor$p.value,HDT55_cor$p.value,CON_cor$p.value,LC_cor$p.value,ME_cor$p.value)
cor_data_2<-rbind(BDC_cor$estimate,HDT55_cor$estimate,CON_cor$estimate,LC_cor$estimate,ME_cor$estimate)
cor_data_3<-c("BDC","HDT55","CON","LC","ME")
cor_data_4<-c("BED REST","BED REST","POST-VIRAL","POST-VIRAL","POST-VIRAL")
cor_data<-as.data.frame(cbind(cor_data_1,cor_data_2, cor_data_3, cor_data_4))
colnames(cor_data)<-c("p_value","r","Session","Group")
cor_data$p_value<-as.numeric(cor_data$p_value)
cor_data$p_value<-signif(cor_data$p_value,2)
cor_data$r<-as.numeric(cor_data$r)
cor_data$r<-round(cor_data$r,2)
cor_data<-add_significance(cor_data, p.col = "p_value" ,cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
cor_data$comparison<-"Oxphos_v_VO2"

cor_pvals<-rbind(cor_pvals,cor_data)


# Make graph --------------------------------------------------------------


lab_set_a<-data.frame(x1=c(80,80),y1=c(70,65),Group=c("BED REST","BED REST"),
                      labs=c(paste("PRE: r=",format(round(cor_data[cor_data$Session=="BDC","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="BDC","p_value"],3),nsmall=3)),
                             paste("POST: r=",format(round(cor_data[cor_data$Session=="HDT55","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="HDT55","p_value"],3),nsmall=3))))


lab_set_b<-data.frame(x1=c(80,80,80),y1=c(70,65,60),Group=c("POST-VIRAL","POST-VIRAL","POST-VIRAL"),
                      labs=c( paste("CON: r=",format(round(cor_data[cor_data$Session=="CON","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="CON","p_value"],3),nsmall=3)),
                              paste("LC: r=",format(round(cor_data[cor_data$Session=="LC","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="LC","p_value"],3),nsmall=3)),
                              paste("ME: r=",format(round(cor_data[cor_data$Session=="ME","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="ME","p_value"],3),nsmall=3))))

VO2_v_Oxphos_a<-
  ggplot(data=plot_data[plot_data$Group=="BED REST",], aes(y=VO2_rel, x=Oxphos))+
  geom_jitter(size=7.5, shape=21,aes(fill=fct_relevel(Session, "BDC","HDT55","CON","LC","ME")),width = 0, height=0, stroke=2)+ 
  scale_fill_manual(values = c("CON"="white", "LC"=colours3[3], "ME"=colours3[4],"BDC"="white","HDT55"=colours3[2]),
                    labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON","LC"="LC","ME"="ME"))+
  scale_colour_manual(values = c("CON"="white", "LC"=colours3[3], "ME"=colours3[4],"BDC"="white","HDT55"=colours3[2]),
                      labels=c("BDC"="Bedrest-PRE", "HDT55"="POST","CON"="CON","LC"="LC","ME"="ME"))+
  geom_smooth(data=plot_data[plot_data$Session=="BDC",],method="lm",formula = y~x, se=FALSE, linetype= "solid",colour="black", size=2, fullrange=F)+
  geom_smooth(data=plot_data[plot_data$Session=="HDT55",],method="lm",formula = y~x, se=FALSE, linetype= "solid",fill="white",colour=colours3[2], size=2, fullrange=F)+
  xlab(bquote(atop("Oxidative Phosphorylation","(pmol"*~s^-1*~mg^-1*")")))+
  ylab(expression("V\U0307" ~O[2][max]*"(mL"~min^-1*~kg^-1*")"))+
  scale_y_continuous(limits = c(0,75), breaks=seq(0,60,20), expand=c(0,0))+
  scale_x_continuous(limits=c(0,160), breaks=seq(30,150,30), expand=c(0,0))+
  theme(
    aspect.ratio = 1/1.2,
    legend.position = "none",
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.line=element_line(colour="black", size=line_size/2),
    axis.ticks = element_line(colour="black", size = line_size/2),
    text=element_text(size=base_size, family = "Arial"),
    panel.background  = element_rect(fill="white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      size = rel((title_text_rel_size + base_size) / base_size),
      hjust = 0.5
    ),
    axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
    axis.title.y = element_text(angle = 90, vjust = 1,size = rel((title_text_rel_size + base_size) / base_size),),
    axis.text.y = element_text( size = rel((title_text_rel_size + base_size) / base_size)
    ),
    axis.text.x=element_text(angle=0, size=rel((title_text_rel_size + base_size) / base_size)),
    strip.background = element_blank(),
    strip.text = element_blank())+
  guides(y=guide_axis(cap="upper"))+
  geom_text(data=lab_set_a,aes(y=y1,x=x1, group=Group,label=labs), family="Arial", fontface="bold", size=base_size/4, inherit.aes = TRUE )



ggsave(plot=VO2_v_Oxphos_a,
       filename = paste(output_folder,"/VO2_v_Oxphos_a_",date,".png", sep = ""),
       device="png",  width = 20, height = 16, units = "in")


VO2_v_Oxphos_b<-
  ggplot(data=plot_data[plot_data$Group=="POST-VIRAL",], aes(y=VO2_rel, x=Oxphos))+
  geom_jitter(size=7.5, shape=21,aes(fill=fct_relevel(Session, "BDC","HDT55","CON","LC","ME")),width = 0, height=0, stroke=2)+  
  scale_fill_manual(values = c("CON"="white", "LC"=colours3[3], "ME"=colours3[4],"BDC"="white","HDT55"=colours3[2]),
                    labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON","LC"="LC","ME"="ME"))+
  scale_colour_manual(values = c("CON"="white", "LC"=colours3[3], "ME"=colours3[4],"BDC"="white","HDT55"=colours3[2]),
                      labels=c("BDC"="Bedrest-PRE", "HDT55"="POST","CON"="CON","LC"="LC","ME"="ME"))+
  geom_smooth(data=plot_data[plot_data$Session=="CON",],method="lm",formula = y~x, se=FALSE, linetype= "solid",colour="black", size=2, fullrange=F)+
  xlab(bquote(atop("Oxidative Phosphorylation","(pmol"*~s^-1*~mg^-1*")")))+
  ylab(expression("V\U0307" ~O[2][max]*"(mL"~min^-1*~kg^-1*")"))+
  scale_y_continuous(limits = c(0,75), breaks=seq(0,60,20), expand=c(0,0))+
  scale_x_continuous(limits=c(0,160), breaks=seq(30,150,30), expand=c(0,0))+
  theme(
    aspect.ratio = 1/1.2,
    legend.position = "none",
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.line=element_line(colour="black", size=line_size/2),
    axis.ticks = element_line(colour="black", size = line_size/2),
    text=element_text(size=base_size, family = "Arial"),
    panel.background  = element_rect(fill="white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      size = rel((title_text_rel_size + base_size) / base_size),
      hjust = 0.5
    ),
    axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
    axis.title.y = element_text(angle = 90, vjust = 1,size = rel((title_text_rel_size + base_size) / base_size),),
    axis.text.y = element_text( size = rel((title_text_rel_size + base_size) / base_size)
    ),
    axis.text.x=element_text(angle=0, size=rel((title_text_rel_size + base_size) / base_size)),
    strip.background = element_blank(),
    strip.text = element_blank())+
  guides(y=guide_axis(cap="upper"))+
  geom_text(data=lab_set_b,aes(y=y1,x=x1, group=Group,label=labs), family="Arial", fontface="bold", size=base_size/4, inherit.aes = TRUE )


ggsave(plot=VO2_v_Oxphos_b,
       filename = paste(output_folder,"/VO2_v_Oxphos_b_",date,".png", sep = ""),
       device="png",  width = 20, height = 16, units = "in")


# Figure 4 ----------------------------------------------------------------

# Capillary fibre ratio ---------------------------------------------------

# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","CF"],na.rm=TRUE); sd(data[data$Session=="BDC","CF"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","CF"],na.rm=TRUE);sd(data[data$Session=="HDT55","CF"],na.rm=TRUE)
# mean(data[data$Session=="CON","CF"],na.rm=TRUE);sd(data[data$Session=="CON","CF"],na.rm=TRUE)
# mean(data[data$Session=="LC","CF"],na.rm=TRUE);sd(data[data$Session=="LC","CF"],na.rm=TRUE)
# mean(data[data$Session=="ME","CF"],na.rm=TRUE);sd(data[data$Session=="ME","CF"],na.rm=TRUE)


remove<-data[is.na(data$CF)==TRUE ,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","CF"]);qqline(test2[test2$Session=="BDC","CF"])
# qqnorm(test2[test2$Session=="HDT55","CF"]);qqline(test2[test2$Session=="HDT55","CF"])
# qqnorm(test2[test2$Session=="CON","CF"]);qqline(test2[test2$Session=="CON","CF"])
# qqnorm(test2[test2$Session=="LC","CF"]);qqline(test2[test2$Session=="LC","CF"])
# qqnorm(test2[test2$Session=="ME","CF"]);qqline(test2[test2$Session=="ME","CF"])

model <-lme(CF~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
all_pvals["ANOVA","CF"]<-a$`Pr(>Chisq)`

# post-hoc testing
a<-
  test2[test2$Group=="POST-VIRAL",] %>% tukey_hsd(CF~Session)

b<-
  t_test(test2[test2$Group=="BED REST",], CF~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))


all_pvals["BDC-HDT55","CF"]<-b$p
all_pvals["CON-ME","CF"]<-a$p.adj[2]
all_pvals["CON-LC","CF"]<-a$p.adj[1]
all_pvals["LC-ME","CF"]<-a$p.adj[3]
all_pvals["BED_REST_test","CF"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","CF"]<-"anova_tukey_post_hoc"



# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55",paste(format(round(b$p,3), drop0trailing=F)),"BED REST",
  "CON","ME","p<0.001","POST-VIRAL",
  "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

CF_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=CF, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("Capillary:Fibre Ratio"))+
  scale_y_continuous(limits=c(0,5), breaks=seq(0,4,1),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(3), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        # axis.ticks = element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90,hjust=0.25, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=CF_a,
       filename = paste(output_folder,"/CF_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

CF_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=CF, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("Capillary:Fibre Ratio"))+
  scale_y_continuous(limits=c(0,5), breaks=seq(0,4,1),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(4.6), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(4.2), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[4,], label = "p.adj",y.position=c(3.2), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90,hjust=0.25, margin = margin(t = 0, r = 0, b = 1.5, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=CF_b,
       filename = paste(output_folder,"/CF_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")



# capillary density -------------------------------------------------------

# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","CD"],na.rm=TRUE); sd(data[data$Session=="BDC","CD"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","CD"],na.rm=TRUE);sd(data[data$Session=="HDT55","CD"],na.rm=TRUE)
# mean(data[data$Session=="CON","CD"],na.rm=TRUE);sd(data[data$Session=="CON","CD"],na.rm=TRUE)
# mean(data[data$Session=="LC","CD"],na.rm=TRUE);sd(data[data$Session=="LC","CD"],na.rm=TRUE)
# mean(data[data$Session=="ME","CD"],na.rm=TRUE);sd(data[data$Session=="ME","CD"],na.rm=TRUE)


remove<-data[is.na(data$CD)==TRUE ,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","CD"]);qqline(test2[test2$Session=="BDC","CD"])
# qqnorm(test2[test2$Session=="HDT55","CD"]);qqline(test2[test2$Session=="HDT55","CD"])
# qqnorm(test2[test2$Session=="CON","CD"]);qqline(test2[test2$Session=="CON","CD"])
# qqnorm(test2[test2$Session=="LC","CD"]);qqline(test2[test2$Session=="LC","CD"])
# qqnorm(test2[test2$Session=="ME","CD"]);qqline(test2[test2$Session=="ME","CD"])

model <-lme(CD~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
all_pvals["ANOVA","CD"]<-a$`Pr(>Chisq)`

# post-hoc testing
a<-
  test2[test2$Group=="POST-VIRAL",] %>% tukey_hsd(CD~Session)


b<-
  t_test(test2[test2$Group=="BED REST",], CD~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))



all_pvals["BDC-HDT55","CD"]<-b$p
all_pvals["CON-ME","CD"]<-a$p.adj[2]
all_pvals["CON-LC","CD"]<-a$p.adj[1]
all_pvals["LC-ME","CD"]<-a$p.adj[3]
all_pvals["BED_REST_test","CD"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","CD"]<-"anova_tukey_post_hoc"



# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55",paste(format(round(b$p,3), drop0trailing=F)),"BED REST",
  "CON","ME",paste("p<0.001"),"POST-VIRAL",
  "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

CD_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=CD, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(atop("Capillary Density" ,"(Capillaries" ~mm^-2*")")))+
  scale_y_continuous(limits=c(0,850), breaks=seq(0,750,250),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(550), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=CD_a,
       filename = paste(output_folder,"/CD_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

CD_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=CD, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(atop("Capillary Density" ,"(Capillaries" ~mm^-2*")")))+
  scale_y_continuous(limits=c(0,850), breaks=seq(0,750,250),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(725), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(675), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[4,], label = "p.adj",y.position=c(625), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=CD_b,
       filename = paste(output_folder,"/CD_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")



# capillary fibre ratio vs FCSA -------------------------------------------

remove<-data[is.na(data$FCSA)==TRUE | is.na(data$CF)==TRUE  ,"Subject"]
plot_data<-data[!data$Subject %in% remove, ]

# df<-list(plot_data[plot_data$Session=="BDC",],plot_data[plot_data$Session=="HDT55",])
# cocor(~FCSA+CF| FCSA  + CF, df )
# 
# df<-list(plot_data[plot_data$Session=="CON",],plot_data[plot_data$Session=="LC",])
# cocor(~FCSA+CF| FCSA  + CF, df )
# 
# df<-list(plot_data[plot_data$Session=="CON",],plot_data[plot_data$Session=="ME",])
# cocor(~FCSA+CF| FCSA  + CF, df )
# 
# df<-list(plot_data[plot_data$Session=="LC",],plot_data[plot_data$Session=="ME",])
# cocor(~FCSA+CF| FCSA  + CF, df )


y<-plot_data$CF
x<-plot_data$FCSA

Session <- plot_data$Session

model<-lm(CF~FCSA*Session, data=plot_data[plot_data$Group=="BED REST",])
model<-lm(CF~FCSA*Session, data=plot_data[plot_data$Group=="POST-VIRAL",])

BDC_cor<-cor.test(x[Session == "BDC"], y[Session == "BDC"], method="pearson")
HDT55_cor<-cor.test(x[Session == "HDT55"], y[Session == "HDT55"], method="pearson")
CON_cor<- cor.test(x[Session == "CON"], y[Session == "CON"], method="pearson")
LC_cor<- cor.test(x[Session == "LC"], y[Session == "LC"], method="pearson")
ME_cor<- cor.test(x[Session == "ME"], y[Session == "ME"], method="pearson")
cor_data_1<-rbind(BDC_cor$p.value,HDT55_cor$p.value,CON_cor$p.value,LC_cor$p.value,ME_cor$p.value)
cor_data_2<-rbind(BDC_cor$estimate,HDT55_cor$estimate,CON_cor$estimate,LC_cor$estimate,ME_cor$estimate)
cor_data_3<-c("BDC","HDT55","CON","LC","ME")
cor_data_4<-c("BED REST","BED REST","POST-VIRAL","POST-VIRAL","POST-VIRAL")
cor_data<-as.data.frame(cbind(cor_data_1,cor_data_2, cor_data_3, cor_data_4))
colnames(cor_data)<-c("p_value","r","Session","Group")
cor_data$p_value<-as.numeric(cor_data$p_value)
cor_data$p_value<-signif(cor_data$p_value,2)
cor_data$r<-as.numeric(cor_data$r)
cor_data$r<-round(cor_data$r,2)
cor_data<-add_significance(cor_data, p.col = "p_value" ,cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
cor_data$comparison<-"CF_v_FCSA"

summary(lm(CF~FCSA, data=plot_data[plot_data$Session=="CON",]))
summary(lm(CF~FCSA, data=plot_data[plot_data$Session=="LC",]))
summary(lm(CF~FCSA, data=plot_data[plot_data$Session=="ME",]))


cor_pvals<-as.data.frame(matrix(ncol = 6))
colnames(cor_pvals)<-colnames(cor_data)
cor_pvals<-cor_data


# Make graph --------------------------------------------------------------
lab_set_a<-data.frame(x1=c(3000,3000),y1=c(3.6,3.4),Group=c("BED REST","BED REST"),
                      labs=c(paste("PRE: r=",format(round(cor_data[cor_data$Session=="BDC","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="BDC","p_value"],3),nsmall=3)),
                             paste("POST: r=",format(round(cor_data[cor_data$Session=="HDT55","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="HDT55","p_value"],3),nsmall=3))))



lab_set_b<-data.frame(x1=c(3000,3000,3000),y1=c(3.8,3.6,3.4),Group=c("POST-VIRAL","POST-VIRAL","POST-VIRAL"),
                      labs=c(paste("CON: r=",format(round(cor_data[cor_data$Session=="CON","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="CON","p_value"],3),nsmall=3)),
                             paste("LC: r=",format(round(cor_data[cor_data$Session=="LC","r"],3),nsmall=2),", p<0.001"),
                             paste("ME: r=",format(round(cor_data[cor_data$Session=="ME","r"],3),nsmall=2),", p<0.001")))
CF_v_FCSA_a<-
  ggplot(data=plot_data[plot_data$Group=="BED REST",], aes(y=CF, x=FCSA))+
  geom_jitter(size=7.5, shape=21,aes(fill=fct_relevel(Session, "BDC","HDT55","CON","LC","ME")),width = 0, height=0, stroke=2)+  
  scale_fill_manual(values = c("CON"="white", "LC"=colours3[3], "ME"=colours3[4],"BDC"="white","HDT55"=colours3[2]),
                    labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON","LC"="LC","ME"="ME"))+
  scale_colour_manual(values = c("CON"="white", "LC"=colours3[3], "ME"=colours3[4],"BDC"="white","HDT55"=colours3[2]),
                      labels=c("BDC"="Bedrest-PRE", "HDT55"="POST","CON"="CON","LC"="LC","ME"="ME"))+
  geom_smooth(data=plot_data[plot_data$Session=="BDC",],method="lm",formula = y~x, se=FALSE, linetype= "solid",fill="white",colour="black", size=2, fullrange=F)+
  geom_smooth(data=plot_data[plot_data$Session=="HDT55",],method="lm",formula = y~x, se=FALSE, linetype= "dashed",fill="white",colour="black", size=2, fullrange=F)+
  xlab(expression("FCSA("*mu*m^2*")"))+
  ylab(expression("Capillary : Fibre Ratio"))+
  scale_x_continuous(limits=c(0,7500), breaks=c(2500,5000,7500), expand=c(0,0))+
  scale_y_continuous(limits = c(0,4.5), breaks=seq(0,4,1), expand=c(0,0))+
  theme(
    aspect.ratio = 1/1.2,
    legend.position = "none",
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.line=element_line(colour="black", size=line_size/2),
    axis.ticks = element_line(colour="black", size = line_size/2),
    text=element_text(size=base_size, family = "Arial"),
    panel.background  = element_rect(fill="white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      size = rel((title_text_rel_size + base_size) / base_size),
      hjust = 0.5
    ),
    axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
    axis.title.y = element_text(angle = 90, vjust = 1,size = rel((title_text_rel_size + base_size) / base_size),),
    axis.text.y = element_text( size = rel((title_text_rel_size + base_size) / base_size)
    ),
    axis.text.x=element_text(angle=0, size=rel((title_text_rel_size + base_size) / base_size)),
    strip.background = element_blank(),
    strip.text = element_blank())+
  guides(y=guide_axis(cap="upper"))+
  geom_text(data=lab_set_a,aes(y=y1,x=x1, group=Group,label=labs), family="Arial", fontface="bold", size=base_size/4, inherit.aes = TRUE )


ggsave(plot=CF_v_FCSA_a,
       filename = paste(output_folder,"/CF_v_FCSA_a_",date,".png", sep = ""),
       device="png",  width = 20, height = 16, units = "in")



CF_v_FCSA_b<-
  ggplot(data=plot_data[plot_data$Group=="POST-VIRAL",], aes(y=CF, x=FCSA))+
  geom_jitter(size=7.5, shape=21,aes(fill=fct_relevel(Session, "BDC","HDT55","CON","LC","ME")),width = 0, height=0, stroke=2)+
  scale_fill_manual(values = c("CON"="white", "LC"=colours3[3], "ME"=colours3[4],"BDC"="white","HDT55"=colours3[2]),
                    labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON","LC"="LC","ME"="ME"))+
  scale_colour_manual(values = c("CON"="white", "LC"=colours3[3], "ME"=colours3[4],"BDC"="white","HDT55"=colours3[2]),
                      labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON","LC"="LC","ME"="ME"))+
  geom_smooth(data=plot_data[plot_data$Session=="CON",],method="lm",formula = y~x, se=FALSE, linetype= "solid",fill="white",colour="black", size=2, fullrange=F)+
  geom_smooth(data=plot_data[plot_data$Session=="LC",],method="lm",formula = y~x, se=FALSE, linetype= "solid",fill="white",colour=colours3[3], size=2, fullrange=F)+
  geom_smooth(data=plot_data[plot_data$Session=="ME",],method="lm",formula = y~x, se=FALSE, linetype= "solid",fill="white",colour=colours3[4], size=2, fullrange=F)+
  xlab(expression("FCSA("*mu*m^2*")"))+
  ylab(expression("Capillary : Fibre Ratio"))+
  scale_x_continuous(limits=c(0,7500), breaks=c(2500,5000,7500), expand=c(0,0))+
  scale_y_continuous(limits = c(0,4.5), breaks=seq(0,4,1), expand=c(0,0))+
  theme(
    aspect.ratio = 1/1.2,
    legend.position = "none",
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.line=element_line(colour="black", size=line_size/2),
    axis.ticks = element_line(colour="black", size = line_size/2),
    text=element_text(size=base_size, family = "Arial"),
    panel.background  = element_rect(fill="white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      size = rel((title_text_rel_size + base_size) / base_size),
      hjust = 0.5
    ),
    axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
    axis.title.y = element_text(angle = 90, vjust = 1,size = rel((title_text_rel_size + base_size) / base_size),),
    axis.text.y = element_text( size = rel((title_text_rel_size + base_size) / base_size)
    ),
    axis.text.x=element_text(angle=0, size=rel((title_text_rel_size + base_size) / base_size)),
    strip.background = element_blank(),
    strip.text = element_blank())+
  guides(y=guide_axis(cap="upper"))+
  geom_text(data=lab_set_b,aes(y=y1,x=x1, group=Group,label=labs), family="Arial", fontface="bold", size=base_size/4, inherit.aes = TRUE )


ggsave(plot=CF_v_FCSA_b,
       filename = paste(output_folder,"/CF_v_FCSA_b_",date,".png", sep = ""),
       device="png",  width = 20, height = 16, units = "in")


# Supplemental 1 ----------------------------------------------------------


selection<-c("K","N","W","G","U","B","T","R1","J","C","L","E","F","D","H","V","M",
             "P110151","P110150","P110106","P110089","P110088","P110096","P110139","P110061","P110140","P110109","P110090","P110097","P110110",
             "P110099","P110132","P110152","P110098")

BC_data<-data[data$Subject %in% selection & data$Session!="HDT55",]

BC_pvals<-as.data.frame(matrix(ncol=90, nrow=1))

colnames(BC_pvals)<-colnames(BC_data)
rownames(BC_pvals)<-c("CON-BDC")

median(BC_data[BC_data$Session=="BDC","Age"],na.rm=TRUE);quantile(BC_data[BC_data$Session=="BDC","Age"], na.rm=TRUE)
median(BC_data[BC_data$Session=="CON","Age"],na.rm=TRUE);quantile(BC_data[BC_data$Session=="CON","Age"], na.rm=TRUE)

median(BC_data[BC_data$Session=="BDC","Height"],na.rm=TRUE);quantile(BC_data[BC_data$Session=="BDC","Height"], na.rm=TRUE)
median(BC_data[BC_data$Session=="CON","Height"],na.rm=TRUE);quantile(BC_data[BC_data$Session=="CON","Height"], na.rm=TRUE)

median(BC_data[BC_data$Session=="BDC","Weight"],na.rm=TRUE);quantile(BC_data[BC_data$Session=="BDC","Weight"], na.rm=TRUE)
median(BC_data[BC_data$Session=="CON","Weight"],na.rm=TRUE);quantile(BC_data[BC_data$Session=="CON","Weight"], na.rm=TRUE)

median(BC_data[BC_data$Session=="CON","Steps"],na.rm=TRUE);quantile(BC_data[BC_data$Session=="CON","Steps"], na.rm=TRUE)

length(BC_data[BC_data$Session=="BDC" & BC_data$Sex=="Female","Age"])/length(BC_data[BC_data$Session=="BDC","Age"])
length(BC_data[BC_data$Session=="CON" & BC_data$Sex=="Female","Age"])/length(BC_data[BC_data$Session=="CON","Age"])

length(BC_data$Age)

BC_test<-as.data.frame(box_cox_transform(BC_data))
BC_test[,c("Subject","Session","Group", "Sex")]<-BC_data[,c("Subject","Session","Group","Sex")]

remove<-BC_test[is.na(BC_test$Age)==TRUE,"Subject"]
test2<-BC_test[!BC_test$Subject %in% remove, ]

qqnorm(test2[test2$Session=="CON","Age"]);qqline(test2[test2$Session=="CON","Age"]) 
qqnorm(test2[test2$Session=="CON","Height"]);qqline(test2[test2$Session=="CON","Height"]) 
qqnorm(test2[test2$Session=="CON","Weight"]);qqline(test2[test2$Session=="CON","Weight"]) 


t_test(test2,Age~Session, paired=F)
kruskal.test(Weight~Session, test2)
kruskal.test(Height~Session, test2)


# BC VO2 ------------------------------------------------------------------
remove<-BC_data[is.na(BC_data$VO2_abs)==TRUE   ,"Subject"]
BC_test2<-BC_test[!BC_test$Subject %in% remove,]
# shapiro.test(BC_test2[BC_test2$Session=="BDC","VO2_abs"])
# shapiro.test(BC_test2[BC_test2$Session=="CON","VO2_abs"]) 
# qqnorm(BC_test2[BC_test2$Session=="BDC","VO2_abs"]);qqline(BC_test2[BC_test2$Session=="BDC","VO2_abs"]) 
# qqnorm(BC_test2[BC_test2$Session=="CON","VO2_abs"]);qqline(BC_test2[BC_test2$Session=="CON","VO2_abs"]) 

a<-
  wilcox.test(BC_test2[BC_test2$Session=="BDC","VO2_abs"],test2[test2$Session=="CON","VO2_abs"],paired=F)

BC_pvals["CON-BDC","VO2_abs"]<-a$p.value


stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","CON",paste(format(round(a$p.value[1],3), drop0trailing=F)),"BED REST")


plot_data<-BC_data[!BC_data$Subject %in% remove,]

VO2_BC<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"BDC","CON"),y=VO2_abs, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("BDC"=colours3[2], "CON"="white"))+
  ylab(expression("V\U0307" ~O[2]*"max (L"~min^-1*")"))+
  scale_y_continuous(limits=c(0,4.9), breaks=seq(0,4,1),expand=c(0,0))+
  scale_x_discrete(labels=c( "BDC"=expression(atop("PRE","BED REST")),"CON"="CON"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(4.3), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        # axis.ticks = element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size-10)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=VO2_BC,
       filename = paste(output_folder,"/VO2_BC_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")



# BC GET ------------------------------------------------------------------
remove<-BC_data[is.na(BC_data$GET)==TRUE   ,"Subject"]
BC_test2<-BC_test[!BC_test$Subject %in% remove,]
# shapiro.test(BC_test2[BC_test2$Session=="BDC","GET"])
# shapiro.test(BC_test2[BC_test2$Session=="CON","GET"])
# qqnorm(BC_test2[BC_test2$Session=="BDC","GET"]);qqline(BC_test2[BC_test2$Session=="BDC","GET"]) 
# qqnorm(BC_test2[BC_test2$Session=="CON","GET"]);qqline(BC_test2[BC_test2$Session=="CON","GET"]) 

model<-lme(GET~Session, data=BC_test2, random = ~ 1|Subject, na.action = na.omit, control="optim")
car::Anova(model)
a<-
  t.test(BC_test2[BC_test2$Session=="CON","GET"], BC_test2[BC_test2$Session=="BDC","GET"], paired=FALSE)

BC_pvals["CON-BDC","GET"]<-a$p.value


stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","CON",paste(format(round(a$p.value[1],3), drop0trailing=F)),"BED REST")


plot_data<-BC_data[!BC_data$Subject %in% remove,]

GET_BC<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"BDC","CON"),y=GET, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("BDC"=colours3[2], "CON"="white"))+
  ylab(expression("GET (L"~min^-1*")"))+
  scale_y_continuous(limits=c(0,3.4), breaks=seq(0,3,1),expand=c(0,0))+
  scale_x_discrete(labels=c( "BDC"=expression(atop("PRE","BED REST")),"CON"="CON"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(3.15), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size-10)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=GET_BC,
       filename = paste(output_folder,"/GET_BC_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")



# BC GET % ----------------------------------------------------------------
remove<-BC_data[is.na(BC_data$GET_perc)==TRUE   ,"Subject"]
BC_test2<-BC_test[!BC_test$Subject %in% remove,]
# shapiro.test(BC_test2[BC_test2$Session=="BDC","GET_perc"])
# shapiro.test(BC_test2[BC_test2$Session=="CON","GET_perc"])
# qqnorm(BC_test2[BC_test2$Session=="BDC","GET_perc"]);qqline(BC_test2[BC_test2$Session=="BDC","GET_perc"]) 
# qqnorm(BC_test2[BC_test2$Session=="CON","GET_perc"]);qqline(BC_test2[BC_test2$Session=="CON","GET_perc"]) 

model<-lme(GET_perc~Session, data=BC_test2, random = ~ 1|Subject, na.action = na.omit, control="optim")
car::Anova(model)
a<-
  t.test(BC_test2[BC_test2$Session=="CON","GET_perc"], BC_test2[BC_test2$Session=="BDC","GET_perc"], paired=FALSE)

BC_pvals["CON-BDC","GET_perc"]<-a$p.value


stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","CON",paste(format(round(a$p.value[1],3), nsmall=3,drop0trailing=F)),"BED REST")


plot_data<-BC_data[!BC_data$Subject %in% remove,]

GET_perc_BC<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"BDC","CON"),y=GET_perc*100, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("BDC"=colours3[2], "CON"="white"))+
  ylab(expression("GET (% of V\U0307"~O[2]*"max)"))+
  scale_y_continuous(limits=c(0,110), breaks=c(0,40,60,80,100),expand=c(0,0))+
  scale_x_discrete(labels=c( "BDC"=expression(atop("PRE","BED REST")),"CON"="CON"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(90), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(30),
              space=c(0.5),
              scales=c(20),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size-10)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=GET_perc_BC,
       filename = paste(output_folder,"/GET_perc_BC_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

# BC SDH ------------------------------------------------------------------
remove<-BC_data[is.na(BC_data$SDH)==TRUE   ,"Subject"]
BC_test2<-BC_test[!BC_test$Subject %in% remove,]
# shapiro.test(BC_test2[BC_test2$Session=="BDC","SDH"]) 
# shapiro.test(BC_test2[BC_test2$Session=="CON","SDH"])
# qqnorm(BC_test2[BC_test2$Session=="BDC","SDH"]);qqline(BC_test2[BC_test2$Session=="BDC","SDH"]) 
# qqnorm(BC_test2[BC_test2$Session=="CON","SDH"]);qqline(BC_test2[BC_test2$Session=="CON","SDH"]) 



a<-
  wilcox.test(BC_test2[BC_test2$Session=="BDC","SDH"],test2[test2$Session=="CON","SDH"],paired=F)


BC_pvals["CON-BDC","SDH"]<-a$p.value


stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","CON",paste(format(round(a$p.value[1],3), drop0trailing=F)),"BED REST")


plot_data<-BC_data[!BC_data$Subject %in% remove,]

SDH_BC<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"BDC","CON"),y=SDH*10^5, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("BDC"=colours3[2], "CON"="white"))+
  ylab(expression(atop("SDH Activity", "(\U0394"~A[660]*~mu*m^-1*""*~s^-1*~10^-5*")")))+
  scale_y_continuous(limits=c(0,1.8), breaks=seq(0,1.5,0.5),expand=c(0,0))+
  scale_x_discrete(labels=c( "BDC"=expression(atop("PRE","BED REST")),"CON"="CON"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(1.65), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        # axis.ticks = element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size-10)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=SDH_BC,
       filename = paste(output_folder,"/SDH_BC_BC_",date,".png", sep = ""),
       device="png",  width = 10, height = 14, units = "in")



# BC Oxphos ---------------------------------------------------------------
remove<-BC_data[is.na(BC_data$Oxphos)==TRUE | BC_data$membrane_intact>1.1  ,"Subject"]
BC_test2<-BC_test[!BC_test$Subject %in% remove,]
# shapiro.test(BC_test2[BC_test2$Session=="BDC","Oxphos"])
# shapiro.test(BC_test2[BC_test2$Session=="CON","Oxphos"])
# qqnorm(BC_test2[BC_test2$Session=="BDC","Oxphos"]);qqline(BC_test2[BC_test2$Session=="BDC","Oxphos"]) 
# qqnorm(BC_test2[BC_test2$Session=="CON","Oxphos"]);qqline(BC_test2[BC_test2$Session=="CON","Oxphos"]) 

model<-lme(Oxphos~Session, data=BC_test2, random = ~ 1|Subject, na.action = na.omit, control="optim")
car::Anova(model)
a<-
  t.test(BC_test2[BC_test2$Session=="CON","Oxphos"], BC_test2[BC_test2$Session=="BDC","Oxphos"], paired=FALSE)

BC_pvals["CON-BDC","Oxphos"]<-a$p.value

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","CON",paste(format(round(a$p.value[1],3), drop0trailing=F)),"BED REST")


plot_data<-BC_data[!BC_data$Subject %in% remove,]

Oxphos_BC<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"BDC","CON"),y=Oxphos, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("BDC"=colours3[2], "CON"="white"))+
  ylab(bquote(atop("Oxidative Phosphorylation","(pmol"*~s^-1*~mg^-1*")")))+
  scale_y_continuous(limits=c(0,180), breaks=seq(0,150,50),expand=c(0,0))+
  scale_x_discrete(labels=c( "BDC"=expression(atop("PRE","BED REST")),"CON"="CON"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(150), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size-10)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=Oxphos_BC,
       filename = paste(output_folder,"/Oxphos_BC_",date,".png", sep = ""),
       device="png",  width = 10, height = 14, units = "in")


# BC Cap: fibre -----------------------------------------------------------
remove<-BC_data[is.na(BC_data$CF)==TRUE   ,"Subject"]
BC_test2<-BC_test[!BC_test$Subject %in% remove,]
# shapiro.test(BC_test2[BC_test2$Session=="BDC","CF"])
# shapiro.test(BC_test2[BC_test2$Session=="CON","CF"])
# qqnorm(BC_test2[BC_test2$Session=="BDC","CF"]);qqline(BC_test2[BC_test2$Session=="BDC","CF"]) 
# qqnorm(BC_test2[BC_test2$Session=="CON","CF"]);qqline(BC_test2[BC_test2$Session=="CON","CF"]) 

model<-lme(CF~Session, data=BC_test2, random = ~ 1|Subject, na.action = na.omit, control="optim")
car::Anova(model)
a<-
  t.test(BC_test2[BC_test2$Session=="CON","CF"], BC_test2[BC_test2$Session=="BDC","CF"], paired=FALSE)

BC_pvals["CON-BDC","CF"]<-a$p.value

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","CON",paste(format(round(a$p.value[1],3), drop0trailing=F)),"BED REST")


plot_data<-BC_data[!BC_data$Subject %in% remove,]

CF_BC<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"BDC","CON"),y=CF, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("BDC"=colours3[2], "CON"="white"))+
  ylab(expression("Capillary:Fibre Ratio"))+
  scale_y_continuous(limits=c(0,5), breaks=seq(0,4,1),expand=c(0,0))+
  scale_x_discrete(labels=c( "BDC"=expression(atop("PRE","BED REST")),"CON"="CON"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(3.2), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size-10)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=CF_BC,
       filename = paste(output_folder,"/CF_BC_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")




# BC FCSA -----------------------------------------------------------------
remove<-BC_data[is.na(BC_data$FCSA)==TRUE   ,"Subject"]
BC_test2<-BC_test[!BC_test$Subject %in% remove,]
# shapiro.test(BC_test2[BC_test2$Session=="BDC","FCSA"])
# shapiro.test(BC_test2[BC_test2$Session=="CON","FCSA"])
# qqnorm(BC_test2[BC_test2$Session=="BDC","FCSA"]);qqline(BC_test2[BC_test2$Session=="BDC","FCSA"]) 
# qqnorm(BC_test2[BC_test2$Session=="CON","FCSA"]);qqline(BC_test2[BC_test2$Session=="CON","FCSA"]) 

a<-
  wilcox.test(BC_test2[BC_test2$Session=="BDC","FCSA"],test2[test2$Session=="CON","FCSA"],paired=F)


BC_pvals["CON-BDC","FCSA"]<-a$p.value

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","CON",paste(format(round(a$p.value[1],3), drop0trailing=F)),"BED REST")


plot_data<-BC_data[!BC_data$Subject %in% remove,]

FCSA_BC<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"BDC","CON"),y=FCSA, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("BDC"=colours3[2], "CON"="white"))+
  ylab(expression("FCSA ("*mu*m^2*")"))+
  scale_y_continuous(limits=c(0,9000), breaks=seq(0,8000,2000),expand=c(0,0))+
  scale_x_discrete(labels=c( "BDC"=expression(atop("PRE","BED REST")),"CON"="CON"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(8100), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text(angle=0, hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size-10)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=FCSA_BC,
       filename = paste(output_folder,"/FCSA_BC_",date,".png", sep = ""),
       device="png",  width = 9.2, height = 14, units = "in")



# BC VO2 vs SDH -----------------------------------------------------------
remove<-BC_data[is.na(BC_data$SDH)==TRUE | is.na(BC_data$VO2_rel)==TRUE  ,"Subject"]
plot_data<-BC_data[!BC_data$Subject %in% remove, ]

df<-list(plot_data[plot_data$Session=="BDC",],plot_data[plot_data$Session=="CON",])
cocor(~VO2_rel+SDH| VO2_rel+SDH, df )

x<-BC_data$SDH
y<-BC_data$VO2_rel
Session<-BC_data$Session
BDC_cor<-cor.test(x[Session == "BDC"], y[Session == "BDC"], method="pearson")
CON_cor<- cor.test(x[Session == "CON"], y[Session == "CON"], method="pearson")
cor_data_1<-rbind(BDC_cor$p.value,CON_cor$p.value)
cor_data_2<-rbind(BDC_cor$estimate,CON_cor$estimate)
cor_data_3<-c("BDC","CON")
cor_data_4<-c("BED REST","POST-VIRAL")
cor_data<-as.data.frame(cbind(cor_data_1,cor_data_2, cor_data_3, cor_data_4))
colnames(cor_data)<-c("p_value","r","Session","Group")
cor_data$p_value<-as.numeric(cor_data$p_value)
cor_data$p_value<-signif(cor_data$p_value,2)
cor_data$r<-as.numeric(cor_data$r)
cor_data$r<-round(cor_data$r,2)
cor_data<-add_significance(cor_data, p.col = "p_value" ,cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
cor_data$comparison<-"SDH_v_VO2_BC"

cor_pvals<-rbind(cor_pvals,cor_data)

lab_set<-data.frame(x1=c(0.9,0.9),y1=c(65,60),Group=c("BED REST","POST-VIRAL"),
                    labs=c(paste("PRE-Bedrest: r=",format(round(cor_data[cor_data$Session=="BDC","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="BDC","p_value"],3),nsmall=3)),
                           paste("CON: r=",format(round(cor_data[cor_data$Session=="CON","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="CON","p_value"],3),nsmall=3))))


VO2_v_SDH_BC<-
  ggplot(data=plot_data, aes(y=VO2_rel, x=SDH*10^5))+
  geom_jitter(size=7.5,aes( shape=fct_relevel(Session, "BDC","CON"),fill=fct_relevel(Session, "BDC","CON")),width = 0, height=0, stroke=2)+
  scale_shape_manual(values=c("BDC"=21,"CON"=24),
                     labels=c("BDC"="PRE-Bed Rest","CON"="CON"))+
  scale_fill_manual(values = c("BDC"=colours3[2], "CON"="white"),
                    labels=c("BDC"="PRE-Bed Rest","CON"="CON"))+
  scale_colour_manual(values = c("BDC"="white","CON"="white"),
                      labels=c("BDC"="PRE-Bed Rest","CON"="CON"))+
  # geom_smooth(data=plot_data[plot_data$Session=="BDC",],method="lm", se=FALSE, linetype= "solid",fill="white",colour="black", size=2)+
  geom_smooth(data=plot_data[plot_data$Session=="CON",],method="lm",formula = y~x, se=FALSE, linetype= "solid",fill="white",colour="black", size=2, fullrange=F)+
  xlab(expression(atop("SDH Activity", "(\U0394"~A[660]*~mu*m^-1*""*~s^-1*~10^-5*")")))+
  ylab(expression("V\U0307" ~O[2][max]*"(mL"~min^-1*~kg^-1*")"))+
  scale_x_continuous(limits=c(0,2.2), breaks=seq(0.5,2.0,0.5), expand=c(0,0))+
  scale_y_continuous(limits = c(0,80), breaks=seq(0,60,20), expand=c(0,0))+
  theme(
    aspect.ratio = 1/1.2,
    legend.position = "none",
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.line=element_line(colour="black", size=line_size/2),
    axis.ticks = element_line(colour="black", size = line_size/2),
    text=element_text(size=base_size, family = "Arial"),
    panel.background  = element_rect(fill="white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      size = rel((title_text_rel_size + base_size) / base_size),
      hjust = 0.5
    ),
    axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
    axis.title.y = element_text(angle = 90, vjust = 1,size = rel((title_text_rel_size + base_size) / base_size),),
    axis.text.y = element_text( size = rel((title_text_rel_size + base_size) / base_size)
    ),
    axis.text.x=element_text(angle=0, size=rel((title_text_rel_size + base_size) / base_size)),
    strip.background = element_blank(),
    strip.text = element_blank())+
  guides(y=guide_axis(cap="upper"))+
  geom_text(data=lab_set,aes(y=y1,x=x1, group=Group,label=labs), family="Arial", fontface="bold", size=base_size/4, inherit.aes = TRUE )


ggsave(plot=VO2_v_SDH_BC,
       filename = paste(output_folder,"/VO2_v_SDH_BC_",date,".png", sep = ""),
       device="png",  width = 20, height = 16, units = "in")



# BC VO2 vs Oxphos --------------------------------------------------------
remove<-BC_data[is.na(BC_data$Oxphos)==TRUE | is.na(BC_data$VO2_rel)==TRUE & BC_data$membrane_intact>1.1  ,"Subject"]
plot_data<-BC_data[!BC_data$Subject %in% remove, ]

df<-list(plot_data[plot_data$Session=="BDC",],plot_data[plot_data$Session=="CON",])
cocor(~VO2_rel+Oxphos| VO2_rel+Oxphos, df )

x<-BC_data$Oxphos
y<-BC_data$VO2_rel
Session<-BC_data$Session
BDC_cor<-cor.test(x[Session == "BDC"], y[Session == "BDC"], method="pearson")
CON_cor<- cor.test(x[Session == "CON"], y[Session == "CON"], method="pearson")
cor_data_1<-rbind(BDC_cor$p.value,CON_cor$p.value)
cor_data_2<-rbind(BDC_cor$estimate,CON_cor$estimate)
cor_data_3<-c("BDC","CON")
cor_data_4<-c("BED REST","POST-VIRAL")
cor_data<-as.data.frame(cbind(cor_data_1,cor_data_2, cor_data_3, cor_data_4))
colnames(cor_data)<-c("p_value","r","Session","Group")
cor_data$p_value<-as.numeric(cor_data$p_value)
cor_data$p_value<-signif(cor_data$p_value,2)
cor_data$r<-as.numeric(cor_data$r)
cor_data$r<-round(cor_data$r,2)
cor_data<-add_significance(cor_data, p.col = "p_value" ,cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
cor_data$comparison<-"Oxphos_v_VO2_BC"

cor_pvals<-rbind(cor_pvals,cor_data)

lab_set<-data.frame(x1=c(80,80),y1=c(70,65),Group=c("BED REST","POST-VIRAL"),
                    labs=c(paste("PRE-Bedrest: r=",format(round(cor_data[cor_data$Session=="BDC","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="BDC","p_value"],3),nsmall=3)),
                           paste("CON: r=",format(round(cor_data[cor_data$Session=="CON","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="CON","p_value"],3),nsmall=3))))


VO2_v_Oxphos_BC<-
  ggplot(data=plot_data, aes(y=VO2_rel, x=Oxphos))+
  geom_jitter(size=7.5,aes( shape=fct_relevel(Session, "BDC","CON"),fill=fct_relevel(Session, "BDC","CON")),width = 0, height=0, stroke=2)+
  scale_shape_manual(values=c("BDC"=21,"CON"=24),
                     labels=c("BDC"="PRE","CON"="CON"))+
  scale_fill_manual(values = c("BDC"=colours3[2], "CON"="white"),
                    labels=c("BDC"="PRE-Bed Rest","CON"="CON"))+
  scale_colour_manual(values = c("BDC"="white","CON"="white"),
                      labels=c("BDC"="PRE-Bed Rest","CON"="CON"))+
  geom_smooth(data=plot_data[plot_data$Session=="BDC",],method="lm",formula = y~x, se=FALSE, linetype= "solid",fill="white",colour="black", size=2, fullrange=F)+
  geom_smooth(data=plot_data[plot_data$Session=="CON",],method="lm",formula = y~x, se=FALSE, linetype= "dashed",fill="white",colour="black", size=2, fullrange=F)+
  xlab(bquote(atop("Oxidative Phosphorylation","(pmol"*~s^-1*~mg^-1*")")))+
  ylab(expression("V\U0307" ~O[2][max]*"(mL"~min^-1*~kg^-1*")"))+
  scale_y_continuous(limits = c(0,75), breaks=seq(0,60,20), expand=c(0,0))+
  scale_x_continuous(limits=c(0,160), breaks=seq(30,150,30), expand=c(0,0))+
  theme(
    aspect.ratio = 1/1.2,
    legend.position = "none",
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.line=element_line(colour="black", size=line_size/2),
    axis.ticks = element_line(colour="black", size = line_size/2),
    text=element_text(size=base_size, family = "Arial"),
    panel.background  = element_rect(fill="white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      size = rel((title_text_rel_size + base_size) / base_size),
      hjust = 0.5
    ),
    axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
    axis.title.y = element_text(angle = 90, vjust = 1,size = rel((title_text_rel_size + base_size) / base_size),),
    axis.text.y = element_text( size = rel((title_text_rel_size + base_size) / base_size)
    ),
    axis.text.x=element_text(angle=0, size=rel((title_text_rel_size + base_size) / base_size)),
    strip.background = element_blank(),
    strip.text = element_blank())+
  guides(y=guide_axis(cap="upper"))+
  geom_text(data=lab_set,aes(y=y1,x=x1, group=Group,label=labs), family="Arial", fontface="bold", size=base_size/4, inherit.aes = TRUE )


ggsave(plot=VO2_v_Oxphos_BC,
       filename = paste(output_folder,"/VO2_v_Oxphos_BC_",date,".png", sep = ""),
       device="png",  width = 20, height = 16, units = "in")




# BC CF vs FCSA -----------------------------------------------------------
remove<-BC_data[is.na(BC_data$FCSA)==TRUE | is.na(BC_data$CF)==TRUE  ,"Subject"]
plot_data<-BC_data[!BC_data$Subject %in% remove, ]

df<-list(plot_data[plot_data$Session=="BDC",],plot_data[plot_data$Session=="CON",])
cocor(~FCSA+CF| FCSA  + CF, df )


x<-plot_data$CF
y<-plot_data$FCSA
Session<-plot_data$Session
BDC_cor<-cor.test(x[Session == "BDC"], y[Session == "BDC"], method="pearson")
CON_cor<- cor.test(x[Session == "CON"], y[Session == "CON"], method="pearson")
cor_data_1<-rbind(BDC_cor$p.value,CON_cor$p.value)
cor_data_2<-rbind(BDC_cor$estimate,CON_cor$estimate)
cor_data_3<-c("BDC","CON")
cor_data_4<-c("BED REST","POST-VIRAL")
cor_data<-as.data.frame(cbind(cor_data_1,cor_data_2, cor_data_3, cor_data_4))
colnames(cor_data)<-c("p_value","r","Session","Group")
cor_data$p_value<-as.numeric(cor_data$p_value)
cor_data$p_value<-signif(cor_data$p_value,2)
cor_data$r<-as.numeric(cor_data$r)
cor_data$r<-round(cor_data$r,2)
cor_data<-add_significance(cor_data, p.col = "p_value" ,cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
cor_data$comparison<-"CF_v_FCSA_BC"

cor_pvals<-rbind(cor_pvals,cor_data)

lab_set<-data.frame(x1=c(3000,3000),y1=c(3.6,3.4),Group=c("BED REST","POST-VIRAL"),
                    labs=c(paste("PRE-Bedrest: r=",format(round(cor_data[cor_data$Session=="BDC","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="BDC","p_value"],3),nsmall=3)),
                           paste("CON: r=",format(round(cor_data[cor_data$Session=="CON","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="CON","p_value"],3),nsmall=3))))


CF_v_FCSA_BC<-
  ggplot(data=plot_data, aes(y=CF, x=FCSA))+
  geom_jitter(size=7.5,aes(shape=fct_relevel(Session, "BDC","CON"),fill=fct_relevel(Session, "BDC","CON")),width = 0, height=0, stroke=2)+
  scale_shape_manual(values=c("BDC"=21,"CON"=24),
                     labels=c("BDC"="PRE","CON"="CON"))+
  scale_fill_manual(values = c("BDC"=colours3[2], "CON"="white"),
                    labels=c("BDC"="PRE-Bed Rest","CON"="CON"))+
  scale_colour_manual(values = c("BDC"="white","CON"="white"),
                      labels=c("BDC"="PRE-Bed Rest","CON"="CON"))+
  geom_smooth(data=plot_data[plot_data$Session=="BDC",],method="lm",formula = y~x, se=FALSE, linetype= "solid",fill="white",colour="black", size=2, fullrange=F)+
  geom_smooth(data=plot_data[plot_data$Session=="CON",],method="lm",formula = y~x, se=FALSE, linetype= "dashed",fill="white",colour="black", size=2, fullrange=F)+
  xlab(expression("FCSA("*mu*m^2*")"))+
  ylab(expression("Capillary : Fibre Ratio"))+
  scale_x_continuous(limits=c(0,7500), breaks=c(2500,5000,7500), expand=c(0,0))+
  scale_y_continuous(limits = c(0,4.5), breaks=seq(0,4,1), expand=c(0,0))+
  theme(
    aspect.ratio = 1/1.2,
    legend.position = "none",
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.line=element_line(colour="black", size=line_size/2),
    axis.ticks = element_line(colour="black", size = line_size/2),
    text=element_text(size=base_size, family = "Arial"),
    panel.background  = element_rect(fill="white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      size = rel((title_text_rel_size + base_size) / base_size),
      hjust = 0.5
    ),
    axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
    axis.title.y = element_text(angle = 90, vjust = 1,size = rel((title_text_rel_size + base_size) / base_size),),
    axis.text.y = element_text( size = rel((title_text_rel_size + base_size) / base_size)
    ),
    axis.text.x=element_text(angle=0, size=rel((title_text_rel_size + base_size) / base_size)),
    strip.background = element_blank(),
    strip.text = element_blank())+
  guides(y=guide_axis(cap="upper"))+
  geom_text(data=lab_set,aes(y=y1,x=x1, group=Group,label=labs), family="Arial", fontface="bold", size=base_size/4, inherit.aes = TRUE )


ggsave(plot=CF_v_FCSA_BC,
       filename = paste(output_folder,"/CF_v_FCSA_BC_",date,".png", sep = ""),
       device="png",  width = 20, height = 16, units = "in")


# Supplemental 2 ----------------------------------------------------------

# Peak Power --------------------------------------------------------------
# Stats -------------------------------------------------------------------

mean(data[data$Session=="BDC","Peak_power"],na.rm=TRUE); sd(data[data$Session=="BDC","Peak_power"],na.rm=TRUE)
mean(data[data$Session=="HDT55","Peak_power"],na.rm=TRUE);sd(data[data$Session=="HDT55","Peak_power"],na.rm=TRUE)
mean(data[data$Session=="CON","Peak_power"],na.rm=TRUE);sd(data[data$Session=="CON","Peak_power"],na.rm=TRUE)
mean(data[data$Session=="LC","Peak_power"],na.rm=TRUE);sd(data[data$Session=="LC","Peak_power"],na.rm=TRUE)
mean(data[data$Session=="ME","Peak_power"],na.rm=TRUE);sd(data[data$Session=="ME","Peak_power"],na.rm=TRUE)


remove<-data[is.na(data$Peak_power)==TRUE,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","Peak_power"]);qqline(test2[test2$Session=="BDC","Peak_power"])
# qqnorm(test2[test2$Session=="HDT55","Peak_power"]);qqline(test2[test2$Session=="HDT55","Peak_power"])
# qqnorm(test2[test2$Session=="CON","Peak_power"]);qqline(test2[test2$Session=="CON","Peak_power"])
# qqnorm(test2[test2$Session=="LC","Peak_power"]);qqline(test2[test2$Session=="LC","Peak_power"])
# qqnorm(test2[test2$Session=="ME","Peak_power"]);qqline(test2[test2$Session=="ME","Peak_power"])

a<-kruskal.test(Peak_power~Session, data=test2[test2$Group=="POST-VIRAL",])
all_pvals["ANOVA","Peak_power"]<-a$p.value

a<-
  pairwise.wilcox.test(test2[test2$Group=="POST-VIRAL","Peak_power"], test2[test2$Group=="POST-VIRAL","Session"], paired=FALSE, p.adjust.method="hommel")

b<-
  wilcox.test(test2[test2$Session=="BDC","Peak_power"],test2[test2$Session=="HDT55","Peak_power"], paired=TRUE)

all_pvals["BDC-HDT55","Peak_power"]<-b$p.value
all_pvals["CON-ME","Peak_power"]<-a$p.value[2]
all_pvals["CON-LC","Peak_power"]<-a$p.value[1]
all_pvals["LC-ME","Peak_power"]<-a$p.value[4]
all_pvals["BED_REST_test","Peak_power"]<-"wilcoxon"
all_pvals["POST_VIRAL_test","Peak_power"]<-"kruskal_wilcoxon_post_hoc"

# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55","p<0.001","BED REST",
  "CON","ME","p<0.001","POST-VIRAL",
  "CON","LC","p<0.001","POST-VIRAL",
  "LC","ME",paste(format(round(a$p.value[4],3), drop0trailing=F)),"POST-VIRAL")

plot_data<-plot_data<-data[!data$Subject %in% remove,]

Power_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=Peak_power, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("Peak Power (W)"))+
  scale_y_continuous(limits=c(0,500), breaks=seq(0,400,100),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(375), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 40, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=Power_a,
       filename = paste(output_folder,"/Power_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

Power_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=Peak_power, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("Peak Power (W)"))+
  scale_y_continuous(limits=c(0,500), breaks=seq(0,400,100),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(450), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(400), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[4,], label = "p.adj",y.position=c(350), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 40, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=Power_b,
       filename = paste(output_folder,"/Power_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

# GET (% of VO2 max) ------------------------------------------------------
# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","GET_perc"],na.rm=TRUE); sd(data[data$Session=="BDC","GET_perc"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","GET_perc"],na.rm=TRUE);sd(data[data$Session=="HDT55","GET_perc"],na.rm=TRUE)
# mean(data[data$Session=="CON","GET_perc"],na.rm=TRUE);sd(data[data$Session=="CON","GET_perc"],na.rm=TRUE)
# mean(data[data$Session=="LC","GET_perc"],na.rm=TRUE);sd(data[data$Session=="LC","GET_perc"],na.rm=TRUE)
# mean(data[data$Session=="ME","GET_perc"],na.rm=TRUE);sd(data[data$Session=="ME","GET_perc"],na.rm=TRUE)

remove<-data[is.na(data$GET_perc)==TRUE ,"Subject"]
# remove<-data[is.na(data$GET_perc)==TRUE | data$Subject=="H","Subject"]
remove2<-c("P110117","P110118","P110127")
test2<-test[!test$Subject %in% remove & !test$Subject %in% remove2,]

# qqnorm(test2[test2$Session=="BDC","GET_perc"]);qqline(test2[test2$Session=="BDC","GET_perc"])
# qqnorm(test2[test2$Session=="HDT55","GET_perc"]);qqline(test2[test2$Session=="HDT55","GET_perc"])
# qqnorm(test2[test2$Session=="CON","GET_perc"]);qqline(test2[test2$Session=="CON","GET_perc"])
# qqnorm(test2[test2$Session=="LC","GET_perc"]);qqline(test2[test2$Session=="LC","GET_perc"])
# qqnorm(test2[test2$Session=="ME","GET_perc"]);qqline(test2[test2$Session=="ME","GET_perc"])


a<-kruskal.test(GET_perc~Session, test2[test2$Group=="POST-VIRAL",])
all_pvals["ANOVA","GET_perc"]<-a$p.value
a<-pairwise.wilcox.test(test2[test2$Group=="POST-VIRAL","GET_perc"], test2[test2$Group=="POST-VIRAL","Session"], p.adjust.method="BH")

b<-
  t_test(test2[test2$Group=="BED REST",], GET_perc~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))



all_pvals["BDC-HDT55","GET_perc"]<-b$p
all_pvals["CON-ME","GET_perc"]<-a$p.value[2]
all_pvals["CON-LC","GET_perc"]<-a$p.value[1]
all_pvals["LC-ME","GET_perc"]<-a$p.value[4]
all_pvals["BED_REST_test","GET_perc"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","GET_perc"]<-"kruskal_wilcoxon_post_hoc"



# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55",paste(format(round(b$p,3), drop0trailing=F)),"BED REST",
  "CON","ME","p<0.001","POST-VIRAL",
  "CON","LC",paste(format(round(a$p.value[1],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME","p<0.001","POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

GET_perc_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=GET_perc*100, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("GET (% of V\U0307"~O[2]*"max)"))+
  scale_y_continuous(limits=c(0,110), breaks=c(0,40,60,80,100),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(100), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(30),
              space=c(0.5),
              scales=c(20),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=GET_perc_a,
       filename = paste(output_folder,"/GET_perc_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

GET_perc_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=GET_perc*100, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("GET (% of V\U0307"~O[2]*"max)"))+
  scale_y_continuous(limits=c(0,110), breaks=c(0,40,60,80,100),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(105), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(90), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[4,], label = "p.adj",y.position=c(95), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(30),
              space=c(0.5),
              scales=c(20),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=GET_perc_b,
       filename = paste(output_folder,"/GET_perc_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")
# VE/VO2 max --------------------------------------------------------------

# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","EqO2_max"],na.rm=TRUE); sd(data[data$Session=="BDC","EqO2_max"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","EqO2_max"],na.rm=TRUE);sd(data[data$Session=="HDT55","EqO2_max"],na.rm=TRUE)
# mean(data[data$Session=="CON","EqO2_max"],na.rm=TRUE);sd(data[data$Session=="CON","EqO2_max"],na.rm=TRUE)
# mean(data[data$Session=="LC","EqO2_max"],na.rm=TRUE);sd(data[data$Session=="LC","EqO2_max"],na.rm=TRUE)
# mean(data[data$Session=="ME","EqO2_max"],na.rm=TRUE);sd(data[data$Session=="ME","EqO2_max"],na.rm=TRUE)


remove<-data[is.na(data$EqO2_max)==TRUE,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","EqO2_max"]);qqline(test2[test2$Session=="BDC","EqO2_max"])
# qqnorm(test2[test2$Session=="HDT55","EqO2_max"]);qqline(test2[test2$Session=="HDT55","EqO2_max"])
# qqnorm(test2[test2$Session=="CON","EqO2_max"]);qqline(test2[test2$Session=="CON","EqO2_max"])
# qqnorm(test2[test2$Session=="LC","EqO2_max"]);qqline(test2[test2$Session=="LC","EqO2_max"])
# qqnorm(test2[test2$Session=="ME","EqO2_max"]);qqline(test2[test2$Session=="ME","EqO2_max"])


model <-lme(EqO2_max~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
all_pvals["ANOVA","EqO2_max"]<-a$`Pr(>Chisq)`

b<-
  t_test(test2[test2$Group=="BED REST",], EqO2_max~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","EqO2_max"]<-b$p
# all_pvals["CON-ME","EqO2_max"]<-a$p.adj[2]
# all_pvals["CON-LC","EqO2_max"]<-a$p.adj[1]
# all_pvals["LC-ME","EqO2_max"]<-a$p.adj[3]
all_pvals["BED_REST_test","EqO2_max"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","EqO2_max"]<-"anova_no_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55","p<0.001","BED REST",
  "CON","ME",paste(round(a$`Pr(>Chisq)`,3)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

EqO2_max_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=EqCO2_max, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("V\U0307"~E[max]*"/V\U0307"~O[2][max]*""))+
  scale_y_continuous(limits=c(0,90), breaks=seq(0,80,20),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(50), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(7),
              space=c(0.5),
              scales=c(20),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05))
  )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=EqO2_max_a,
       filename = paste(output_folder,"/EqO2_max_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

EqO2_max_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=EqO2_max, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("V\U0307"~E[max]*"/V\U0307"~O[2][max]*""))+
  scale_y_continuous(limits=c(0,90), breaks=seq(0,80,20),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(85), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(7),
              space=c(0.5),
              scales=c(20),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=EqO2_max_b,
       filename = paste(output_folder,"/EqO2_max_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

# VE/VCO2 slope -----------------------------------------------------------

# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","VE_VCO2_slope"],na.rm=TRUE); sd(data[data$Session=="BDC","VE_VCO2_slope"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","VE_VCO2_slope"],na.rm=TRUE);sd(data[data$Session=="HDT55","VE_VCO2_slope"],na.rm=TRUE)
# mean(data[data$Session=="CON","VE_VCO2_slope"],na.rm=TRUE);sd(data[data$Session=="CON","VE_VCO2_slope"],na.rm=TRUE)
# mean(data[data$Session=="LC","VE_VCO2_slope"],na.rm=TRUE);sd(data[data$Session=="LC","VE_VCO2_slope"],na.rm=TRUE)
# mean(data[data$Session=="ME","VE_VCO2_slope"],na.rm=TRUE);sd(data[data$Session=="ME","VE_VCO2_slope"],na.rm=TRUE)


remove<-data[is.na(data$VE_VCO2_slope)==TRUE,"Subject"]
remove2<-c("P110117","P110118","P110127")
test2<-test[!test$Subject %in% remove & !test$Subject %in% remove2,]

# qqnorm(test2[test2$Session=="BDC","VE_VCO2_slope"]);qqline(test2[test2$Session=="BDC","VE_VCO2_slope"])
# qqnorm(test2[test2$Session=="HDT55","VE_VCO2_slope"]);qqline(test2[test2$Session=="HDT55","VE_VCO2_slope"])
# qqnorm(test2[test2$Session=="CON","VE_VCO2_slope"]);qqline(test2[test2$Session=="CON","VE_VCO2_slope"])
# qqnorm(test2[test2$Session=="LC","VE_VCO2_slope"]);qqline(test2[test2$Session=="LC","VE_VCO2_slope"])
# qqnorm(test2[test2$Session=="ME","VE_VCO2_slope"]);qqline(test2[test2$Session=="ME","VE_VCO2_slope"])


model <-lme(VE_VCO2_slope~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")

a<-car::Anova(model)
all_pvals["ANOVA","VE_VCO2_slope"]<-a$`Pr(>Chisq)`

b<-
  t_test(test2[test2$Group=="BED REST",], VE_VCO2_slope~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","VE_VCO2_slope"]<-b$p
# all_pvals["CON-ME","VE_VCO2_slope"]<-a$p.adj[2]
# all_pvals["CON-LC","VE_VCO2_slope"]<-a$p.adj[1]
# all_pvals["LC-ME","VE_VCO2_slope"]<-a$p.adj[3]
all_pvals["BED_REST_test","VE_VCO2_slope"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","VE_VCO2_slope"]<-"anova_no_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55","p<0.001","BED REST",
  "CON","ME",paste(format(round(a$`Pr(>Chisq)`,3),nsmall=3)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove & !data$Subject %in% remove2,]

VE_VCO2_slope_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=VE_VCO2_slope, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("V\U0307 E-V\U0307"~CO[2]*"slope"))+
  scale_y_continuous(limits=c(0,70), breaks=seq(0,60,20),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(58), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=VE_VCO2_slope_a,
       filename = paste(output_folder,"/VE_VCO2_slope_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

VE_VCO2_slope_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=VE_VCO2_slope, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("V\U0307 E-V\U0307"~CO[2]*"slope"))+
  scale_y_continuous(limits=c(0,70), breaks=seq(0,60,20),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(50), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=VE_VCO2_slope_b,
       filename = paste(output_folder,"/VE_VCO2_slope_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")
# O2 Pulse ----------------------------------------------------------------


# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","O2_pulse"],na.rm=TRUE); sd(data[data$Session=="BDC","O2_pulse"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","O2_pulse"],na.rm=TRUE);sd(data[data$Session=="HDT55","O2_pulse"],na.rm=TRUE)
# mean(data[data$Session=="CON","O2_pulse"],na.rm=TRUE);sd(data[data$Session=="CON","O2_pulse"],na.rm=TRUE)
# mean(data[data$Session=="LC","O2_pulse"],na.rm=TRUE);sd(data[data$Session=="LC","O2_pulse"],na.rm=TRUE)
# mean(data[data$Session=="ME","O2_pulse"],na.rm=TRUE);sd(data[data$Session=="ME","O2_pulse"],na.rm=TRUE)


remove<-data[is.na(data$O2_pulse)==TRUE,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","O2_pulse"]);qqline(test2[test2$Session=="BDC","O2_pulse"])
# qqnorm(test2[test2$Session=="HDT55","O2_pulse"]);qqline(test2[test2$Session=="HDT55","O2_pulse"])
# qqnorm(test2[test2$Session=="CON","O2_pulse"]);qqline(test2[test2$Session=="CON","O2_pulse"])
# qqnorm(test2[test2$Session=="LC","O2_pulse"]);qqline(test2[test2$Session=="LC","O2_pulse"])
# qqnorm(test2[test2$Session=="ME","O2_pulse"]);qqline(test2[test2$Session=="ME","O2_pulse"])

model <-lme(O2_pulse~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
all_pvals["ANOVA","O2_pulse"]<-a$`Pr(>Chisq)`

# post-hoc testing 
a<-
  test2[test2$Group=="POST-VIRAL",] %>% tukey_hsd(O2_pulse~Session)

b<-
  wilcox.test(test2[test2$Session=="BDC","O2_pulse"], test2[test2$Session=="HDT55","O2_pulse"],paired=TRUE)

all_pvals["BDC-HDT55","O2_pulse"]<-b$p.value
all_pvals["CON-ME","O2_pulse"]<-a$p.adj[2]
all_pvals["CON-LC","O2_pulse"]<-a$p.adj[1]
all_pvals["LC-ME","O2_pulse"]<-a$p.adj[3]
all_pvals["BED_REST_test","O2_pulse"]<-"wilcoxon"
all_pvals["POST_VIRAL_test","O2_pulse"]<-"anova_tukey_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55",paste("p<0.001"),"BED REST",
  "CON","ME",paste("p<0.001"),"POST-VIRAL",
  "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

O2_pulse_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=O2_pulse, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(~O[2]*"-pulse (L"~beat^-1*")"))+
  scale_y_continuous(limits=c(0,30), breaks=seq(0,25,5),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(25), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))




ggsave(plot=O2_pulse_a,
       filename = paste(output_folder,"/O2_pulse_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

O2_pulse_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=O2_pulse, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(~O[2]*"-pulse (L"~beat^-1*")"))+
  scale_y_continuous(limits=c(0,30), breaks=seq(0,25,5),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(27), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(25), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[4,], label = "p.adj",y.position=c(23), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=O2_pulse_b,
       filename = paste(output_folder,"/O2_pulse_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")
# Max HR ------------------------------------------------------------------

# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","HR_max"],na.rm=TRUE); sd(data[data$Session=="BDC","HR_max"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","HR_max"],na.rm=TRUE);sd(data[data$Session=="HDT55","HR_max"],na.rm=TRUE)
# mean(data[data$Session=="CON","HR_max"],na.rm=TRUE);sd(data[data$Session=="CON","HR_max"],na.rm=TRUE)
# mean(data[data$Session=="LC","HR_max"],na.rm=TRUE);sd(data[data$Session=="LC","HR_max"],na.rm=TRUE)
# mean(data[data$Session=="ME","HR_max"],na.rm=TRUE);sd(data[data$Session=="ME","HR_max"],na.rm=TRUE)


remove<-data[is.na(data$HR_max)==TRUE,"Subject"]
test2<-test[!test$Subject %in% remove ,]
# 
# qqnorm(test2[test2$Session=="BDC","HR_max"]);qqline(test2[test2$Session=="BDC","HR_max"])
# qqnorm(test2[test2$Session=="HDT55","HR_max"]);qqline(test2[test2$Session=="HDT55","HR_max"])
# qqnorm(test2[test2$Session=="CON","HR_max"]);qqline(test2[test2$Session=="CON","HR_max"])
# qqnorm(test2[test2$Session=="LC","HR_max"]);qqline(test2[test2$Session=="LC","HR_max"])
# qqnorm(test2[test2$Session=="ME","HR_max"]);qqline(test2[test2$Session=="ME","HR_max"])

model <-lme(HR_max~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
all_pvals["ANOVA","HR_max"]<-a$`Pr(>Chisq)`


b<-
  wilcox.test(test2[test2$Session=="BDC","HR_max"], test2[test2$Session=="HDT55","HR_max"],paired=TRUE)

all_pvals["BDC-HDT55","HR_max"]<-b$p.value
# all_pvals["CON-ME","HR_max"]<-a$p.adj[2]
# all_pvals["CON-LC","HR_max"]<-a$p.adj[1]
# all_pvals["LC-ME","HR_max"]<-a$p.adj[3]
all_pvals["BED_REST_test","HR_max"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","HR_max"]<-"anova_no_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55","p<0.001","BED REST",
  "CON","ME",paste(round(a$`Pr(>Chisq)`,3)),"POST-VIRAL")
# "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
# "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

HR_max_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=HR_max, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(atop("Maximum Heart Rate", "(Beats"*~min^-1*")")))+
  scale_y_continuous(limits=c(0,225), breaks=c(0,120,150,180,210),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(210), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(90),
              space=c(0.5),
              scales=c(30),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))



ggsave(plot=HR_max_a,
       filename = paste(output_folder,"/HR_max_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

HR_max_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=HR_max, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(atop("Maximum Heart Rate", "(Beats"*~min^-1*")")))+
  scale_y_continuous(limits=c(0,225), breaks=c(0,120,150,180,210),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(220), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(90),
              space=c(0.5),
              scales=c(30),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=HR_max_b,
       filename = paste(output_folder,"/HR_max_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")


# Adjusted Heart Rate Reserve --------------------------------------------------------------------

# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","AHRR"],na.rm=TRUE); sd(data[data$Session=="BDC","AHRR"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","AHRR"],na.rm=TRUE);sd(data[data$Session=="HDT55","AHRR"],na.rm=TRUE)
# mean(data[data$Session=="CON","AHRR"],na.rm=TRUE);sd(data[data$Session=="CON","AHRR"],na.rm=TRUE)
# mean(data[data$Session=="LC","AHRR"],na.rm=TRUE);sd(data[data$Session=="LC","AHRR"],na.rm=TRUE)
# mean(data[data$Session=="ME","AHRR"],na.rm=TRUE);sd(data[data$Session=="ME","AHRR"],na.rm=TRUE)


remove<-data[is.na(data$AHRR)==TRUE,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","AHRR"]);qqline(test2[test2$Session=="BDC","AHRR"])
# qqnorm(test2[test2$Session=="HDT55","AHRR"]);qqline(test2[test2$Session=="HDT55","AHRR"])
# qqnorm(test2[test2$Session=="CON","AHRR"]);qqline(test2[test2$Session=="CON","AHRR"])
# qqnorm(test2[test2$Session=="LC","AHRR"]);qqline(test2[test2$Session=="LC","AHRR"])
# qqnorm(test2[test2$Session=="ME","AHRR"]);qqline(test2[test2$Session=="ME","AHRR"])


# data2<-data[!data$Subject %in% remove,]
# length(data2[data2$Session=="BDC" & data2$AHRR<0.8,"Subject"])
# length(data2[data2$Session=="HDT55" & data2$AHRR<0.8,"Subject"])
# length(data2[data2$Session=="CON" & data2$AHRR<0.8,"Subject"])
# length(data2[data2$Session=="LC" & data2$AHRR<0.8,"Subject"])
# length(data2[data2$Session=="ME" & data2$AHRR<0.8,"Subject"])



a<-
  kruskal.test(AHRR~Session, data=test2[test2$Group=="POST-VIRAL",])
all_pvals["ANOVA","AHRR"]<-a$p.value

b<-
  wilcox.test(test2[test2$Session=="BDC","AHRR"], test2[test2$Session=="HDT55","AHRR"],paired=TRUE)

all_pvals["BDC-HDT55","AHRR"]<-b$p.value
# all_pvals["CON-ME","AHRR"]<-a$p.value[2]
# all_pvals["CON-LC","AHRR"]<-a$p.value[1]
# all_pvals["LC-ME","AHRR"]<-a$p.value[4]
all_pvals["BED_REST_test","AHRR"]<-"wilcoxon"
all_pvals["POST_VIRAL_test","AHRR"]<-"kruskal_no_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55",paste(format(round(b$p.value,3), drop0trailing=F)),"BED REST",
  "CON","ME",paste(round(a$p.value,3)),"POST-VIRAL")
# "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
# "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

AHRR_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=AHRR*100, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  geom_hline(yintercept=80, linetype ="dashed", colour="black", size=line_size/2)+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(atop("Adjusted Heart","Rate Reserve (%)")))+
  scale_y_continuous(limits=c(0,190), breaks=seq(0,175,25),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(175), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=AHRR_a,
       filename = paste(output_folder,"/AHRR_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

AHRR_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=AHRR*100, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  geom_hline(yintercept=80, linetype ="dashed", colour="black", size=line_size/2)+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(atop("Adjusted Heart","Rate Reserve (%)")))+
  scale_y_continuous(limits=c(0,190), breaks=seq(0,175,25),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(180), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=AHRR_b,
       filename = paste(output_folder,"/AHRR_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

# VO2-Heart rate slope ----------------------------------------------------
# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","VO2_HR_slope"],na.rm=TRUE); sd(data[data$Session=="BDC","VO2_HR_slope"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","VO2_HR_slope"],na.rm=TRUE);sd(data[data$Session=="HDT55","VO2_HR_slope"],na.rm=TRUE)
# mean(data[data$Session=="CON","VO2_HR_slope"],na.rm=TRUE);sd(data[data$Session=="CON","VO2_HR_slope"],na.rm=TRUE)
# mean(data[data$Session=="LC","VO2_HR_slope"],na.rm=TRUE);sd(data[data$Session=="LC","VO2_HR_slope"],na.rm=TRUE)
# mean(data[data$Session=="ME","VO2_HR_slope"],na.rm=TRUE);sd(data[data$Session=="ME","VO2_HR_slope"],na.rm=TRUE)


remove<-data[is.na(data$VO2_HR_slope)==TRUE,"Subject"]
remove2<-c("P110117","P110118","P110127")
test2<-test[!test$Subject %in% remove & !test$Subject %in% remove2,]

# qqnorm(test2[test2$Session=="BDC","VO2_HR_slope"]);qqline(test2[test2$Session=="BDC","VO2_HR_slope"])
# qqnorm(test2[test2$Session=="HDT55","VO2_HR_slope"]);qqline(test2[test2$Session=="HDT55","VO2_HR_slope"])
# qqnorm(test2[test2$Session=="CON","VO2_HR_slope"]);qqline(test2[test2$Session=="CON","VO2_HR_slope"])
# qqnorm(test2[test2$Session=="LC","VO2_HR_slope"]);qqline(test2[test2$Session=="LC","VO2_HR_slope"])
# qqnorm(test2[test2$Session=="ME","VO2_HR_slope"]);qqline(test2[test2$Session=="ME","VO2_HR_slope"])

model <-lme(VO2_HR_slope~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")

a<-car::Anova(model)
all_pvals["ANOVA","VO2_HR_slope"]<-a$`Pr(>Chisq)`
# post-hoc testing 
a<-
  test2[test2$Group=="POST-VIRAL",] %>% tukey_hsd(VO2_HR_slope~Session)

b<-
  t_test(test2[test2$Group=="BED REST",], VO2_HR_slope~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","VO2_HR_slope"]<-b$p
all_pvals["CON-ME","VO2_HR_slope"]<-a$p.adj[2]
all_pvals["CON-LC","VO2_HR_slope"]<-a$p.adj[1]
all_pvals["LC-ME","VO2_HR_slope"]<-a$p.adj[3]
all_pvals["BED_REST_test","VO2_HR_slope"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","VO2_HR_slope"]<-"anova_tukey_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55",paste(format(round(b$p,3), drop0trailing=F)),"BED REST",
  "CON","ME",paste(format(round(a$p.adj[2],3), drop0trailing=F)),"POST-VIRAL",
  "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

VO2_HR_slope_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=VO2_HR_slope, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("V\U0307" ~O[2]*"-Heart Rate slope"))+
  scale_y_continuous(limits=c(0,110), breaks=seq(0,90,30),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(60), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=VO2_HR_slope_a,
       filename = paste(output_folder,"/VO2_HR_slope_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

VO2_HR_slope_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=VO2_HR_slope, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("V\U0307" ~O[2]*"-Heart Rate slope"))+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  scale_y_continuous(limits=c(0,110), breaks=seq(0,90,30),expand=c(0,0))+
  scale_y_continuous(limits=c(0,120), breaks=seq(0,90,30),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(100), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(95), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[4,], label = "p.adj",y.position=c(90), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=VO2_HR_slope_b,
       filename = paste(output_folder,"/VO2_HR_slope_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

# VE VCO2 slope submax -------------------------------------------------------

# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","VE_VCO2_slope_submax"],na.rm=TRUE); sd(data[data$Session=="BDC","VE_VCO2_slope_submax"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","VE_VCO2_slope_submax"],na.rm=TRUE);sd(data[data$Session=="HDT55","VE_VCO2_slope_submax"],na.rm=TRUE)
# mean(data[data$Session=="CON","VE_VCO2_slope_submax"],na.rm=TRUE);sd(data[data$Session=="CON","VE_VCO2_slope_submax"],na.rm=TRUE)
# mean(data[data$Session=="LC","VE_VCO2_slope_submax"],na.rm=TRUE);sd(data[data$Session=="LC","VE_VCO2_slope_submax"],na.rm=TRUE)
# mean(data[data$Session=="ME","VE_VCO2_slope_submax"],na.rm=TRUE);sd(data[data$Session=="ME","VE_VCO2_slope_submax"],na.rm=TRUE)


remove<-data[is.na(data$VE_VCO2_slope_submax)==TRUE,"Subject"]
remove2<-c("P110117","P110118","P110127")
test2<-test[!test$Subject %in% remove & !test$Subject %in% remove2,]

# qqnorm(test2[test2$Session=="BDC","VE_VCO2_slope_submax"]);qqline(test2[test2$Session=="BDC","VE_VCO2_slope_submax"])
# qqnorm(test2[test2$Session=="HDT55","VE_VCO2_slope_submax"]);qqline(test2[test2$Session=="HDT55","VE_VCO2_slope_submax"])
# qqnorm(test2[test2$Session=="CON","VE_VCO2_slope_submax"]);qqline(test2[test2$Session=="CON","VE_VCO2_slope_submax"])
# qqnorm(test2[test2$Session=="LC","VE_VCO2_slope_submax"]);qqline(test2[test2$Session=="LC","VE_VCO2_slope_submax"])
# qqnorm(test2[test2$Session=="ME","VE_VCO2_slope_submax"]);qqline(test2[test2$Session=="ME","VE_VCO2_slope_submax"])


model <-lme(VE_VCO2_slope_submax~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")

a<-car::Anova(model)
all_pvals["ANOVA","VE_VCO2_slope_submax"]<-a$`Pr(>Chisq)`
# post-hoc testing 
a<-
  test2[test2$Group=="POST-VIRAL",] %>% tukey_hsd(VE_VCO2_slope_submax~Session)
       
b<-
  t_test(test2[test2$Group=="BED REST",], VE_VCO2_slope_submax~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","VE_VCO2_slope_submax"]<-b$p
all_pvals["CON-ME","VE_VCO2_slope_submax"]<-a$p.adj[2]
all_pvals["CON-LC","VE_VCO2_slope_submax"]<-a$p.adj[1]
all_pvals["LC-ME","VE_VCO2_slope_submax"]<-a$p.adj[3]
all_pvals["BED_REST_test","VE_VCO2_slope_submax"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","VE_VCO2_slope_submax"]<-"anova_tukey_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55","p<0.001","BED REST",
  "CON","ME",paste(format(round(a$p.adj[2],3), drop0trailing=F)),"POST-VIRAL",
  "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.adj[3],3), nsmall=2)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

VE_VCO2_slope_submax_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=VE_VCO2_slope_submax, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(atop("V\U0307 E-V\U0307"~CO[2]*"slope", "(Submaximal)")))+
  scale_y_continuous(limits=c(0,70), breaks=seq(0,60,20),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(58), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=VE_VCO2_slope_submax_a,
       filename = paste(output_folder,"/VE_VCO2_slope_submax_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

VE_VCO2_slope_submax_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=VE_VCO2_slope_submax, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(atop("V\U0307 E-V\U0307"~CO[2]*"slope", "(Submaximal)")))+
  scale_y_continuous(limits=c(0,70), breaks=seq(0,60,20),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(60), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(55), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[4,], label = "p.adj",y.position=c(50), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=VE_VCO2_slope_submax_b,
       filename = paste(output_folder,"/VE_VCO2_slope_submax_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

# resting heart rate -------------------------------------------------------

# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","HR_rest"],na.rm=TRUE); sd(data[data$Session=="BDC","HR_rest"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","HR_rest"],na.rm=TRUE);sd(data[data$Session=="HDT55","HR_rest"],na.rm=TRUE)
# mean(data[data$Session=="CON","HR_rest"],na.rm=TRUE);sd(data[data$Session=="CON","HR_rest"],na.rm=TRUE)
# mean(data[data$Session=="LC","HR_rest"],na.rm=TRUE);sd(data[data$Session=="LC","HR_rest"],na.rm=TRUE)
# mean(data[data$Session=="ME","HR_rest"],na.rm=TRUE);sd(data[data$Session=="ME","HR_rest"],na.rm=TRUE)


remove<-data[is.na(data$HR_rest)==TRUE,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","HR_rest"]);qqline(test2[test2$Session=="BDC","HR_rest"])
# qqnorm(test2[test2$Session=="HDT55","HR_rest"]);qqline(test2[test2$Session=="HDT55","HR_rest"])
# qqnorm(test2[test2$Session=="CON","HR_rest"]);qqline(test2[test2$Session=="CON","HR_rest"])
# qqnorm(test2[test2$Session=="LC","HR_rest"]);qqline(test2[test2$Session=="LC","HR_rest"])
# qqnorm(test2[test2$Session=="ME","HR_rest"]);qqline(test2[test2$Session=="ME","HR_rest"])

model <-lme(HR_rest~ Session, data=test2[test2$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
all_pvals["ANOVA","HR_rest"]<-a$`Pr(>Chisq)`
# post-hoc testing 
a<-
  test2[test2$Group=="POST-VIRAL",] %>% tukey_hsd(HR_rest~Session)


b<-
  wilcox.test(test2[test2$Session=="BDC","HR_rest"], test2[test2$Session=="HDT55","HR_rest"],paired=TRUE)

all_pvals["BDC-HDT55","HR_rest"]<-b$p.value
all_pvals["CON-ME","HR_rest"]<-a$p.adj[2]
all_pvals["CON-LC","HR_rest"]<-a$p.adj[1]
all_pvals["LC-ME","HR_rest"]<-a$p.adj[3]
all_pvals["BED_REST_test","HR_rest"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","HR_rest"]<-"anova_no_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55","p<0.001","BED REST",
  "CON","ME",paste(format(round(a$p.adj[2],3), drop0trailing=F)),"POST-VIRAL",
  "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

HR_rest_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=HR_rest, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(atop("Baseline Heart Rate", "(Beats"*~min^-1*")")))+
  scale_y_continuous(limits=c(0,170), breaks=c(0,50,75,100,125,150),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(140), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(35),
              space=c(0.5),
              scales=c(30),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))



ggsave(plot=HR_rest_a,
       filename = paste(output_folder,"/HR_rest_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

HR_rest_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=HR_rest, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(atop("Baseline Heart Rate", "(Beats"*~min^-1*")")))+
  scale_y_continuous(limits=c(0,170), breaks=c(0,50,75,100,125,150),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(160), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(150), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[4,], label = "p.adj",y.position=c(140), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(35),
              space=c(0.5),
              scales=c(30),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=HR_rest_b,
       filename = paste(output_folder,"/HR_rest_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")




# Supplemental 4 ----------------------------------------------------------


# Oxphos normalized to SDH ------------------------------------------------


# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","Oxphos_norm"],na.rm=TRUE); sd(data[data$Session=="BDC","Oxphos_norm"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","Oxphos_norm"],na.rm=TRUE);sd(data[data$Session=="HDT55","Oxphos_norm"],na.rm=TRUE)
# mean(data[data$Session=="CON","Oxphos_norm"],na.rm=TRUE);sd(data[data$Session=="CON","Oxphos_norm"],na.rm=TRUE)
# mean(data[data$Session=="LC","Oxphos_norm"],na.rm=TRUE);sd(data[data$Session=="LC","Oxphos_norm"],na.rm=TRUE)
# mean(data[data$Session=="ME","Oxphos_norm"],na.rm=TRUE);sd(data[data$Session=="ME","Oxphos_norm"],na.rm=TRUE)


remove<-data[is.na(data$Oxphos_norm)==TRUE | data$membrane_intact>1.1 ,"Subject"]
resp_data<-data[!data$Subject %in% remove, ]
resp_test<-as.data.frame(box_cox_transform(resp_data))
resp_test[,c("Subject","Session","Group", "Sex")]<-resp_data[,c("Subject","Session","Group","Sex")]
resp_test_norm<-normality_test(resp_test)
test2<-resp_test[!resp_test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","Oxphos_norm"]);qqline(test2[test2$Session=="BDC","Oxphos_norm"])
# qqnorm(test2[test2$Session=="HDT55","Oxphos_norm"]);qqline(test2[test2$Session=="HDT55","Oxphos_norm"])
# qqnorm(test2[test2$Session=="CON","Oxphos_norm"]);qqline(test2[test2$Session=="CON","Oxphos_norm"])
# qqnorm(test2[test2$Session=="LC","Oxphos_norm"]);qqline(test2[test2$Session=="LC","Oxphos_norm"])
# qqnorm(test2[test2$Session=="ME","Oxphos_norm"]);qqline(test2[test2$Session=="ME","Oxphos_norm"])


a<-kruskal.test(Oxphos_norm~Session, test2[test2$Group=="POST-VIRAL",])
all_pvals["ANOVA","Oxphos_norm"]<-a$p.value
#post hoc testing
a<-pairwise.wilcox.test(test2[test2$Group=="POST-VIRAL","Oxphos_norm"], test2[test2$Group=="POST-VIRAL","Session"], p.adjust.method="BH")


b<-
  wilcox.test(test2[test2$Session=="BDC","Oxphos_norm"],test2[test2$Session=="HDT55","Oxphos_norm"],paired=TRUE)


all_pvals["BDC-HDT55","Oxphos_norm"]<-b$p.value
all_pvals["CON-ME","Oxphos_norm"]<-a$p.value[2]
all_pvals["CON-LC","Oxphos_norm"]<-a$p.value[1]
all_pvals["LC-ME","Oxphos_norm"]<-a$p.value[4]
all_pvals["BED_REST_test","Oxphos_norm"]<-"wilcoxon"
all_pvals["POST_VIRAL_test","Oxphos_norm"]<-"kruskal_wilcoxon_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55",paste(format(round(b$p.value,3), drop0trailing=F)),"BED REST",
  "CON","ME",paste(format(round(a$p.value[2],3), drop0trailing=F)),"POST-VIRAL",
  "CON","LC",paste(format(round(a$p.value[1],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.value[4],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

Oxphos_norm_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=Oxphos_norm/10^6, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(bquote(atop("Oxidative Phosphorylation","(pmol"*~s^-1*~mg^-1*~SDH^-1*~10^-6*")")))+
  scale_y_continuous(limits=c(0,21), breaks=seq(0,20,5),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(15), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=Oxphos_norm_a,
       filename = paste(output_folder,"/Oxphos_norm_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")


Oxphos_norm_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=Oxphos_norm/10^6, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(bquote(atop("Oxidative Phosphorylation","(pmol"*~s^-1*~mg^-1*~SDH^-1*~10^-6*")")))+
  scale_y_continuous(limits=c(0,21), breaks=seq(0,21,5),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(20), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(17), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[4,], label = "p.adj",y.position=c(12), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=Oxphos_norm_b,
       filename = paste(output_folder,"/Oxphos_norm_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")
# E/L coupling efficiency -------------------------------------------------


# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","E_L_coup"],na.rm=TRUE); sd(data[data$Session=="BDC","E_L_coup"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","E_L_coup"],na.rm=TRUE);sd(data[data$Session=="HDT55","E_L_coup"],na.rm=TRUE)
# mean(data[data$Session=="CON","E_L_coup"],na.rm=TRUE);sd(data[data$Session=="CON","E_L_coup"],na.rm=TRUE)
# mean(data[data$Session=="LC","E_L_coup"],na.rm=TRUE);sd(data[data$Session=="LC","E_L_coup"],na.rm=TRUE)
# mean(data[data$Session=="ME","E_L_coup"],na.rm=TRUE);sd(data[data$Session=="ME","E_L_coup"],na.rm=TRUE)


remove<-data[is.na(data$E_L_coup)==TRUE | data$membrane_intact>1.1 ,"Subject"]
resp_data<-data[!data$Subject %in% remove, ]
resp_test<-as.data.frame(box_cox_transform(resp_data))
resp_test[,c("Subject","Session","Group", "Sex")]<-resp_data[,c("Subject","Session","Group","Sex")]
resp_test_norm<-normality_test(resp_test)
test2<-resp_test[!resp_test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","E_L_coup"]);qqline(test2[test2$Session=="BDC","E_L_coup"])
# qqnorm(test2[test2$Session=="HDT55","E_L_coup"]);qqline(test2[test2$Session=="HDT55","E_L_coup"])
# qqnorm(test2[test2$Session=="CON","E_L_coup"]);qqline(test2[test2$Session=="CON","E_L_coup"])
# qqnorm(test2[test2$Session=="LC","E_L_coup"]);qqline(test2[test2$Session=="LC","E_L_coup"])
# qqnorm(test2[test2$Session=="ME","E_L_coup"]);qqline(test2[test2$Session=="ME","E_L_coup"])


a<-
  kruskal.test(E_L_coup~Session, test2[test2$Group=="POST-VIRAL",])
all_pvals["ANOVA","E_L_coup"]<-a$p.value


b<-
  wilcox.test(test2[test2$Session=="BDC","E_L_coup"],test2[test2$Session=="HDT55","E_L_coup"],paired=TRUE)

all_pvals["BDC-HDT55","E_L_coup"]<-b$p.value
# all_pvals["CON-ME","E_L_coup"]<-a$p.value[2]
# all_pvals["CON-LC","E_L_coup"]<-a$p.value[1]
# all_pvals["LC-ME","E_L_coup"]<-a$p.value[4]
all_pvals["BED_REST_test","E_L_coup"]<-"wilcoxon"
all_pvals["POST_VIRAL_test","E_L_coup"]<-"kruskal_no_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55",paste(format(round(b$p.value,3), drop0trailing=F)),"BED REST",
  "CON","ME",paste(round(a$p.value,3)),"POST-VIRAL")
# "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
# "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

E_L_coupling_efficiency_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=E_L_coup*100, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(bquote(atop("E/L Coupling Efficiency (AU)")))+
  scale_y_continuous(limits=c(0,110), breaks=c(0,60,70,80,90,100),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(100), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(50),
              space=c(1),
              scales=c(50),
              which=c(1),
              expand = expansion(mult = c(0,0)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))



ggsave(plot=E_L_coupling_efficiency_a,
       filename = paste(output_folder,"/E_L_coupling_efficiency_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")


E_L_coupling_efficiency_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=E_L_coup*100, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(bquote(atop("E/L Coupling Efficiency (AU)")))+
  scale_y_continuous(limits=c(0,110), breaks=c(0,60,70,80,90,100),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(100), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(50),
              space=c(1),
              scales=c(50),
              which=c(1),
              expand = expansion(mult = c(0,0)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=E_L_coupling_efficiency_b,
       filename = paste(output_folder,"/E_L_coupling_efficiency_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

# PS_PNS ------------------------------------------------------------------


# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","PS_PNS"],na.rm=TRUE); sd(data[data$Session=="BDC","PS_PNS"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","PS_PNS"],na.rm=TRUE);sd(data[data$Session=="HDT55","PS_PNS"],na.rm=TRUE)
# mean(data[data$Session=="CON","PS_PNS"],na.rm=TRUE);sd(data[data$Session=="CON","PS_PNS"],na.rm=TRUE)
# mean(data[data$Session=="LC","PS_PNS"],na.rm=TRUE);sd(data[data$Session=="LC","PS_PNS"],na.rm=TRUE)
# mean(data[data$Session=="ME","PS_PNS"],na.rm=TRUE);sd(data[data$Session=="ME","PS_PNS"],na.rm=TRUE)



remove<-data[is.na(data$PS_PNS)==TRUE | data$membrane_intact>1.1 ,"Subject"]
resp_data<-data[!data$Subject %in% remove, ]
resp_test<-as.data.frame(box_cox_transform(resp_data))
resp_test[,c("Subject","Session","Group", "Sex")]<-resp_data[,c("Subject","Session","Group","Sex")]
resp_test_norm<-normality_test(resp_test)
test2<-resp_test[!resp_test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","PS_PNS"]);qqline(test2[test2$Session=="BDC","PS_PNS"])
# qqnorm(test2[test2$Session=="HDT55","PS_PNS"]);qqline(test2[test2$Session=="HDT55","PS_PNS"])
# qqnorm(test2[test2$Session=="CON","PS_PNS"]);qqline(test2[test2$Session=="CON","PS_PNS"])
# qqnorm(test2[test2$Session=="LC","PS_PNS"]);qqline(test2[test2$Session=="LC","PS_PNS"])
# qqnorm(test2[test2$Session=="ME","PS_PNS"]);qqline(test2[test2$Session=="ME","PS_PNS"])

a<-
  kruskal.test(PS_PNS~Session, test2[test2$Group=="POST-VIRAL",])
all_pvals["ANOVA","PS_PNS"]<-a$p.value


b<-
  wilcox.test(test2[test2$Session=="BDC","PS_PNS"],test2[test2$Session=="HDT55","PS_PNS"],paired=TRUE)

all_pvals["BDC-HDT55","PS_PNS"]<-b$p.value
# all_pvals["CON-ME","PS_PNS"]<-a$p.value[2]
# all_pvals["CON-LC","PS_PNS"]<-a$p.value[1]
# all_pvals["LC-ME","PS_PNS"]<-a$p.value[4]
all_pvals["BED_REST_test","PS_PNS"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","PS_PNS"]<-"kruskal_no_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55",paste(format(round(b$p.value,3), drop0trailing=F)),"BED REST",
  "CON","ME",paste(round(a$p.value,3)),"POST-VIRAL")
# "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
# "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

PS_PNS_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=PS_PNS, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(bquote("PS/PNS (AU)"))+
  scale_y_continuous(limits=c(0,110), breaks=c(0,20,40,60,80,100),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(80), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(15),
              space=c(1),
              scales=c(15),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        # axis.ticks = element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))



ggsave(plot=PS_PNS_a,
       filename = paste(output_folder,"/PS_PNS_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")


PS_PNS_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=PS_PNS, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(bquote("PS/PNS (AU)"))+
  scale_y_continuous(limits=c(0,110), breaks=c(0,20,40,60,80,100),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(100), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(15),
              space=c(1),
              scales=c(15),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=PS_PNS_b,
       filename = paste(output_folder,"/PS_PNS_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")




# PN_PNS ------------------------------------------------------------------


# Stats -------------------------------------------------------------------

# mean(data[data$Session=="BDC","PN_PNS"],na.rm=TRUE); sd(data[data$Session=="BDC","PN_PNS"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","PN_PNS"],na.rm=TRUE);sd(data[data$Session=="HDT55","PN_PNS"],na.rm=TRUE)
# mean(data[data$Session=="CON","PN_PNS"],na.rm=TRUE);sd(data[data$Session=="CON","PN_PNS"],na.rm=TRUE)
# mean(data[data$Session=="LC","PN_PNS"],na.rm=TRUE);sd(data[data$Session=="LC","PN_PNS"],na.rm=TRUE)
# mean(data[data$Session=="ME","PN_PNS"],na.rm=TRUE);sd(data[data$Session=="ME","PN_PNS"],na.rm=TRUE)


remove<-data[is.na(data$PN_PNS)==TRUE | data$membrane_intact>1.1 ,"Subject"]
resp_data<-data[!data$Subject %in% remove, ]
resp_test<-as.data.frame(box_cox_transform(resp_data))
resp_test[,c("Subject","Session","Group", "Sex")]<-resp_data[,c("Subject","Session","Group","Sex")]
resp_test_norm<-normality_test(resp_test)
test2<-resp_test[!resp_test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","PN_PNS"]);qqline(test2[test2$Session=="BDC","PN_PNS"])
# qqnorm(test2[test2$Session=="HDT55","PN_PNS"]);qqline(test2[test2$Session=="HDT55","PN_PNS"])
# qqnorm(test2[test2$Session=="CON","PN_PNS"]);qqline(test2[test2$Session=="CON","PN_PNS"])
# qqnorm(test2[test2$Session=="LC","PN_PNS"]);qqline(test2[test2$Session=="LC","PN_PNS"])
# qqnorm(test2[test2$Session=="ME","PN_PNS"]);qqline(test2[test2$Session=="ME","PN_PNS"])

a<-
  kruskal.test(PN_PNS~Session, test2[test2$Group=="POST-VIRAL",])
all_pvals["ANOVA","PN_PNS"]<-a$p.value

b<-
  wilcox.test(test2[test2$Session=="BDC","PN_PNS"],test2[test2$Session=="HDT55","PN_PNS"],paired=TRUE)

all_pvals["BDC-HDT55","PN_PNS"]<-b$p.value
# all_pvals["CON-ME","PN_PNS"]<-a$p.value[2]
# all_pvals["CON-LC","PN_PNS"]<-a$p.value[1]
# all_pvals["LC-ME","PN_PNS"]<-a$p.value[4]
all_pvals["BED_REST_test","PN_PNS"]<-"wilcoxon"
all_pvals["POST_VIRAL_test","PN_PNS"]<-"kruskal_no_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55",paste(format(round(b$p.value,3), drop0trailing=F)),"BED REST",
  "CON","ME",paste(round(a$p.value,3)),"POST-VIRAL")
# "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
# "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

PN_PNS_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=PN_PNS, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(bquote(atop("PN/PNS (AU)")))+
  scale_y_continuous(limits=c(0,110), breaks=seq(0,100,20),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(80), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))



ggsave(plot=PN_PNS_a,
       filename = paste(output_folder,"/PN_PNS_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")


PN_PNS_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=PN_PNS, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(bquote(atop("PN/PNS (AU)")))+
  scale_y_continuous(limits=c(0,110), breaks=seq(0,100,20),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(100), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=PN_PNS_b,
       filename = paste(output_folder,"/PN_PNS_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")


# Supplemental 5 ----------------------------------------------------------

# capillary size ----------------------------------------------------------


# Stats -------------------------------------------------------------------
# 
# mean(data[data$Session=="BDC","Cap_size"],na.rm=TRUE); sd(data[data$Session=="BDC","Cap_size"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","Cap_size"],na.rm=TRUE);sd(data[data$Session=="HDT55","Cap_size"],na.rm=TRUE)
# mean(data[data$Session=="CON","Cap_size"],na.rm=TRUE);sd(data[data$Session=="CON","Cap_size"],na.rm=TRUE)
# mean(data[data$Session=="LC","Cap_size"],na.rm=TRUE);sd(data[data$Session=="LC","Cap_size"],na.rm=TRUE)
# mean(data[data$Session=="ME","Cap_size"],na.rm=TRUE);sd(data[data$Session=="ME","Cap_size"],na.rm=TRUE)


remove<-data[is.na(data$Cap_size)!=FALSE ,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","Cap_size"]);qqline(test2[test2$Session=="BDC","Cap_size"])
# qqnorm(test2[test2$Session=="HDT55","Cap_size"]);qqline(test2[test2$Session=="HDT55","Cap_size"])
# qqnorm(test2[test2$Session=="CON","Cap_size"]);qqline(test2[test2$Session=="CON","Cap_size"])
# qqnorm(test2[test2$Session=="LC","Cap_size"]);qqline(test2[test2$Session=="LC","Cap_size"])
# qqnorm(test2[test2$Session=="ME","Cap_size"]);qqline(test2[test2$Session=="ME","Cap_size"])

a<-kruskal.test(Cap_size~Session, test2[test2$Group=="POST-VIRAL",])
all_pvals["ANOVA","Cap_size"]<-a$p.value

b<-
  t_test(test2[test2$Group=="BED REST",], Cap_size~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))



all_pvals["BDC-HDT55","Cap_size"]<-b$p
# all_pvals["CON-ME","Cap_size"]<-a$p.value[2]
# all_pvals["CON-LC","Cap_size"]<-a$p.value[1]
# all_pvals["LC-ME","Cap_size"]<-a$p.value[4]
all_pvals["BED_REST_test","Cap_size"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","Cap_size"]<-"kruskal_no_post_hoc"



# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55",paste("p<0.001"),"BED REST",
  "CON","ME",paste(round(a$p.value,3)),"POST-VIRAL")
# "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
# "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

Cap_size_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=Cap_size, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("Capillary Area ("*mu*m^2*")"))+
  scale_y_continuous(limits=c(0,170), breaks=seq(0,150,30),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(75), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))



ggsave(plot=Cap_size_a,
       filename = paste(output_folder,"/Cap_size_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

Cap_size_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=Cap_size, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("Capillary Area ("*mu*m^2*")"))+
  scale_y_continuous(limits=c(0,170), breaks=seq(0,150,30),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(160), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=Cap_size_b,
       filename = paste(output_folder,"/Cap_size_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")



# Supplemental 6 ----------------------------------------------------------

# Myoglobin ---------------------------------------------------------------


# Stats -------------------------------------------------------------------
# 
# mean(data[data$Session=="BDC","Myoglobin"],na.rm=TRUE); sd(data[data$Session=="BDC","Myoglobin"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","Myoglobin"],na.rm=TRUE);sd(data[data$Session=="HDT55","Myoglobin"],na.rm=TRUE)
# mean(data[data$Session=="CON","Myoglobin"],na.rm=TRUE);sd(data[data$Session=="CON","Myoglobin"],na.rm=TRUE)
# mean(data[data$Session=="LC","Myoglobin"],na.rm=TRUE);sd(data[data$Session=="LC","Myoglobin"],na.rm=TRUE)
# mean(data[data$Session=="ME","Myoglobin"],na.rm=TRUE);sd(data[data$Session=="ME","Myoglobin"],na.rm=TRUE)


remove<-data[is.na(data$Myoglobin)!=FALSE ,"Subject"]
test2<-test[!test$Subject %in% remove ,]

# qqnorm(test2[test2$Session=="BDC","Myoglobin"]);qqline(test2[test2$Session=="BDC","Myoglobin"])
# qqnorm(test2[test2$Session=="HDT55","Myoglobin"]);qqline(test2[test2$Session=="HDT55","Myoglobin"])
# qqnorm(test2[test2$Session=="CON","Myoglobin"]);qqline(test2[test2$Session=="CON","Myoglobin"])
# qqnorm(test2[test2$Session=="LC","Myoglobin"]);qqline(test2[test2$Session=="LC","Myoglobin"])
# qqnorm(test2[test2$Session=="ME","Myoglobin"]);qqline(test2[test2$Session=="ME","Myoglobin"])


a<-kruskal.test(Myoglobin~Session, test2[test2$Group=="POST-VIRAL",])
all_pvals["ANOVA","Myoglobin"]<-a$p.value


b<-
  t_test(test2[test2$Group=="BED REST",], Myoglobin~Session, paired=TRUE)

b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))



all_pvals["BDC-HDT55","Myoglobin"]<-b$p
# all_pvals["CON-ME","Myoglobin"]<-a$p.adj[2]
# all_pvals["CON-LC","Myoglobin"]<-a$p.adj[1]
# all_pvals["LC-ME","Myoglobin"]<-a$p.adj[3]
all_pvals["BED_REST_test","Myoglobin"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","Myoglobin"]<-"kruskal_no_post_hoc"


# Make graph --------------------------------------------------------------

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "BDC","HDT55",paste(format(round(b$p,3), drop0trailing=F)),"BED REST",
  "CON","ME",paste(round(a$p.value,3)),"POST-VIRAL")
# "CON","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"POST-VIRAL",
# "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-data[!data$Subject %in% remove,]

Myo_a<-
  ggplot(plot_data[plot_data$Group=="BED REST",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=Myoglobin, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="POST-VIRAL","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(bquote("Myoglobin Content (A.U.)"))+
  scale_y_continuous(limits=c(0,180), breaks=c(0,80,120,160),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(160), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(50),
              space=c(1),
              scales=c(50),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=Myo_a,
       filename = paste(output_folder,"/Myo_a_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

Myo_b<-
  ggplot(plot_data[plot_data$Group=="POST-VIRAL",], 
         aes(
           x=fct_relevel(Session, "BDC","HDT55","CON","LC","ME"),
           y=Myoglobin, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0,width=1.5/length(unique(plot_data[plot_data$Group=="BED REST","Session"])))+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(bquote("Myoglobin Content (A.U.)"))+
  scale_y_continuous(limits=c(0,180), breaks=c(0,80,120,160),expand=c(0,0))+
  scale_x_discrete(labels=c("BDC"="PRE", "HDT55"="POST","CON"="CON ","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(170), bracket.size=2,
                     label.size=14,tip.length = c(0.02,0.02),
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(50),
              space=c(1),
              scales=c(50),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=Myo_b,
       filename = paste(output_folder,"/Myo_b_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")


# Supplemental 7 ----------------------------------------------------------

# HDT55 v LC v ME ---------------------------------------------------------

selection<-c("K","N","G","P","Q1","U","T","L","C","E","F","D","H",
             "P110068","P110050","P110087","P110017","P110058","P110086","P110078","P110079","P110056","P110066","P110053","P110010","P110012",
             "P110118","P110138","P110117","P110141","P110128","P110129","P110144","P110142","P110122","P110120","P110134","P110143","P110124")

HLM_data<-data[data$Subject %in% selection & data$Session!="BDC",]

HLM_pvals<-as.data.frame(matrix(ncol=90, nrow=4))

colnames(HLM_pvals)<-colnames(HLM_data)
rownames(HLM_pvals)<-c("ANOVA","HDT55-ME","HDT55-LC","LC-ME")

median(HLM_data[HLM_data$Session=="HDT55","Age"],na.rm=TRUE);quantile(HLM_data[HLM_data$Session=="HDT55","Age"], na.rm=TRUE)
median(HLM_data[HLM_data$Session=="LC","Age"],na.rm=TRUE);quantile(HLM_data[HLM_data$Session=="LC","Age"], na.rm=TRUE)
median(HLM_data[HLM_data$Session=="ME","Age"],na.rm=TRUE);quantile(HLM_data[HLM_data$Session=="ME","Age"], na.rm=TRUE)

median(HLM_data[HLM_data$Session=="HDT55","Height"],na.rm=TRUE);quantile(HLM_data[HLM_data$Session=="HDT55","Height"], na.rm=TRUE)
median(HLM_data[HLM_data$Session=="LC","Height"],na.rm=TRUE);quantile(HLM_data[HLM_data$Session=="LC","Height"], na.rm=TRUE)
median(HLM_data[HLM_data$Session=="ME","Height"],na.rm=TRUE);quantile(HLM_data[HLM_data$Session=="ME","Height"], na.rm=TRUE)


median(HLM_data[HLM_data$Session=="HDT55","Weight"],na.rm=TRUE);quantile(HLM_data[HLM_data$Session=="HDT55","Weight"], na.rm=TRUE)
median(HLM_data[HLM_data$Session=="LC","Weight"],na.rm=TRUE);quantile(HLM_data[HLM_data$Session=="LC","Weight"], na.rm=TRUE)
median(HLM_data[HLM_data$Session=="ME","Weight"],na.rm=TRUE);quantile(HLM_data[HLM_data$Session=="ME","Weight"], na.rm=TRUE)


median(HLM_data[HLM_data$Session=="LC","Steps"],na.rm=TRUE);quantile(HLM_data[HLM_data$Session=="LC","Steps"], na.rm=TRUE)
median(HLM_data[HLM_data$Session=="ME","Steps"],na.rm=TRUE);quantile(HLM_data[HLM_data$Session=="ME","Steps"], na.rm=TRUE)


median(HLM_data[HLM_data$Session=="LC","Sx_duration"],na.rm=TRUE);quantile(HLM_data[HLM_data$Session=="LC","Sx_duration"], na.rm=TRUE)
median(HLM_data[HLM_data$Session=="ME","Sx_duration"],na.rm=TRUE);quantile(HLM_data[HLM_data$Session=="ME","Sx_duration"], na.rm=TRUE)


length(HLM_data[HLM_data$Session=="HDT55" & HLM_data$Sex=="Female","Age"])#/length(HLM_data[HLM_data$Session=="HDT55","Age"])
length(HLM_data[HLM_data$Session=="LC" & HLM_data$Sex=="Female","Age"])#/length(HLM_data[HLM_data$Session=="LC","Age"])
length(HLM_data[HLM_data$Session=="ME" & HLM_data$Sex=="Female","Age"])#/length(HLM_data[HLM_data$Session=="ME","Age"])

HLM_test<-as.data.frame(box_cox_transform(HLM_data))
HLM_test[,c("Subject","Session","Group", "Sex")]<-HLM_data[,c("Subject","Session","Group","Sex")]

remove<-HLM_test[is.na(HLM_test$Age)==TRUE,"Subject"]
test2<-HLM_test[!HLM_test$Subject %in% remove, ]

# qqnorm(test2[test2$Session=="ME","Age"]);qqline(test2[test2$Session=="ME","Age"]) 
# qqnorm(test2[test2$Session=="ME","Height"]);qqline(test2[test2$Session=="ME","Height"]) 
# qqnorm(test2[test2$Session=="ME","Weight"]);qqline(test2[test2$Session=="ME","Weight"]) 
# qqnorm(test2[test2$Session=="ME","Sx_duration"]);qqline(test2[test2$Session=="ME","Sx_duration"]) 
# qqnorm(test2[test2$Session=="ME","Steps"]);qqline(test2[test2$Session=="ME","Steps"]) 

kruskal.test(Age~Session, test2)
kruskal.test(Weight~Session, test2)
kruskal.test(Height~Session, test2)
wilcox.test(Sx_duration~Session, test2[test2$Session!="HDT55",])
wilcox.test(Steps~Session, test2[test2$Session!="HDT55",])



# HLM VO2 -----------------------------------------------------------------
remove<-HLM_data[is.na(HLM_data$VO2_abs)==TRUE   ,"Subject"]
HLM_test2<-HLM_test[!HLM_test$Subject %in% remove,]
# shapiro.test(HLM_test2[HLM_test2$Session=="HDT55","VO2_abs"])
# shapiro.test(HLM_test2[HLM_test2$Session=="LC","VO2_abs"])
# shapiro.test(HLM_test2[HLM_test2$Session=="ME","VO2_abs"])
# qqnorm(HLM_test2[HLM_test2$Session=="HDT55","VO2_abs"]);qqline(HLM_test2[HLM_test2$Session=="HDT55","VO2_rel"])
# qqnorm(HLM_test2[HLM_test2$Session=="LC","VO2_abs"]);qqline(HLM_test2[HLM_test2$Session=="LC","VO2_rel"])
# qqnorm(HLM_test2[HLM_test2$Session=="ME","VO2_abs"]);qqline(HLM_test2[HLM_test2$Session=="ME","VO2_rel"])

a<-
  kruskal.test(VO2_abs~Session, data=HLM_test)

HLM_pvals["ANOVA","VO2_abs"]<-a$p.value

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  # "HDT55","LC",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"BED REST",
  "HDT55","ME",paste(format(round(a$p.value[1],3), drop0trailing=F)),"POST-VIRAL")
# "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-HLM_data[!HLM_data$Subject %in% remove,]

VO2_HLM<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"HDT55","LC","ME"),y=VO2_abs, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("V\U0307" ~O[2]*"max (L"~min^-1*")"))+
  scale_y_continuous(limits=c(0,4.5), breaks=seq(0,4,1),expand=c(0,0))+
  scale_x_discrete(labels=c( "HDT55"="POST","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(4), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=VO2_HLM,
       filename = paste(output_folder,"/VO2_HLM_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")

# HLM GET -----------------------------------------------------------------
remove<-HLM_data[is.na(HLM_data$GET)==TRUE   ,"Subject"]
HLM_test2<-HLM_test[!HLM_test$Subject %in% remove,]
# shapiro.test(HLM_test2[HLM_test2$Session=="HDT55","GET"])
# shapiro.test(HLM_test2[HLM_test2$Session=="LC","GET"])  
# shapiro.test(HLM_test2[HLM_test2$Session=="ME","GET"])
# qqnorm(HLM_test2[HLM_test2$Session=="HDT55","GET"]);qqline(HLM_test2[HLM_test2$Session=="HDT55","GET"])
# qqnorm(HLM_test2[HLM_test2$Session=="LC","GET"]);qqline(HLM_test2[HLM_test2$Session=="LC","GET"])
# qqnorm(HLM_test2[HLM_test2$Session=="ME","GET"]);qqline(HLM_test2[HLM_test2$Session=="ME","GET"])

model<-lme(GET~Session, data=HLM_test2, random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
HLM_pvals["ANOVA","GET"]<-a$`Pr(>Chisq)`

a<-
  HLM_test2 %>% tukey_hsd(GET~Session)

HLM_pvals["HDT55-ME","GET"]<-a$p.adj[2]
HLM_pvals["HDT55-LC","GET"]<-a$p.adj[1]
HLM_pvals["LC-ME","GET"]<-a$p.adj[3]


stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "HDT55","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"BED REST",
  "HDT55","ME",paste(format(round(a$p.adj[2],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-HLM_data[!HLM_data$Subject %in% remove,]

GET_HLM<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"HDT55","LC","ME"),y=GET, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("GET (L"~min^-1*")"))+
  scale_y_continuous(limits=c(0,3.8), breaks=seq(0,3,1),expand=c(0,0))+
  scale_x_discrete(labels=c( "HDT55"="POST","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(3.6), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(3.2), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(2.8), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=GET_HLM,
       filename = paste(output_folder,"/GET_HLM_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")



# HLM GET % ---------------------------------------------------------------
remove<-HLM_data[is.na(HLM_data$GET_perc)==TRUE & HLM_data$GET_perc>0   ,"Subject"]
HLM_test2<-HLM_test[!HLM_test$Subject %in% remove,]
# shapiro.test(HLM_test2[HLM_test2$Session=="HDT55","GET_perc"])
# shapiro.test(HLM_test2[HLM_test2$Session=="LC","GET_perc"]) 
# shapiro.test(HLM_test2[HLM_test2$Session=="ME","GET_perc"])
# qqnorm(HLM_test2[HLM_test2$Session=="HDT55","GET_perc"]);qqline(HLM_test2[HLM_test2$Session=="HDT55","GET_perc"])
# qqnorm(HLM_test2[HLM_test2$Session=="LC","GET_perc"]);qqline(HLM_test2[HLM_test2$Session=="LC","GET_perc"])
# qqnorm(HLM_test2[HLM_test2$Session=="ME","GET_perc"]);qqline(HLM_test2[HLM_test2$Session=="ME","GET_perc"])


a<-
  kruskal.test(GET_perc~Session, data=HLM_test2)
HLM_pvals["ANOVA","GET_perc"]<-a$p.value

a<-
  pairwise.wilcox.test(HLM_test2[,"GET_perc"], HLM_test2[,"Session"], p.adjust.method="BH")

HLM_pvals["HDT55-ME","GET_perc"]<-a$p.value[2]
HLM_pvals["HDT55-LC","GET_perc"]<-a$p.value[1]
HLM_pvals["LC-ME","GET_perc"]<-a$p.value[4]


stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "HDT55","LC",paste(format(round(a$p.value[1],3), drop0trailing=F)),"BED REST",
  "HDT55","ME",paste("p<0.001"),"POST-VIRAL",
  # "HDT55","ME",paste(format(round(a$p.value[2],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.value[4],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-HLM_data[!HLM_data$Subject %in% remove,]

GET_perc_HLM<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"HDT55","LC","ME"),y=GET_perc*100, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("GET (% of V\U0307"~O[2]*"max)"))+
  scale_y_continuous(limits=c(0,130), breaks=c(0,40,60,80,100),expand=c(0,0))+
  scale_x_discrete(labels=c( "HDT55"="POST","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(120), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(110), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(100), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  scale_y_cut(breaks=c(30),
              space=c(0.5),
              scales=c(20),
              which=c(1),
              expand = expansion(mult = c(0.02,0.05)) )+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=GET_perc_HLM,
       filename = paste(output_folder,"/GET_perc_HLM_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")


# HLM SDH -----------------------------------------------------------------
remove<-HLM_data[is.na(HLM_data$SDH)==TRUE   ,"Subject"]
HLM_test2<-HLM_test[!HLM_test$Subject %in% remove,]
# shapiro.test(HLM_test2[HLM_test2$Session=="HDT55","SDH"])
# shapiro.test(HLM_test2[HLM_test2$Session=="LC","SDH"])
# shapiro.test(HLM_test2[HLM_test2$Session=="ME","SDH"])
# qqnorm(HLM_test2[HLM_test2$Session=="HDT55","SDH"]);qqline(HLM_test2[HLM_test2$Session=="HDT55","SDH"])
# qqnorm(HLM_test2[HLM_test2$Session=="LC","SDH"]);qqline(HLM_test2[HLM_test2$Session=="LC","SDH"])
# qqnorm(HLM_test2[HLM_test2$Session=="ME","SDH"]);qqline(HLM_test2[HLM_test2$Session=="ME","SDH"])

model<-lme(SDH~Session, data=HLM_test2, random = ~ 1|Subject, na.action = na.omit, control="optim")

a<-car::Anova(model)
HLM_pvals["ANOVA","SDH"]<-a$`Pr(>Chisq)`
a<-
  HLM_test2 %>% tukey_hsd(SDH~Session)
HLM_pvals["HDT55-ME","SDH"]<-a$p.adj[2]
HLM_pvals["HDT55-LC","SDH"]<-a$p.adj[1]
HLM_pvals["LC-ME","SDH"]<-a$p.adj[3]


stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "HDT55","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"BED REST",
  "HDT55","ME",paste(format(round(a$p.adj[2],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-HLM_data[!HLM_data$Subject %in% remove,]

SDH_HLM<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"HDT55","LC","ME"),y=SDH*10^5, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression(atop("SDH Activity", "(\U0394"~A[660]*~mu*m^-1*""*~s^-1*~10^-5*")")))+
  scale_y_continuous(limits=c(0,2.2), breaks=seq(0,1.5,0.5),expand=c(0,0))+
  scale_x_discrete(labels=c( "HDT55"="POST","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(2.0), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(1.8), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(1.6), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=SDH_HLM,
       filename = paste(output_folder,"/SDH_HLM_",date,".png", sep = ""),
       device="png",  width = 10, height = 14, units = "in")



# HLM Oxphos --------------------------------------------------------------
remove<-HLM_data[is.na(HLM_data$Oxphos)==TRUE | HLM_data$membrane_intact>1.1  ,"Subject"]
HLM_test2<-HLM_test[!HLM_test$Subject %in% remove,]
# shapiro.test(HLM_test2[HLM_test2$Session=="HDT55","Oxphos"])
# shapiro.test(HLM_test2[HLM_test2$Session=="LC","Oxphos"]) 
# shapiro.test(HLM_test2[HLM_test2$Session=="ME","Oxphos"])
# qqnorm(HLM_test2[HLM_test2$Session=="HDT55","Oxphos"]);qqline(HLM_test2[HLM_test2$Session=="HDT55","Oxphos"])
# qqnorm(HLM_test2[HLM_test2$Session=="LC","Oxphos"]);qqline(HLM_test2[HLM_test2$Session=="LC","Oxphos"])
# qqnorm(HLM_test2[HLM_test2$Session=="ME","Oxphos"]);qqline(HLM_test2[HLM_test2$Session=="ME","Oxphos"])

a<-kruskal.test(Oxphos~Session, data=HLM_test2)
HLM_pvals["ANOVA","Oxphos"]<-a$p.value

a<-
  pairwise.wilcox.test(HLM_test2[,"Oxphos"], HLM_test2[,"Session"], p.adjust.method="BH")

HLM_pvals["HDT55-ME","Oxphos"]<-a$p.value[2]
HLM_pvals["HDT55-LC","Oxphos"]<-a$p.value[1]
HLM_pvals["LC-ME","Oxphos"]<-a$p.value[4]


stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "HDT55","LC",paste(format(round(a$p.value[1],3), drop0trailing=F)),"BED REST",
  "HDT55","ME",paste(format(round(a$p.value[2],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.value[4],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-HLM_data[!HLM_data$Subject %in% remove,]

Oxphos_HLM<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"HDT55","LC","ME"),y=Oxphos, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(bquote(atop("Oxidative Phosphorylation","(pmol"*~s^-1*~mg^-1*")")))+
  scale_y_continuous(limits=c(0,180), breaks=seq(0,150,50),expand=c(0,0))+
  scale_x_discrete(labels=c( "HDT55"="POST","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(150), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(135), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(120), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=Oxphos_HLM,
       filename = paste(output_folder,"/Oxphos_HLM_",date,".png", sep = ""),
       device="png",  width = 10, height = 14, units = "in")



# HLM Cap : Fibre ---------------------------------------------------------
remove<-HLM_data[is.na(HLM_data$CF)==TRUE   ,"Subject"]
HLM_test2<-HLM_test[!HLM_test$Subject %in% remove,]
# shapiro.test(HLM_test2[HLM_test2$Session=="HDT55","CF"])
# shapiro.test(HLM_test2[HLM_test2$Session=="LC","CF"])
# shapiro.test(HLM_test2[HLM_test2$Session=="ME","CF"])
# qqnorm(HLM_test2[HLM_test2$Session=="HDT55","CF"]);qqline(HLM_test2[HLM_test2$Session=="HDT55","CF"])
# qqnorm(HLM_test2[HLM_test2$Session=="LC","CF"]);qqline(HLM_test2[HLM_test2$Session=="LC","CF"])
# qqnorm(HLM_test2[HLM_test2$Session=="ME","CF"]);qqline(HLM_test2[HLM_test2$Session=="ME","CF"])


model<-lme(CF~Session, data=HLM_test2, random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
HLM_pvals["ANOVA","CF"]<-a$`Pr(>Chisq)`

a<-
  HLM_test2 %>% tukey_hsd(CF~Session)

HLM_pvals["HDT55-ME","CF"]<-a$p.adj[2]
HLM_pvals["HDT55-LC","CF"]<-a$p.adj[1]
HLM_pvals["LC-ME","CF"]<-a$p.adj[3]



stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "HDT55","LC",paste(format(round(a$p.adj[1],3), drop0trailing=F)),"BED REST",
  "HDT55","ME",paste(format(round(a$p.adj[2],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.adj[3],3), drop0trailing=F)),"POST-VIRAL")


plot_data<-HLM_data[!HLM_data$Subject %in% remove,]

CF_HLM<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"HDT55","LC","ME"),y=CF, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("Capillary:Fibre Ratio"))+
  scale_y_continuous(limits=c(0,5), breaks=seq(0,4,1),expand=c(0,0))+
  scale_x_discrete(labels=c( "HDT55"="POST","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(4.0), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(3.6), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(3.2), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=CF_HLM,
       filename = paste(output_folder,"/CF_HLM_",date,".png", sep = ""),
       device="png",  width = 9, height = 14, units = "in")




# HLM FCSA ----------------------------------------------------------------
remove<-HLM_data[is.na(HLM_data$FCSA)==TRUE   ,"Subject"]
HLM_test2<-HLM_test[!HLM_test$Subject %in% remove,]
# shapiro.test(HLM_test2[HLM_test2$Session=="HDT55","FCSA"]) 
# shapiro.test(HLM_test2[HLM_test2$Session=="LC","FCSA"]) 
# shapiro.test(HLM_test2[HLM_test2$Session=="ME","FCSA"])
# qqnorm(HLM_test2[HLM_test2$Session=="HDT55","FCSA"]);qqline(HLM_test2[HLM_test2$Session=="HDT55","FCSA"])
# qqnorm(HLM_test2[HLM_test2$Session=="LC","FCSA"]);qqline(HLM_test2[HLM_test2$Session=="LC","FCSA"])
# qqnorm(HLM_test2[HLM_test2$Session=="ME","FCSA"]);qqline(HLM_test2[HLM_test2$Session=="ME","FCSA"])

a<-kruskal.test(FCSA~Session, data=HLM_test2)
HLM_pvals["ANOVA","FCSA"]<-a$p.value

a<-
  pairwise.wilcox.test(HLM_test2[,"FCSA"], HLM_test2[,"Session"], p.adjust.method="BH")

HLM_pvals["HDT55-ME","FCSA"]<-a$p.value[2]
HLM_pvals["HDT55-LC","FCSA"]<-a$p.value[1]
HLM_pvals["LC-ME","FCSA"]<-a$p.value[4]

stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj, ~Group,
  "HDT55","LC",paste(format(round(a$p.value[1],3), drop0trailing=F)),"BED REST",
  "HDT55","ME",paste(format(round(a$p.value[2],3), drop0trailing=F)),"POST-VIRAL",
  "LC","ME",paste(format(round(a$p.value[4],3), nsmall=3)),"POST-VIRAL")


plot_data<-HLM_data[!HLM_data$Subject %in% remove,]

FCSA_HLM<-
  ggplot(plot_data, 
         aes(x=fct_relevel(Session,"HDT55","LC","ME"),y=FCSA, group=Session,fill=Session))+
  geom_boxplot(linewidth=2, outlier.shape = NA, coef=0)+
  geom_point(size=6,stroke=2 , shape=21,position = position_jitter(width=0.2, height=0, seed=1), fill="white", colour="black")+
  scale_fill_manual(values = c("HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]))+
  ylab(expression("FCSA ("*mu*m^2*")"))+
  scale_y_continuous(limits=c(0,9000), breaks=seq(0,8000,2000),expand=c(0,0))+
  scale_x_discrete(labels=c( "HDT55"="POST","LC"="LC","ME"="ME"))+
  stat_pvalue_manual(data=stat_test[2,], label = "p.adj",y.position=c(8600), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[1,], label = "p.adj",y.position=c(8000), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  stat_pvalue_manual(data=stat_test[3,], label = "p.adj",y.position=c(7400), bracket.size=2,
                     label.size=14,tip.length = 0.02,
                     linetype="solid", inherit.aes=FALSE)+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
        axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 30, b = 0, l = 0)), # for atop functions export as 9.5x14in
        axis.title.x = element_blank(),
        axis.text.y = element_text( hjust=1,size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( vjust=0,size =rel((base_size+title_text_rel_size)/base_size))
  )+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=FCSA_HLM,
       filename = paste(output_folder,"/FCSA_HLM_",date,".png", sep = ""),
       device="png",  width = 10, height = 14, units = "in")


# HLM VO2 vs SDH ---------------------------------------------------------
remove<-HLM_data[is.na(HLM_data$SDH)==TRUE|is.na(HLM_data$VO2_rel)==TRUE   ,"Subject"]
plot_data<-HLM_data[!HLM_data$Subject %in% remove,]

df<-list(plot_data[plot_data$Session=="LC",],plot_data[plot_data$Session=="HDT55",])
cocor(~VO2_rel+SDH| VO2_rel+SDH, df )
df<-list(plot_data[plot_data$Session=="ME",],plot_data[plot_data$Session=="HDT55",])
cocor(~VO2_rel+SDH| VO2_rel+SDH, df )
df<-list(plot_data[plot_data$Session=="LC",],plot_data[plot_data$Session=="ME",])
cocor(~VO2_rel+SDH| VO2_rel+SDH, df )


x<-plot_data$SDH
y<-plot_data$VO2_rel
Session<-plot_data$Session
HDT55_cor<-cor.test(x[Session == "HDT55"], y[Session == "HDT55"], method="pearson")
LC_cor<- cor.test(x[Session == "LC"], y[Session == "LC"], method="pearson")
ME_cor<- cor.test(x[Session == "ME"], y[Session == "ME"], method="pearson")
cor_data_1<-rbind(HDT55_cor$p.value,LC_cor$p.value,ME_cor$p.value)
cor_data_2<-rbind(HDT55_cor$estimate,LC_cor$estimate,ME_cor$estimate)
cor_data_3<-c("HDT55","LC","ME")
cor_data_4<-c("BED REST","POST-VIRAL","POST-VIRAL")
cor_data<-as.data.frame(cbind(cor_data_1,cor_data_2, cor_data_3, cor_data_4))
colnames(cor_data)<-c("p_value","r","Session","Group")
cor_data$p_value<-as.numeric(cor_data$p_value)
cor_data$p_value<-signif(cor_data$p_value,2)
cor_data$r<-as.numeric(cor_data$r)
cor_data$r<-round(cor_data$r,2)
cor_data<-add_significance(cor_data, p.col = "p_value" ,cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
cor_data$comparison<-"SDH_v_VO2_HLM"

cor_pvals<-rbind(cor_pvals,cor_data)

lab_set<-data.frame(x1=c(0.9,0.9,0.9),y1=c(67,60,53),Group=c("BED REST","POST-VIRAL","POST-VIRAL"),
                    labs=c(paste("POST-Bed Rest: r=",format(round(cor_data[cor_data$Session=="HDT55","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="HDT55","p_value"],3),nsmall=3), sep=""),
                           paste("LC: r=",format(round(cor_data[cor_data$Session=="LC","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="LC","p_value"],3),nsmall=3), sep=""),
                           paste("ME: r=",format(round(cor_data[cor_data$Session=="ME","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="ME","p_value"],3),nsmall=3), sep="")))


VO2_v_SDH_HLM<-
  ggplot(data=plot_data, aes(y=VO2_rel, x=SDH*10^5))+
  geom_jitter(size=7.5,aes( shape=Session,fill=Session),width = 0, height=0, stroke=2)+
  scale_shape_manual(values=c("HDT55"=21,"LC"=24, "ME"=24),
                     labels=c("HDT55"="POST\nBed Rest","LC"="LC", "ME"="ME"))+
  scale_fill_manual(values = c("HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]),
                    labels=c("HDT55"="POST\nBed Rest","LC"="LC", "ME"="ME"))+
  scale_colour_manual(values = c("HDT55"="white","LC"="white","ME"="white"),
                      labels=c("HDT55"="POST\nBed Rest","LC"="LC", "ME"="ME"))+
  geom_smooth(data=plot_data[plot_data$Session=="HDT55",],method="lm",formula = y~x, se=FALSE, linetype= "solid",fill="white",colour="black", size=2, fullrange=F)+
  # geom_smooth(data=plot_data[plot_data$Session=="LC",],method="lm", se=FALSE, linetype= "solid",fill="white",colour=colours3[3], size=2)+
  # geom_smooth(data=plot_data[plot_data$Session=="ME",],method="lm", se=FALSE, linetype= "solid",fill="white",colour=colours3[4], size=2)+
  xlab(expression(atop("SDH Activity", "(\U0394"~A[660]*~mu*m^-1*""*~s^-1*~10^-5*")")))+
  ylab(expression("V\U0307" ~O[2][max]*"(mL"~min^-1*~kg^-1*")"))+
  scale_x_continuous(limits=c(0,2.2), breaks=seq(0.5,2.0,0.5), expand=c(0,0))+
  scale_y_continuous(limits = c(0,80), breaks=seq(0,60,20), expand=c(0,0))+
  theme(
    aspect.ratio = 1/1.2,
    legend.position = "none",
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.line=element_line(colour="black", size=line_size/2),
    axis.ticks = element_line(colour="black", size = line_size/2),
    text=element_text(size=base_size, family = "Arial"),
    panel.background  = element_rect(fill="white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      size = rel((title_text_rel_size + base_size) / base_size),
      hjust = 0.5
    ),
    axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
    axis.title.y = element_text(angle = 90, vjust = 1,size = rel((title_text_rel_size + base_size) / base_size),),
    axis.text.y = element_text( size = rel((title_text_rel_size + base_size) / base_size)
    ),
    axis.text.x=element_text(angle=0, size=rel((title_text_rel_size + base_size) / base_size)),
    strip.background = element_blank(),
    strip.text = element_blank())+
  guides(y=guide_axis(cap="upper"))+
  geom_text(data=lab_set,aes(y=y1,x=x1, group=Group,label=labs), family="Arial", fontface="bold", size=base_size/4, inherit.aes = TRUE )


ggsave(plot=VO2_v_SDH_HLM,
       filename = paste(output_folder,"/VO2_v_SDH_HLM_",date,".png", sep = ""),
       device="png",  width = 20, height = 16, units = "in")

# HLM VO2 vs Oxphos -------------------------------------------------------
remove<-HLM_data[is.na(HLM_data$Oxphos)==TRUE|is.na(HLM_data$VO2_rel)==TRUE &HLM_data$membrane_intact>1.1   ,"Subject"]
plot_data<-HLM_data[!HLM_data$Subject %in% remove,]

df<-list(plot_data[plot_data$Session=="LC",],plot_data[plot_data$Session=="HDT55",])
cocor(~Oxphos+VO2_rel| Oxphos+VO2_rel, df )
df<-list(plot_data[plot_data$Session=="ME",],plot_data[plot_data$Session=="HDT55",])
cocor(~Oxphos+VO2_rel| Oxphos+VO2_rel, df )
df<-list(plot_data[plot_data$Session=="LC",],plot_data[plot_data$Session=="ME",])
cocor(~Oxphos+VO2_rel| Oxphos+VO2_rel, df )


x<-plot_data$Oxphos
y<-plot_data$VO2_rel
Session<-plot_data$Session
HDT55_cor<-cor.test(x[Session == "HDT55"], y[Session == "HDT55"], method="pearson")
LC_cor<- cor.test(x[Session == "LC"], y[Session == "LC"], method="pearson")
ME_cor<- cor.test(x[Session == "ME"], y[Session == "ME"], method="pearson")
cor_data_1<-rbind(HDT55_cor$p.value,LC_cor$p.value,ME_cor$p.value)
cor_data_2<-rbind(HDT55_cor$estimate,LC_cor$estimate,ME_cor$estimate)
cor_data_3<-c("HDT55","LC","ME")
cor_data_4<-c("BED REST","POST-VIRAL","POST-VIRAL")
cor_data<-as.data.frame(cbind(cor_data_1,cor_data_2, cor_data_3, cor_data_4))
colnames(cor_data)<-c("p_value","r","Session","Group")
cor_data$p_value<-as.numeric(cor_data$p_value)
cor_data$p_value<-signif(cor_data$p_value,2)
cor_data$r<-as.numeric(cor_data$r)
cor_data$r<-round(cor_data$r,2)
cor_data<-add_significance(cor_data, p.col = "p_value" ,cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
cor_data$comparison<-"Oxphos_v_VO2_HLM"

cor_pvals<-rbind(cor_pvals,cor_data)

lab_set<-data.frame(x1=c(80,80,80),y1=c(65,60,55),Group=c("BED REST","POST-VIRAL","POST-VIRAL"),
                    labs=c(paste("POST-Bed Rest: r=",format(round(cor_data[cor_data$Session=="HDT55","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="HDT55","p_value"],3),nsmall=3)),
                           paste("LC: r=",format(round(cor_data[cor_data$Session=="LC","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="LC","p_value"],3),nsmall=3)),
                           paste("ME: r=",format(round(cor_data[cor_data$Session=="ME","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="ME","p_value"],3),nsmall=3))))


VO2_v_oxphos_HLM<-
  ggplot(data=plot_data, aes(y=VO2_rel, x=Oxphos))+
  geom_jitter(size=7.5,aes( shape=Session,fill=Session),width = 0, height=0, stroke=2)+
  scale_shape_manual(values=c("HDT55"=21,"LC"=24, "ME"=24),
                     labels=c("HDT55"="POST\nBed Rest","LC"="LC", "ME"="ME"))+
  scale_fill_manual(values = c("HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]),
                    labels=c("HDT55"="POST\nBed Rest","LC"="LC", "ME"="ME"))+
  scale_colour_manual(values = c("HDT55"="white","LC"="white","ME"="white"),
                      labels=c("HDT55"="POST\nBed Rest","LC"="LC", "ME"="ME"))+
  geom_smooth(data=plot_data[plot_data$Session=="HDT55",],method="lm", se=FALSE, linetype= "solid",fill="white",colour="black", size=2)+
  # geom_smooth(data=plot_data[plot_data$Session=="LC",],method="lm", se=FALSE, linetype= "solid",fill="white",colour=colours3[3], size=2)+
  # geom_smooth(data=plot_data[plot_data$Session=="ME",],method="lm", se=FALSE, linetype= "solid",fill="white",colour=colours3[4], size=2)+
  xlab(bquote(atop("Oxidative Phosphorylation","(pmol"*~s^-1*~mg^-1*")")))+
  ylab(expression("V\U0307" ~O[2][max]*"(mL"~min^-1*~kg^-1*")"))+
  scale_y_continuous(limits = c(0,75), breaks=seq(0,60,20), expand=c(0,0))+
  scale_x_continuous(limits=c(0,160), breaks=seq(30,150,30), expand=c(0,0))+
  theme(
    aspect.ratio = 1/1.2,
    legend.position = "none",
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.line=element_line(colour="black", size=line_size/2),
    axis.ticks = element_line(colour="black", size = line_size/2),
    text=element_text(size=base_size, family = "Arial"),
    panel.background  = element_rect(fill="white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      size = rel((title_text_rel_size + base_size) / base_size),
      hjust = 0.5
    ),
    axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
    axis.title.y = element_text(angle = 90, vjust = 1,size = rel((title_text_rel_size + base_size) / base_size),),
    axis.text.y = element_text( size = rel((title_text_rel_size + base_size) / base_size)
    ),
    axis.text.x=element_text(angle=0, size=rel((title_text_rel_size + base_size) / base_size)),
    strip.background = element_blank(),
    strip.text = element_blank())+
  guides(y=guide_axis(cap="upper"))+
  geom_text(data=lab_set,aes(y=y1,x=x1, group=Group,label=labs), family="Arial", fontface="bold", size=base_size/4, inherit.aes = TRUE )


ggsave(plot=VO2_v_oxphos_HLM,
       filename = paste(output_folder,"/VO2_v_oxphos_HLM_",date,".png", sep = ""),
       device="png",  width = 20, height = 16, units = "in")

# HLM CF vs FCSA ----------------------------------------------------------
remove<-HLM_data[is.na(HLM_data$CF)==TRUE|is.na(HLM_data$FCSA)==TRUE   ,"Subject"]
plot_data<-HLM_data[!HLM_data$Subject %in% remove,]

df<-list(plot_data[plot_data$Session=="LC",],plot_data[plot_data$Session=="HDT55",])
cocor(~FCSA+CF| FCSA  + CF, df )
df<-list(plot_data[plot_data$Session=="ME",],plot_data[plot_data$Session=="HDT55",])
cocor(~FCSA+CF| FCSA  + CF, df )
df<-list(plot_data[plot_data$Session=="LC",],plot_data[plot_data$Session=="ME",])
cocor(~FCSA+CF| FCSA  + CF, df )


x<-plot_data$CF
y<-plot_data$FCSA
Session<-plot_data$Session
HDT55_cor<-cor.test(x[Session == "HDT55"], y[Session == "HDT55"], method="pearson")
LC_cor<- cor.test(x[Session == "LC"], y[Session == "LC"], method="pearson")
ME_cor<- cor.test(x[Session == "ME"], y[Session == "ME"], method="pearson")
cor_data_1<-rbind(HDT55_cor$p.value,LC_cor$p.value,ME_cor$p.value)
cor_data_2<-rbind(HDT55_cor$estimate,LC_cor$estimate,ME_cor$estimate)
cor_data_3<-c("HDT55","LC","ME")
cor_data_4<-c("BED REST","POST-VIRAL","POST-VIRAL")
cor_data<-as.data.frame(cbind(cor_data_1,cor_data_2, cor_data_3, cor_data_4))
colnames(cor_data)<-c("p_value","r","Session","Group")
cor_data$p_value<-as.numeric(cor_data$p_value)
cor_data$p_value<-signif(cor_data$p_value,2)
cor_data$r<-as.numeric(cor_data$r)
cor_data$r<-round(cor_data$r,2)
cor_data<-add_significance(cor_data, p.col = "p_value" ,cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
cor_data$comparison<-"CF_v_FCSA_HLM"

cor_pvals<-rbind(cor_pvals,cor_data)

lab_set<-data.frame(x1=c(3400,3400,3400),y1=c(3.6,3.3,3.0),Group=c("BED REST","POST-VIRAL","POST-VIRAL"),
                    labs=c(paste("POST-Bed Rest: r=",format(round(cor_data[cor_data$Session=="HDT55","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="HDT55","p_value"],3),nsmall=3)),
                           paste("LC: r=",format(round(cor_data[cor_data$Session=="LC","r"],3),nsmall=2),", p<0.001"),
                           paste("ME: r=",format(round(cor_data[cor_data$Session=="ME","r"],3),nsmall=2),", p=",format(round(cor_data[cor_data$Session=="ME","p_value"],3),nsmall=3))))


CF_v_FCSA_HLM<-
  ggplot(data=plot_data, aes(y=CF, x=FCSA))+
  geom_jitter(size=7.5,aes( shape=Session,fill=Session),width = 0, height=0, stroke=2)+
  scale_shape_manual(values=c("HDT55"=21,"LC"=24, "ME"=24),
                     labels=c("HDT55"="POST\nBed Rest","LC"="LC", "ME"="ME"))+
  scale_fill_manual(values = c("HDT55"=colours3[2], "LC"=colours3[3], "ME"=colours3[4]),
                    labels=c("HDT55"="POST\nBed Rest","LC"="LC", "ME"="ME"))+
  scale_colour_manual(values = c("HDT55"="white","LC"="white","ME"="white"),
                      labels=c("HDT55"="POST\nBed Rest","LC"="LC", "ME"="ME"))+
  # geom_smooth(data=plot_data[plot_data$Session=="HDT55",],method="lm", se=FALSE, linetype= "solid",fill="white",colour="black", size=2)+
  geom_smooth(data=plot_data[plot_data$Session=="LC",],method="lm",formula = y~x, se=FALSE, linetype= "solid",fill="white",colour=colours3[3], size=2, fullrange=F)+
  geom_smooth(data=plot_data[plot_data$Session=="ME",],method="lm",formula = y~x, se=FALSE, linetype= "dashed",fill="white",colour=colours3[4], size=2, fullrange=F)+
  xlab(expression("FCSA("*mu*m^2*")"))+
  ylab(expression("Capillary : Fibre Ratio"))+
  scale_x_continuous(limits=c(0,7500), breaks=c(2500,5000,7500), expand=c(0,0))+
  scale_y_continuous(limits = c(0,4.5), breaks=seq(0,4,1), expand=c(0,0))+
  theme(
    aspect.ratio = 1/1.2,
    legend.position = "none",
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.line=element_line(colour="black", size=line_size/2),
    axis.ticks = element_line(colour="black", size = line_size/2),
    text=element_text(size=base_size, family = "Arial"),
    panel.background  = element_rect(fill="white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      size = rel((title_text_rel_size + base_size) / base_size),
      hjust = 0.5
    ),
    axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)),
    axis.title.y = element_text(angle = 90, vjust = 1,size = rel((title_text_rel_size + base_size) / base_size),),
    axis.text.y = element_text( size = rel((title_text_rel_size + base_size) / base_size)
    ),
    axis.text.x=element_text(angle=0, size=rel((title_text_rel_size + base_size) / base_size)),
    strip.background = element_blank(),
    strip.text = element_blank())+
  guides(y=guide_axis(cap="upper"))+
  geom_text(data=lab_set,aes(y=y1,x=x1, group=Group,label=labs), family="Arial", fontface="bold", size=base_size/4, inherit.aes = TRUE )


ggsave(plot=CF_v_FCSA_HLM,
       filename = paste(output_folder,"/CF_v_FCSA_HLM_",date,".png", sep = ""),
       device="png",  width = 20, height = 16, units = "in")


# Figure 2C ---------------------------------------------------------------
# Fibre type proportions ---------------------------------------------------

# Stats -------------------------------------------------------------------

data<-as.data.frame(read_xlsx("Manuscript_data_clean_230525.xlsx", sheet="Sheet2"))

data<-data %>%
  mutate(Group=dplyr::recode(Group, "AGBRESA"="BED REST","MUSCLE-ME"="POST-VIRAL"))

data$Percent_IIa_IIx_IIx<-data$Percent_IIa_IIx+data$Percent_IIx

# Check normality ---------------------------------------------------------

op<-par(mfrow=c(1,5), mar=c(2,2,2,2))
sessions<-unique(data$Session)
for (i in 4:length(colnames(data)) ){
  try({
    parameter<-colnames(data[i])
    for (j in 1:length(sessions)){
      session<-sessions[j]
      tmp<- data[data$Session==session,parameter]
      qqnorm(tmp,xlab=parameter, main=paste(session,parameter,sep="_"))
      qqline(tmp)
      
    }
    
  },silent=TRUE)
  
  
}

norm_p_vals<-normality_test(data)
# Fibre Percents
test<-as.data.frame(box_cox_transform(data))
test<-test %>% 
  filter_all(all_vars(!is.infinite(.)))
test[,c("Subject","Session","Group","Sex")]<-data[,c("Subject","Session","Group","Sex")]
test_norm<-normality_test(test)
# remove<-data[is.na(data$AHRR)==TRUE,"Subject"]
# test2<-test[!test$Subject %in% remove,]
# 
# mean(data[data$Session=="BDC","Percent_I"],na.rm=TRUE); sd(data[data$Session=="BDC","Percent_I"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","Percent_I"],na.rm=TRUE); sd(data[data$Session=="HDT55","Percent_I"],na.rm=TRUE)
# mean(data[data$Session=="CON","Percent_I"],na.rm=TRUE); sd(data[data$Session=="CON","Percent_I"],na.rm=TRUE)
# mean(data[data$Session=="LC","Percent_I"],na.rm=TRUE); sd(data[data$Session=="LC","Percent_I"],na.rm=TRUE)
# mean(data[data$Session=="ME","Percent_I"],na.rm=TRUE); sd(data[data$Session=="ME","Percent_I"],na.rm=TRUE)
# 
# mean(data[data$Session=="BDC","Percent_I_IIa"],na.rm=TRUE); sd(data[data$Session=="BDC","Percent_I_IIa"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","Percent_I_IIa"],na.rm=TRUE); sd(data[data$Session=="HDT55","Percent_I_IIa"],na.rm=TRUE)
# mean(data[data$Session=="CON","Percent_I_IIa"],na.rm=TRUE); sd(data[data$Session=="CON","Percent_I_IIa"],na.rm=TRUE)
# mean(data[data$Session=="LC","Percent_I_IIa"],na.rm=TRUE); sd(data[data$Session=="LC","Percent_I_IIa"],na.rm=TRUE)
# mean(data[data$Session=="ME","Percent_I_IIa"],na.rm=TRUE); sd(data[data$Session=="ME","Percent_I_IIa"],na.rm=TRUE)
# 
# mean(data[data$Session=="BDC","Percent_IIa"],na.rm=TRUE); sd(data[data$Session=="BDC","Percent_IIa"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","Percent_IIa"],na.rm=TRUE); sd(data[data$Session=="HDT55","Percent_IIa"],na.rm=TRUE)
# mean(data[data$Session=="CON","Percent_IIa"],na.rm=TRUE); sd(data[data$Session=="CON","Percent_IIa"],na.rm=TRUE)
# mean(data[data$Session=="LC","Percent_IIa"],na.rm=TRUE); sd(data[data$Session=="LC","Percent_IIa"],na.rm=TRUE)
# mean(data[data$Session=="ME","Percent_IIa"],na.rm=TRUE); sd(data[data$Session=="ME","Percent_IIa"],na.rm=TRUE)
# 
# mean(data[data$Session=="BDC","Percent_IIa_IIx"],na.rm=TRUE); sd(data[data$Session=="BDC","Percent_IIa_IIx"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","Percent_IIa_IIx"],na.rm=TRUE); sd(data[data$Session=="HDT55","Percent_IIa_IIx"],na.rm=TRUE)
# mean(data[data$Session=="CON","Percent_IIa_IIx"],na.rm=TRUE); sd(data[data$Session=="CON","Percent_IIa_IIx"],na.rm=TRUE)
# mean(data[data$Session=="LC","Percent_IIa_IIx"],na.rm=TRUE); sd(data[data$Session=="LC","Percent_IIa_IIx"],na.rm=TRUE)
# mean(data[data$Session=="ME","Percent_IIa_IIx"],na.rm=TRUE); sd(data[data$Session=="ME","Percent_IIa_IIx"],na.rm=TRUE)
# 
# mean(data[data$Session=="BDC","Percent_IIx"],na.rm=TRUE); sd(data[data$Session=="BDC","Percent_IIx"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","Percent_IIx"],na.rm=TRUE); sd(data[data$Session=="HDT55","Percent_IIx"],na.rm=TRUE)
# mean(data[data$Session=="CON","Percent_IIx"],na.rm=TRUE); sd(data[data$Session=="CON","Percent_IIx"],na.rm=TRUE)
# mean(data[data$Session=="LC","Percent_IIx"],na.rm=TRUE); sd(data[data$Session=="LC","Percent_IIx"],na.rm=TRUE)
# mean(data[data$Session=="ME","Percent_IIx"],na.rm=TRUE); sd(data[data$Session=="ME","Percent_IIx"],na.rm=TRUE)

# Type I  
# qqnorm(test[test$Session=="BDC","Percent_I"]);qqline(test[test$Session=="BDC","Percent_I"])
# qqnorm(test[test$Session=="HDT55","Percent_I"]);qqline(test[test$Session=="HDT55","Percent_I"])
# qqnorm(test[test$Session=="CON","Percent_I"]);qqline(test[test$Session=="CON","Percent_I"])
# qqnorm(test[test$Session=="LC","Percent_I"]);qqline(test[test$Session=="LC","Percent_I"])
# qqnorm(test[test$Session=="ME","Percent_I"]);qqline(test[test$Session=="ME","Percent_I"])

model <-lme(Percent_I~ Session, data=test[test$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")

a<-car::Anova(model)
all_pvals["ANOVA","Percent_I"]<-a$`Pr(>Chisq)`

a<-
  test[test$Group=="POST-VIRAL",] %>% tukey_hsd(Percent_I~Session)   

b<-
  t_test(test[test$Group=="BED REST",], Percent_I~Session, paired=TRUE)
b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","Percent_I"]<-b$p
all_pvals["CON-ME","Percent_I"]<-a$p.adj[2]
all_pvals["CON-LC","Percent_I"]<-a$p.adj[1]
all_pvals["LC-ME","Percent_I"]<-a$p.adj[3]
all_pvals["BED_REST_test","Percent_I"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","Percent_I"]<-"anova_tukey_post_hoc"


# Type I/IIa  
# qqnorm(test[test$Session=="BDC","Percent_I_IIa"]);qqline(test[test$Session=="BDC","Percent_I_IIa"])
# qqnorm(test[test$Session=="HDT55","Percent_I_IIa"]);qqline(test[test$Session=="HDT55","Percent_I_IIa"])
# qqnorm(test[test$Session=="CON","Percent_I_IIa"]);qqline(test[test$Session=="CON","Percent_I_IIa"])
# qqnorm(test[test$Session=="LC","Percent_I_IIa"]);qqline(test[test$Session=="LC","Percent_I_IIa"])
# qqnorm(test[test$Session=="ME","Percent_I_IIa"]);qqline(test[test$Session=="ME","Percent_I_IIa"])

model <-lme(Percent_I_IIa~ Session, data=test[test$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")

a<-car::Anova(model)
all_pvals["ANOVA","Percent_I_IIa"]<-a$`Pr(>Chisq)`
b<-
  t_test(test[test$Group=="BED REST",], Percent_I_IIa~Session, paired=TRUE)
b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","Percent_I_IIa"]<-b$p
# all_pvals["CON-ME","Percent_I_IIa"]<-a$p.adj[2]
# all_pvals["CON-LC","Percent_I_IIa"]<-a$p.adj[1]
# all_pvals["LC-ME","Percent_I_IIa"]<-a$p.adj[3]
all_pvals["BED_REST_test","Percent_I_IIa"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","Percent_I_IIa"]<-"anova_no_post_hoc"



# Type IIa

# qqnorm(test[test$Session=="BDC","Percent_IIa"]);qqline(test[test$Session=="BDC","Percent_IIa"])
# qqnorm(test[test$Session=="HDT55","Percent_IIa"]);qqline(test[test$Session=="HDT55","Percent_IIa"])
# qqnorm(test[test$Session=="CON","Percent_IIa"]);qqline(test[test$Session=="CON","Percent_IIa"])
# qqnorm(test[test$Session=="LC","Percent_IIa"]);qqline(test[test$Session=="LC","Percent_IIa"])
# qqnorm(test[test$Session=="ME","Percent_IIa"]);qqline(test[test$Session=="ME","Percent_IIa"])

model <-lme(Percent_IIa~ Session, data=test[test$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
all_pvals["ANOVA","Percent_IIa"]<-a$`Pr(>Chisq)`

b<-
  t_test(test[test$Group=="BED REST",], Percent_IIa~Session, paired=TRUE)
b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","Percent_IIa"]<-b$p
# all_pvals["CON-ME","Percent_IIa"]<-a$p.adj[2]
# all_pvals["CON-LC","Percent_IIa"]<-a$p.adj[1]
# all_pvals["LC-ME","Percent_IIa"]<-a$p.adj[3]
all_pvals["BED_REST_test","Percent_IIa"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","Percent_IIa"]<-"anova_no_post_hoc"



# Type IIa/IIx

# qqnorm(test[test$Session=="BDC","Percent_IIa_IIx"]);qqline(test[test$Session=="BDC","Percent_IIa_IIx"])
# qqnorm(test[test$Session=="HDT55","Percent_IIa_IIx"]);qqline(test[test$Session=="HDT55","Percent_IIa_IIx"])
# qqnorm(test[test$Session=="CON","Percent_IIa_IIx"]);qqline(test[test$Session=="CON","Percent_IIa_IIx"])
# qqnorm(test[test$Session=="LC","Percent_IIa_IIx"]);qqline(test[test$Session=="LC","Percent_IIa_IIx"])
# qqnorm(test[test$Session=="ME","Percent_IIa_IIx"]);qqline(test[test$Session=="ME","Percent_IIa_IIx"])

model <-lme(Percent_IIa_IIx~ Session, data=test[test$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
all_pvals["ANOVA","Percent_IIa_IIx"]<-a$`Pr(>Chisq)`
a<-
  test[test$Group=="POST-VIRAL",] %>% tukey_hsd(Percent_IIa_IIx~Session)
b<-
  t_test(test[test$Group=="BED REST",], Percent_IIa_IIx~Session, paired=TRUE)
b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","Percent_IIa_IIx"]<-b$p
all_pvals["CON-ME","Percent_IIa_IIx"]<-a$p.adj[2]
all_pvals["CON-LC","Percent_IIa_IIx"]<-a$p.adj[1]
all_pvals["LC-ME","Percent_IIa_IIx"]<-a$p.adj[3]
all_pvals["BED_REST_test","Percent_IIa_IIx"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","Percent_IIa_IIx"]<-"anova_tukey_post_hoc"


# Type IIx 

# qqnorm(test[test$Session=="BDC","Percent_IIx"]);qqline(test[test$Session=="BDC","Percent_IIx"])
# qqnorm(test[test$Session=="HDT55","Percent_IIx"]);qqline(test[test$Session=="HDT55","Percent_IIx"])
# qqnorm(test[test$Session=="CON","Percent_IIx"]);qqline(test[test$Session=="CON","Percent_IIx"])
# qqnorm(test[test$Session=="LC","Percent_IIx"]);qqline(test[test$Session=="LC","Percent_IIx"])
# qqnorm(test[test$Session=="ME","Percent_IIx"]);qqline(test[test$Session=="ME","Percent_IIx"])

model <-lme(Percent_IIx~ Session, data=test[test$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")
a<-car::Anova(model)
all_pvals["ANOVA","Percent_IIx"]<-a$`Pr(>Chisq)`
b<-
  t_test(test[test$Group=="BED REST",], Percent_IIx~Session, paired=TRUE)
b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","Percent_IIx"]<-b$p
# all_pvals["CON-ME","Percent_IIx"]<-a$p.adj[2]
# all_pvals["CON-LC","Percent_IIx"]<-a$p.adj[1]
# all_pvals["LC-ME","Percent_IIx"]<-a$p.adj[3]
all_pvals["BED_REST_test","Percent_IIx"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","Percent_IIx"]<-"anova_no_post_hoc"



# type IIa/IIx + IIx sum

# qqnorm(test[test$Session=="BDC","Percent_IIa_IIx_IIx"]);qqline(test[test$Session=="BDC","Percent_IIa_IIx_IIx"])
# qqnorm(test[test$Session=="HDT55","Percent_IIa_IIx_IIx"]);qqline(test[test$Session=="HDT55","Percent_IIa_IIx_IIx"])
# qqnorm(test[test$Session=="CON","Percent_IIa_IIx_IIx"]);qqline(test[test$Session=="CON","Percent_IIa_IIx_IIx"])
# qqnorm(test[test$Session=="LC","Percent_IIa_IIx_IIx"]);qqline(test[test$Session=="LC","Percent_IIa_IIx_IIx"])
# qqnorm(test[test$Session=="ME","Percent_IIa_IIx_IIx"]);qqline(test[test$Session=="ME","Percent_IIa_IIx_IIx"])

a<-
  kruskal.test(Percent_IIa_IIx_IIx~Session, data=test[test$Group=="POST-VIRAL",])
all_pvals["ANOVA","Percent_IIa_IIx_IIx"]<-a$p.value

a<-pairwise.wilcox.test(test[test$Group=="POST-VIRAL","Percent_IIa_IIx_IIx"], test[test$Group=="POST-VIRAL","Session"], p.adjust.method="BH")

b<-
  wilcox.test(test[test$Session=="BDC","Percent_IIa_IIx_IIx"], test[test$Session=="HDT55","Percent_IIa_IIx_IIx"],paired=TRUE)

all_pvals["BDC-HDT55","Percent_IIa_IIx_IIx"]<-b$p.value
all_pvals["CON-ME","Percent_IIa_IIx_IIx"]<-a$p.value[2]
all_pvals["CON-LC","Percent_IIa_IIx_IIx"]<-a$p.value[1]
all_pvals["LC-ME","Percent_IIa_IIx_IIx"]<-a$p.value[4]
all_pvals["BED_REST_test","Percent_IIa_IIx_IIx"]<-"wilcoxon"
all_pvals["POST_VIRAL_test","Percent_IIa_IIx_IIx"]<-"kruskal_pw_wilcoxon_post_hoc"





stat_test <- tibble::tribble(
  ~group1, ~group2, ~p.adj.signif, ~Group, ~Fibre_type, ~y.coord, ~x_min,~x_max, ~tip1,~tip2, ~label_size, ~bracket_size, ~colours,
  "CON","LC",paste(format(round(as.numeric(all_pvals["CON-LC","Percent_I"]),3), drop0trailing=F)),"POST-VIRAL",1  , 105, 0.75, 1,   0.01, 0.5,  14, 6,"black",
  "CON","ME","p<0.001","POST-VIRAL",1  ,115 , 0.75, 1.25,0.01, 0.5,  14, 6,"black",
  "LC","ME",paste(format(round(as.numeric(all_pvals["LC-ME","Percent_I"]),3), drop0trailing=F)),"POST-VIRAL",1  , 95, 1, 1.25,     0.01, 0.5,  14, 6,"black",
  "CON","ME","p<0.001","POST-VIRAL",4  , 75, 3.75, 4.25, 0.01, 0.5,  14, 6,"black",
  "CON","LC",paste(format(round(as.numeric(all_pvals["CON-LC","Percent_IIa_IIx_IIx"]),3), drop0trailing=F)),"POST-VIRAL",4  , 60, 3.75, 4, 0.01, 0.5,  14, 6,"black")

data<-as.data.frame(read_xlsx("Manuscript_data_clean_230525.xlsx", sheet="Fig2a"))

data<-data %>%
  mutate(Group=dplyr::recode(Group, "AGBRESA"="BED REST","MUSCLE-ME"="POST-VIRAL"))
data<-data %>%
  mutate(Fibre_type=dplyr::recode(Fibre_type, "I"=1,"I_IIa"=2,"IIa"=3,"IIa_IIx"=5,"IIx"=6,"IIa_IIx_IIx"=4))

data<-data[data$Fibre_type!=5 & data$Fibre_type!=6,]

remove<-data[data$Fibre_type==1 &is.na(data$Percent)==TRUE,"Subject"]

data<-data[!data$Subject %in% remove,]

data[["Session"]] = factor( data[["Session"]], levels = c("BDC", "HDT55", "CON","LC","ME"))

data1<-data[data$Group=="BED REST",]
data2<-data[data$Group=="POST-VIRAL",]

fibre_type_a<-
  ggplot(data=data1, aes(x=Fibre_type, y=Percent))+
  geom_boxplot(data=data1[data1$Fibre_type==1 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.5),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data2[,"Session"])))+
  geom_boxplot(data=data1[data1$Fibre_type==2 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.5),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data2[,"Session"])))+
  geom_boxplot(data=data1[data1$Fibre_type==3 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.5),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data2[,"Session"])))+
  geom_boxplot(data=data1[data1$Fibre_type==4 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.5),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data2[,"Session"])))+
  geom_point(size=5,stroke=2 , shape=21, 
             position = position_jitterdodge(dodge.width = 0.5, jitter.height = 0, jitter.width = 0.2, seed=1), 
             aes(group=factor(Session),colour=factor(Session)), fill="white")+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"="grey42", "LC"=colours3[3], "ME"=colours3[4]))+
  scale_colour_manual(values=c("CON"="black","BDC"="black", "HDT55"="black", "LC"="black", "ME"="black"))+
  ylab(expression("% of Fibers"))+
  xlab("Fiber Type")+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size/2),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)), 
        axis.title.y = element_text(angle = 90, vjust = 1, hjust=0.3), # for atop functions export as 9.5x14in
        axis.title.x = element_text(angle = 0, vjust = 1, hjust=0.5),
        axis.text.y = element_text( size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( size =rel((base_size+axis_text_rel_size-10)/base_size)),
        strip.background = element_blank(),
        strip.text = element_blank())+
  scale_y_continuous(limits=c(0,140), breaks=seq(0,100,25),expand=c(0,0))+
  scale_x_continuous(breaks=c(1,2,3,4), labels = c("I","I/IIa","IIa","IIa/IIx + IIx"))+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=fibre_type_a,
       filename = paste(output_folder,"/fibre_type_a_",date,".png", sep = ""),
       device="png",  width = 17, height = 12, units = "in")

fibre_type_b<-
  ggplot(data=data2, aes(x=Fibre_type, y=Percent))+
  geom_boxplot(data=data2[data2$Fibre_type==1 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.75),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data1[,"Session"])))+
  geom_boxplot(data=data2[data2$Fibre_type==2 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.75),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data1[,"Session"])))+
  geom_boxplot(data=data2[data2$Fibre_type==3 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.75),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data1[,"Session"])))+
  geom_boxplot(data=data2[data2$Fibre_type==4 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.75),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data1[,"Session"])))+
  geom_boxplot(data=data2[data2$Fibre_type==5 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.75),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data1[,"Session"])))+
  geom_point(size=5,stroke=2 , shape=21, 
             position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0, jitter.width = 0.2, seed=1), 
             aes(group=factor(Session),colour=factor(Session)), fill="white")+
  stat_pvalue_manual(data=stat_test, label = "p.adj.signif",y.position="y.coord"
                     ,label.size=8.5,tip.length = 0.01,
                     xmin="x_min",xmax="x_max",
                     linetype="solid", 
                     bracket.size =2, 
                     inherit.aes=FALSE)+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"="grey42", "LC"=colours3[3], "ME"=colours3[4]))+
  scale_colour_manual(values=c("CON"="black","BDC"="black", "HDT55"="black", "LC"="black", "ME"="black"))+
  ylab(expression("% of Fibers"))+
  xlab("Fiber Type")+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size/2),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)), 
        axis.title.y = element_text(angle = 90, vjust = 1, hjust=0.3), # for atop functions export as 9.5x14in
        axis.title.x = element_text(angle = 0, vjust = 1, hjust=0.5),
        axis.text.y = element_text( size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( size =rel((base_size+axis_text_rel_size-10)/base_size)),
        strip.background = element_blank(),
        strip.text = element_blank())+
  scale_y_continuous(limits=c(0,140), breaks=seq(0,100,25),expand=c(0,0))+
  scale_x_continuous(breaks=c(1,2,3,4), labels = c("I","I/IIa","IIa","IIa/IIx + IIx"))+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=fibre_type_b,
       filename = paste(output_folder,"/fibre_type_b_",date,".png", sep = ""),
       device="png",  width = 17, height = 12, units = "in")




# Supplemental 3 ----------------------------------------------------------
# fibre type FCSA ---------------------------------------------------------


data<-as.data.frame(read_xlsx("Manuscript_data_clean_230525.xlsx", sheet="Sheet2"))

data<-data %>%
  mutate(Group=dplyr::recode(Group, "AGBRESA"="BED REST","MUSCLE-ME"="POST-VIRAL"))

data$Percent_IIa_IIx_IIx<-data$Percent_IIa_IIx+data$Percent_IIx

# Check normality ---------------------------------------------------------

op<-par(mfrow=c(1,5), mar=c(2,2,2,2))
sessions<-unique(data$Session)
for (i in 4:length(colnames(data)) ){
  try({
    parameter<-colnames(data[i])
    for (j in 1:length(sessions)){
      session<-sessions[j]
      tmp<- data[data$Session==session,parameter]
      qqnorm(tmp,xlab=parameter, main=paste(session,parameter,sep="_"))
      qqline(tmp)
      
    }
    
  },silent=TRUE)
  
  
}

norm_p_vals<-normality_test(data)
# FCSA

test<-as.data.frame(box_cox_transform(data))
test<-test %>% 
  filter_all(all_vars(!is.infinite(.)))
test[,c("Subject","Session","Group","Sex")]<-data[,c("Subject","Session","Group","Sex")]
test_norm<-normality_test(test)
test2<-test[!test$Subject %in% remove,]


# mean(data[data$Session=="BDC","TypeI_FCSA"],na.rm=TRUE); sd(data[data$Session=="BDC","TypeI_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","TypeI_FCSA"],na.rm=TRUE); sd(data[data$Session=="HDT55","TypeI_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="CON","TypeI_FCSA"],na.rm=TRUE); sd(data[data$Session=="CON","TypeI_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="LC","TypeI_FCSA"],na.rm=TRUE); sd(data[data$Session=="LC","TypeI_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="ME","TypeI_FCSA"],na.rm=TRUE); sd(data[data$Session=="ME","TypeI_FCSA"],na.rm=TRUE)
# 
# mean(data[data$Session=="BDC","TypeI_IIa_FCSA"],na.rm=TRUE); sd(data[data$Session=="BDC","TypeI_IIa_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","TypeI_IIa_FCSA"],na.rm=TRUE); sd(data[data$Session=="HDT55","TypeI_IIa_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="CON","TypeI_IIa_FCSA"],na.rm=TRUE); sd(data[data$Session=="CON","TypeI_IIa_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="LC","TypeI_IIa_FCSA"],na.rm=TRUE); sd(data[data$Session=="LC","TypeI_IIa_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="ME","TypeI_IIa_FCSA"],na.rm=TRUE); sd(data[data$Session=="ME","TypeI_IIa_FCSA"],na.rm=TRUE)
# 
# mean(data[data$Session=="BDC","TypeIIa_FCSA"],na.rm=TRUE); sd(data[data$Session=="BDC","TypeIIa_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","TypeIIa_FCSA"],na.rm=TRUE); sd(data[data$Session=="HDT55","TypeIIa_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="CON","TypeIIa_FCSA"],na.rm=TRUE); sd(data[data$Session=="CON","TypeIIa_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="LC","TypeIIa_FCSA"],na.rm=TRUE); sd(data[data$Session=="LC","TypeIIa_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="ME","TypeIIa_FCSA"],na.rm=TRUE); sd(data[data$Session=="ME","TypeIIa_FCSA"],na.rm=TRUE)
# 
# mean(data[data$Session=="BDC","Type_IIa_IIx_FCSA"],na.rm=TRUE); sd(data[data$Session=="BDC","Type_IIa_IIx_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","Type_IIa_IIx_FCSA"],na.rm=TRUE); sd(data[data$Session=="HDT55","Type_IIa_IIx_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="CON","Type_IIa_IIx_FCSA"],na.rm=TRUE); sd(data[data$Session=="CON","Type_IIa_IIx_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="LC","Type_IIa_IIx_FCSA"],na.rm=TRUE); sd(data[data$Session=="LC","Type_IIa_IIx_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="ME","Type_IIa_IIx_FCSA"],na.rm=TRUE); sd(data[data$Session=="ME","Type_IIa_IIx_FCSA"],na.rm=TRUE)
# 
# mean(data[data$Session=="BDC","TypeIIx_FCSA"],na.rm=TRUE); sd(data[data$Session=="BDC","TypeIIx_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="HDT55","TypeIIx_FCSA"],na.rm=TRUE); sd(data[data$Session=="HDT55","TypeIIx_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="CON","TypeIIx_FCSA"],na.rm=TRUE); sd(data[data$Session=="CON","TypeIIx_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="LC","TypeIIx_FCSA"],na.rm=TRUE); sd(data[data$Session=="LC","TypeIIx_FCSA"],na.rm=TRUE)
# mean(data[data$Session=="ME","TypeIIx_FCSA"],na.rm=TRUE); sd(data[data$Session=="ME","TypeIIx_FCSA"],na.rm=TRUE)

# Type I
# qqnorm(test[test$Session=="BDC","TypeI_FCSA"]);qqline(test[test$Session=="BDC","TypeI_FCSA"])
# qqnorm(test[test$Session=="HDT55","TypeI_FCSA"]);qqline(test[test$Session=="HDT55","TypeI_FCSA"])
# qqnorm(test[test$Session=="CON","TypeI_FCSA"]);qqline(test[test$Session=="CON","TypeI_FCSA"])
# qqnorm(test[test$Session=="LC","TypeI_FCSA"]);qqline(test[test$Session=="LC","TypeI_FCSA"])
# qqnorm(test[test$Session=="ME","TypeI_FCSA"]);qqline(test[test$Session=="ME","TypeI_FCSA"])

model <-lme(TypeI_FCSA~ Session, data=test[test$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")

a<-car::Anova(model)
all_pvals["ANOVA","TypeI_FCSA"]<-a$`Pr(>Chisq)`
# post-hoc testing    
a<-
  test[test$Group=="POST-VIRAL",] %>% tukey_hsd(TypeI_FCSA~Session)   

b<-
  t_test(test[test$Group=="BED REST",], TypeI_FCSA~Session, paired=TRUE)
b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","TypeI_FCSA"]<-b$p
all_pvals["CON-ME","TypeI_FCSA"]<-a$p.adj[2]
all_pvals["CON-LC","TypeI_FCSA"]<-a$p.adj[1]
all_pvals["LC-ME","TypeI_FCSA"]<-a$p.adj[3]
all_pvals["BED_REST_test","TypeI_FCSA"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","TypeI_FCSA"]<-"anova_tukey_post_hoc"


# Type I/IIa
# qqnorm(test[test$Session=="BDC","TypeI_IIa_FCSA"]);qqline(test[test$Session=="BDC","TypeI_IIa_FCSA"])
# qqnorm(test[test$Session=="HDT55","TypeI_IIa_FCSA"]);qqline(test[test$Session=="HDT55","TypeI_IIa_FCSA"])
# qqnorm(test[test$Session=="CON","TypeI_IIa_FCSA"]);qqline(test[test$Session=="CON","TypeI_IIa_FCSA"])
# qqnorm(test[test$Session=="LC","TypeI_IIa_FCSA"]);qqline(test[test$Session=="LC","TypeI_IIa_FCSA"])
# qqnorm(test[test$Session=="ME","TypeI_IIa_FCSA"]);qqline(test[test$Session=="ME","TypeI_IIa_FCSA"])


model <-lme(TypeI_IIa_FCSA~ Session, data=test[test$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")

a<-car::Anova(model)
all_pvals["ANOVA","TypeI_IIa_FCSA"]<-a$`Pr(>Chisq)`

b<-
  t_test(test[test$Group=="BED REST",], TypeI_IIa_FCSA~Session, paired=TRUE)
b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","TypeI_IIa_FCSA"]<-b$p
# all_pvals["CON-ME","TypeI_IIa_FCSA"]<-a$p.adj[2]
# all_pvals["CON-LC","TypeI_IIa_FCSA"]<-a$p.adj[1]
# all_pvals["LC-ME","TypeI_IIa_FCSA"]<-a$p.adj[3]
all_pvals["BED_REST_test","TypeI_IIa_FCSA"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","TypeI_IIa_FCSA"]<-"anova_no_post_hoc"


# Type IIa  
# qqnorm(test[test$Session=="BDC","TypeIIa_FCSA"]);qqline(test[test$Session=="BDC","TypeIIa_FCSA"])
# qqnorm(test[test$Session=="HDT55","TypeIIa_FCSA"]);qqline(test[test$Session=="HDT55","TypeIIa_FCSA"])
# qqnorm(test[test$Session=="CON","TypeIIa_FCSA"]);qqline(test[test$Session=="CON","TypeIIa_FCSA"])
# qqnorm(test[test$Session=="LC","TypeIIa_FCSA"]);qqline(test[test$Session=="LC","TypeIIa_FCSA"])
# qqnorm(test[test$Session=="ME","TypeIIa_FCSA"]);qqline(test[test$Session=="ME","TypeIIa_FCSA"])

model <-lme(TypeIIa_FCSA~ Session, data=test[test$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")

a<-car::Anova(model)
all_pvals["ANOVA","TypeIIa_FCSA"]<-a$`Pr(>Chisq)`
b<-
  t_test(test[test$Group=="BED REST",], TypeIIa_FCSA~Session, paired=TRUE)
b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))


all_pvals["BDC-HDT55","TypeIIa_FCSA"]<-b$p
# all_pvals["CON-ME","TypeIIa_FCSA"]<-a$p.adj[2]
# all_pvals["CON-LC","TypeIIa_FCSA"]<-a$p.adj[1]
# all_pvals["LC-ME","TypeIIa_FCSA"]<-a$p.adj[3]
all_pvals["BED_REST_test","TypeIIa_FCSA"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","TypeIIa_FCSA"]<-"anova_no_post_hoc"

# Type IIa/IIx 
# qqnorm(test[test$Session=="BDC","Type_IIa_IIx_FCSA"]);qqline(test[test$Session=="BDC","Type_IIa_IIx_FCSA"])
# qqnorm(test[test$Session=="HDT55","Type_IIa_IIx_FCSA"]);qqline(test[test$Session=="HDT55","Type_IIa_IIx_FCSA"])
# qqnorm(test[test$Session=="CON","Type_IIa_IIx_FCSA"]);qqline(test[test$Session=="CON","Type_IIa_IIx_FCSA"])
# qqnorm(test[test$Session=="LC","Type_IIa_IIx_FCSA"]);qqline(test[test$Session=="LC","Type_IIa_IIx_FCSA"])
# qqnorm(test[test$Session=="ME","Type_IIa_IIx_FCSA"]);qqline(test[test$Session=="ME","Type_IIa_IIx_FCSA"])

model <-lme(Type_IIa_IIx_FCSA~ Session, data=test[test$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")

a<-car::Anova(model)
all_pvals["ANOVA","Type_IIa_IIx_FCSA"]<-a$`Pr(>Chisq)`

b<-
  t_test(test[test$Group=="BED REST",], Type_IIa_IIx_FCSA~Session, paired=TRUE)
b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","Type_IIa_IIx_FCSA"]<-b$p
# all_pvals["CON-ME","Type_IIa_IIx_FCSA"]<-a$p.adj[2]
# all_pvals["CON-LC","Type_IIa_IIx_FCSA"]<-a$p.adj[1]
# all_pvals["LC-ME","Type_IIa_IIx_FCSA"]<-a$p.adj[3]
all_pvals["BED_REST_test","Type_IIa_IIx_FCSA"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","Type_IIa_IIx_FCSA"]<-"anova_no_post_hoc"


# Type IIx 

# qqnorm(test[test$Session=="BDC","TypeIIx_FCSA"]);qqline(test[test$Session=="BDC","TypeIIx_FCSA"])
# qqnorm(test[test$Session=="HDT55","TypeIIx_FCSA"]);qqline(test[test$Session=="HDT55","TypeIIx_FCSA"])
# qqnorm(test[test$Session=="CON","TypeIIx_FCSA"]);qqline(test[test$Session=="CON","TypeIIx_FCSA"])
# qqnorm(test[test$Session=="LC","TypeIIx_FCSA"]);qqline(test[test$Session=="LC","TypeIIx_FCSA"])
# qqnorm(test[test$Session=="ME","TypeIIx_FCSA"]);qqline(test[test$Session=="ME","TypeIIx_FCSA"])

model <-lme(TypeIIx_FCSA~ Session, data=test[test$Group=="POST-VIRAL" ,], random = ~ 1|Subject, na.action = na.omit, control="optim")

a<-car::Anova(model)
all_pvals["ANOVA","TypeIIx_FCSA"]<-a$`Pr(>Chisq)`

b<-
  t_test(test[test$Group=="BED REST",], TypeIIx_FCSA~Session, paired=TRUE)
b<-add_significance(b, cutpoints = c(0,1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

all_pvals["BDC-HDT55","TypeIIx_FCSA"]<-b$p
# all_pvals["CON-ME","TypeIIx_FCSA"]<-a$p.adj[2]
# all_pvals["CON-LC","TypeIIx_FCSA"]<-a$p.adj[1]
# all_pvals["LC-ME","TypeIIx_FCSA"]<-a$p.adj[3]
all_pvals["BED_REST_test","TypeIIx_FCSA"]<-"paired_t_test"
all_pvals["POST_VIRAL_test","TypeIIx_FCSA"]<-"anova_no_post_hoc"




stat_test_a <- tibble::tribble(
  ~group1, ~group2, ~p.adj.signif, ~Group, ~Fibre_type, ~y.coord, ~x_min,~x_max, ~tip1,~tip2, ~label_size, ~bracket_size, ~colours,
  "BDC","HDT55",paste(format(round(as.numeric(all_pvals["BDC-HDT55","TypeI_FCSA"]),3), drop0trailing=F)),"BED REST",1  , 6700, 0.75, 1.25,0.01, 0.01, 30, 10, "black",
  "BDC","HDT55",paste(format(round(as.numeric(all_pvals["BDC-HDT55","TypeI_IIa_FCSA"]),3), drop0trailing=F)),"BED REST",2  , 8300, 1.75, 2.25,0.01, 0.01, 30, 10, "black",
  "BDC","HDT55",paste(format(round(as.numeric(all_pvals["BDC-HDT55","TypeIIa_FCSA"]),3), drop0trailing=F)),"BED REST",3  , 8500, 2.75, 3.25,0.01, 0.01, 30, 10, "black",
  "BDC","HDT55",paste(format(round(as.numeric(all_pvals["BDC-HDT55","Type_IIa_IIx_FCSA"]),3), drop0trailing=F)),"BED REST",4  , 8300, 3.75, 4.25,0.01, 0.01, 30, 10, "black")


stat_test_b <- tibble::tribble(
  ~group1, ~group2, ~p.adj.signif, ~Group, ~Fibre_type, ~y.coord, ~x_min,~x_max, ~tip1,~tip2, ~label_size, ~bracket_size, ~colours,
  "CON","ME",paste(format(round(as.numeric(all_pvals["CON-ME","TypeI_FCSA"]),3), drop0trailing=F)),"POST-VIRAL",1  , 9200, 0.75, 1.25,0.01, 0.01, 30, 10, "black",
  "CON","LC",paste(format(round(as.numeric(all_pvals["CON-LC","TypeI_FCSA"]),3), drop0trailing=F)),"POST-VIRAL",1  , 8500, 0.75, 1.0,0.01, 0.01, 30, 10, "black",
  "LC","ME",paste(format(round(as.numeric(all_pvals["LC-ME","TypeI_FCSA"]),3), drop0trailing=F)),"POST-VIRAL",1  , 7800, 1.0, 1.25,0.01, 0.01, 30, 10, "black")


data<-as.data.frame(read_xlsx("Manuscript_data_clean_230525.xlsx", sheet="Fig2a"))

data<-data %>%
  mutate(Group=dplyr::recode(Group, "AGBRESA"="BED REST","MUSCLE-ME"="POST-VIRAL"))


data<-data %>%
  mutate(Fibre_type=dplyr::recode(Fibre_type, "I"=1,"I_IIa"=2,"IIa"=3,"IIa_IIx"=4,"IIx"=5))

data[["Session"]] = factor( data[["Session"]], levels = c("BDC", "HDT55", "CON","LC","ME"))

data1<-data[data$Group=="BED REST",]
data2<-data[data$Group=="POST-VIRAL",]

fibre_type_FCSA_a<-
  ggplot(data=data1, aes(x=Fibre_type, y=FCSA))+
  geom_boxplot(data=data1[data1$Fibre_type==1 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.5),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data2[,"Session"])))+
  geom_boxplot(data=data1[data1$Fibre_type==2 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.5),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data2[,"Session"])))+
  geom_boxplot(data=data1[data1$Fibre_type==3 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.5),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data2[,"Session"])))+
  geom_boxplot(data=data1[data1$Fibre_type==4 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.5),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data2[,"Session"])))+
  geom_boxplot(data=data1[data1$Fibre_type==5 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.5),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data2[,"Session"])))+
  geom_point(size=5,stroke=2 , shape=21, 
             position = position_jitterdodge(dodge.width = 0.5, jitter.height = 0, jitter.width = 0.2, seed=1), 
             aes(group=factor(Session),colour=factor(Session)), fill="white")+
  stat_pvalue_manual(data=stat_test_a, label = "p.adj.signif",y.position="y.coord"
                     ,label.size=8.5,tip.length = 0.01,
                     xmin="x_min",xmax="x_max",
                     linetype="solid",
                     bracket.size =2,
                     inherit.aes=FALSE)+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"="grey42", "LC"=colours3[3], "ME"=colours3[4]))+
  scale_colour_manual(values=c("CON"="black","BDC"="black", "HDT55"="black", "LC"="black", "ME"="black"))+
  ylab(expression("FCSA ("*mu*m^2*")"))+
  xlab("Fiber Type")+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size/2),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)), 
        axis.title.y = element_text(angle = 90, vjust = 1, hjust=0.3), # for atop functions export as 9.5x14in
        axis.title.x = element_text(angle = 0, vjust = 1, hjust=0.5),
        axis.text.y = element_text( size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( size =rel((base_size+axis_text_rel_size-10)/base_size)),
        strip.background = element_blank(),
        strip.text = element_blank())+
  scale_y_continuous(limits=c(0,10500), breaks=seq(0,10000,2000),expand=c(0,0))+
  scale_x_continuous(breaks=c(1,2,3,4,5), labels = c("I","I/IIa","IIa","IIa/IIx","IIx"))+
  guides(y=guide_axis(cap="upper"))

ggsave(plot=fibre_type_FCSA_a,
       filename = paste(output_folder,"/fibre_type_FCSA_a_",date,".png", sep = ""),
       device="png",  width = 17, height = 12, units = "in")

fibre_type_FCSA_b<-
  ggplot(data=data2, aes(x=Fibre_type, y=FCSA))+
  geom_boxplot(data=data2[data2$Fibre_type==1 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.75),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data1[,"Session"])))+
  geom_boxplot(data=data2[data2$Fibre_type==2 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.75),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data1[,"Session"])))+
  geom_boxplot(data=data2[data2$Fibre_type==3 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.75),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data1[,"Session"])))+
  geom_boxplot(data=data2[data2$Fibre_type==4 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.75),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data1[,"Session"])))+
  geom_boxplot(data=data2[data2$Fibre_type==5 ,],linewidth=2, outlier.shape = NA, coef=0,
               position = position_dodge(width = 0.75),
               aes(group=factor(Session),fill=factor(Session)),
               width=1.5/length(unique(data1[,"Session"])))+
  geom_point(size=5,stroke=2 , shape=21, 
             position = position_jitterdodge(dodge.width = 0.75, jitter.height = 0, jitter.width = 0.2, seed=1), 
             aes(group=factor(Session),colour=factor(Session)), fill="white")+
  stat_pvalue_manual(data=stat_test_b, label = "p.adj.signif",y.position="y.coord"
                     ,label.size=8.5,tip.length = 0.01,
                     xmin="x_min",xmax="x_max",
                     linetype="solid", 
                     bracket.size =2, 
                     inherit.aes=FALSE)+
  scale_fill_manual(values = c("CON"=colours3[1],"BDC"=colours3[1], "HDT55"="grey42", "LC"=colours3[3], "ME"=colours3[4]))+
  scale_colour_manual(values=c("CON"="black","BDC"="black", "HDT55"="black", "LC"="black", "ME"="black"))+
  ylab(expression("FCSA ("*mu*m^2*")"))+
  xlab("Fiber Type")+
  theme(legend.position = "none",
        axis.line=element_line(colour="black", size = line_size/2),
        axis.ticks = element_blank(),
        axis.ticks.x =element_blank(),
        text=element_text(size=base_size, family = "Arial"),
        panel.background  = element_rect(fill="white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(
          size = rel((title_text_rel_size + base_size) / base_size),
          hjust = 0.5
        ),
        axis.title = element_text(size = rel((title_text_rel_size + base_size) / base_size)), 
        axis.title.y = element_text(angle = 90, vjust = 1, hjust=0.3), # for atop functions export as 9.5x14in
        axis.title.x = element_text(angle = 0, vjust = 1, hjust=0.5),
        axis.text.y = element_text( size =rel((base_size+axis_text_rel_size)/base_size)),
        axis.text.x = element_text( size =rel((base_size+axis_text_rel_size-10)/base_size)),
        strip.background = element_blank(),
        strip.text = element_blank())+
  scale_y_continuous(limits=c(0,10500), breaks=seq(0,10000,2000),expand=c(0,0))+
  scale_x_continuous(breaks=c(1,2,3,4,5), labels = c("I","I/IIa","IIa","IIa/IIx","IIx"))+
  guides(y=guide_axis(cap="upper"))


ggsave(plot=fibre_type_FCSA_b,
       filename = paste(output_folder,"/fibre_type_FCSA_b_",date,".png", sep = ""),
       device="png",  width = 17, height = 12, units = "in")





# save csv files -------------------------------------------------------------

write.csv(all_pvals, file =paste(output_folder,"/all_pvals_",date,".csv", sep = "") )
write.csv(cor_pvals, file =paste(output_folder,"/cor_pvals_",date,".csv", sep = "") )
write.csv(HLM_pvals, file =paste(output_folder,"/HLM_pvals_",date,".csv", sep = "") )
write.csv(BC_pvals, file =paste(output_folder,"/BC_pvals_",date,".csv", sep = "") )
write.csv(test_norm, file =paste(output_folder,"/test_norm_",date,".csv", sep = "") )
