# -----------------------------------------
# droso_trans_prob_FE. R
#   implement transitional probabilities(Markov-like) methods to analysis 
#   drosophilia gag gene network
#   database: FlyEx
# -----------------------------------------
# different section of the codes:
#           I. setup and loading dataset 
#           II. cleaning up data & analysis
#           III. simulation 
#           IV. spacial decomposion
# -----------------------------------------
# Author:Jianhong Chen
# Date: Jan 20th 2019
# -----------------------------------------


# load required packages

library('ggplot2')
library('tidyr')
library('dplyr')
library('tidyverse')
library('purrr')
library('reshape2')

# -----------------------------------------------------------------------------------------------------
# reset all files and created variables
rm(list = ls())
# -----------------------------------------------------------------------------------------------------


#########################################################################################
# I.  setup and loading dataset :
#########################################################################################
# -----------------------------------------------------------------------------------------------------
# load the FlyEx data files
# -----------------------------------------------------------------------------------------------------
# data for gt & kni:
dir1 = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Gap_Gene_of_Drosophila/FlyEx_data/gt_kni14At1-8_w_R3/txt/byEmbryos"
setwd(dir1)
f1 = list.files()
gt_kni_raw = map(f1, read.table, header = T, col.names = c('Nucleus_Num', 'AP', 'DV', 'eve', 'kni', 'gt'))

# data for kr & hb:
dir2 = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Gap_Gene_of_Drosophila/FlyEx_data/kr_hb14At1-8_w_R3/txt/byEmbryos"
setwd(dir2)
f2 = list.files()
kr_hb_raw = map(f2, read.table, header = T, col.names = c('Nucleus_Num', 'AP', 'DV', 'eve', 'kr', 'hb'))

# quick plot of the raw data without any processing
gt_kni_exp = data.frame()
kr_hb_exp = data.frame()
for (i in 1:8) {
  gt_kni_temp = gt_kni_raw[[i]] %>%
    mutate('time' = rep(i, dim(gt_kni_raw[[i]])[1]))
  gt_kni_exp = rbind(gt_kni_exp, gt_kni_temp)
  #gt_kni_exp = rbind(gt_kni_exp, gt_kni_temp)
  kr_hb_temp = kr_hb_raw[[i]] %>%
    mutate('time' = rep(i, dim(kr_hb_raw[[i]])[1]))
  kr_hb_exp = rbind(kr_hb_exp, kr_hb_temp)
}
gt_kni_raw_long = gt_kni_exp %>%
  select(c(AP, kni, gt, time)) %>%
  gather(key = gene, value = intensity, kni:gt)

kr_hb_raw_long = kr_hb_exp %>%
  select(c(AP, kr, hb, time)) %>%
  gather(key = gene, value = intensity, kr:hb)

ggplot(gt_kni_raw_long, aes(x = AP, y = intensity, color = gene))+
  geom_point(alpha = 0.8, shape = ".") + 
  facet_wrap(~time, nrow = 2) +
  ggtitle('Gap Gene Expression Raw FlyEx Data gt-kni') +
  theme(text = element_text(size = 15))

ggplot(kr_hb_raw_long, aes(x = AP, y = intensity, color = gene))+
  geom_point(alpha = 0.8, shape = ".") + 
  facet_wrap(~time, nrow = 2) +
  ggtitle('Gap Gene Expression Raw FlyEx Data kr-hb') +
  theme(text = element_text(size = 15))


#########################################################################################
# II.  cleanup & analysis :
#########################################################################################
# -----------------------------------------------------------------------------------------------------
# clean up data: converting into 100 bins
# -----------------------------------------------------------------------------------------------------
# converting into standarized AP position (100 points)
Strip_Data = function(x, sep,  gen1, gen2){
  # bin the intensity measurements of the gap genes into specified strips by taking the mean value:
  #   inputs: x = data file;
  #           sep = width of the strip
  #           gen1, gen2 = interested genes
  #   output: a dataframe of the binned (averaged) of the gap gene 
  oneAvg = data.frame()
  for (n in seq(0, (100-sep),sep)){
    one = subset(x, (x$AP > n) & (x$AP < (n+sep)))
    oneAvg = rbind(oneAvg, c(mean(one[[gen1]]), mean(one[[gen2]])))
  }
  colnames(oneAvg) = c(gen1, gen2)
  return(oneAvg)
}

gt_kni_bin = map(gt_kni_raw, Strip_Data, sep = 1, gen1 = "gt", gen2 = "kni")
kr_hb_bin = map(kr_hb_raw, Strip_Data, sep = 1, gen1 = "kr", gen2 = "hb")

# finally assemble the binned data into one complete dataframe
droso_gap_bin = data.frame()
for (i in 1:8){
  gap_temp = cbind('AP' = 4:98, 
                    gt_kni_bin[[i]][4:98, ], kr_hb_bin[[i]][4:98, ], 
                   'time' = rep(i, dim(gt_kni_bin[[i]][4:98, ])[1]))
  droso_gap_bin = rbind(droso_gap_bin, gap_temp)
}

# long df for binned intensity value
droso_gap_bin_long = droso_gap_bin %>%
  gather(key = gene, value = intensity, gt:hb)

# quick plot of the binned intensity results
ggplot(droso_gap_bin_long, aes(x = AP, y = intensity, color = gene)) +
  geom_line(alpha = 0.7) + 
  facet_wrap(~time, nrow = 2) +
  ggtitle('Gap Gene Expression Binned FlyEx Data') +
  theme(text = element_text(size = 15))

# -----------------------------------------------------------------------------------------------------
# converting intensity into Boolean Value (1/0 = on/off)
# -----------------------------------------------------------------------------------------------------
gene_mean = map(droso_gap_bin[, 2:5], mean)

droso_Boolean = droso_gap_bin %>%
  mutate(
    'hb_Boolean' = case_when(
      hb > gene_mean['hb'] ~ 1,
      hb < gene_mean['hb'] ~ 0), 
    'kr_Boolean' = case_when(
      kr > gene_mean['kr'] ~ 1,
      kr < gene_mean['kr'] ~ 0),
    'gt_Boolean' = case_when(
      gt > gene_mean['gt'] ~ 1,
      gt < gene_mean['gt'] ~ 0),
    'kni_Boolean' = case_when(
      kni > gene_mean['kni'] ~ 1,
      kni < gene_mean['kni'] ~ 0)
  ) %>%
  select(c(AP, time:kni_Boolean))

droso_Boolean_long = droso_Boolean %>%
  gather(key = gene, value = Boolean, hb_Boolean:kni_Boolean)
  
# quick plot of Boolean Value 
ggplot(droso_Boolean_long, aes(x= AP, y = Boolean, color = gene)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~time, nrow = 2) +
  ggtitle('Gap Gene Expression in Boolean') +
  theme(text = element_text(size = 15))

# -----------------------------------------------------------------------------------------------------
# create gene state: (gt, kr, kni, hb) 
#   calculate the transitional probabilitie between each different state
# Notice: Markov-like behaviro: only the previous timeframe can affects the next one
# -----------------------------------------------------------------------------------------------------

state_name = c('0 0 0 0', '1 0 0 0', '0 1 0 0', '0 0 1 0', '0 0 0 1',
               '1 1 0 0', '1 0 1 0', '1 0 0 1', '0 1 1 0', '0 1 0 1',
               '0 0 1 1', '1 1 1 0', '1 0 1 1', '1 1 0 1', '0 1 1 1',
               '1 1 1 1' ) # works similiar to Python Dictionary object in this case

TP_notation = function(x, dict = state_name) {
  # a function that conver state into 1-16 notation for convienience (built for 'lapply' or 'map' funciton)
  # input:  x: a given state
  #         dict: provided list that acts as Python dictionary object in this case for the given states
  return(which(x == state_name))
}

# 1. create the right dataframe for counting the transitional state
droso_state = droso_Boolean %>%
  mutate('gt' = as.character(gt_Boolean),
         'kr' = as.character(kr_Boolean),
         'kni' = as.character(kni_Boolean),
         'hb' = as.character(hb_Boolean)) %>%
  mutate('state' = paste(gt, kr, kni, hb)) %>%
  select(c(1, 2, 11)) 

droso_state_int = map_int(matrix(unlist(droso_state[, 3])), TP_notation) 
droso_state_nor = data_frame('time' = droso_state[, 2], 'state' = droso_state_int)
# make droso_state_string into a matrix format
droso_state_mat = subset(droso_state_nor, time == '1')
for (t in 2:8){
  droso_state_mat = cbind(droso_state_mat, subset(droso_state_nor, time == t))
}
droso_state_mat = droso_state_mat[, c(2,4,6,8,10,12,14,16)]
row.names(droso_state_mat) = 4:98 # row = AP position
colnames(droso_state_mat) = 1:8 # col = timeframe

# 2. count the number of the transition
trans_counter_mat = matrix(0, 16, 16) # row = initial_state, col = trans_state
for (x in 1:dim(droso_state_mat)[1]){
  for (t in 1:7){
    initial_state = droso_state_mat[x,t]
    trans_state = droso_state_mat[x, t+1]
    trans_counter_mat[initial_state, trans_state] = trans_counter_mat[initial_state, trans_state] + 1
  }
}

# 3. convert the results into probabilities
trans_prob_list = list()
for (i in 1:16){
  trans_prob_temp = map(trans_counter_mat[i, ], function(x){x/sum(trans_counter_mat[i, ])})
  trans_prob_list= append(trans_prob_list, trans_prob_temp)
}
trans_prob_mat = matrix(unlist(trans_prob_list), nrow = 16) # row = transitional state, col = initial state
trans_prob_mat[is.na(trans_prob_mat)] = 0 

# 4. visulizating the reulst
tp_long = melt(trans_prob_mat)
ggplot(tp_long, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  scale_fill_gradient(low="grey90", high="red") +
  xlab('Initial States') +
  ylab("Transitonal States") +
  ggtitle('Transitional Probabilities Matrix') +
  theme(text = element_text(size = 15))

ggplot(tp_long, aes(x = Var1, y = value)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Var2) +
  ylab('Transitional Probabilities') +
  xlab("Transitonal States") +
  ggtitle('Transitional Probabilities Barplot') +
  theme(text = element_text(size = 15))

#########################################################################################
# III.  Simulation:
#########################################################################################
# -----------------------------------------------------------------------------------------------------
# Simulate the gap gene interactions with the calculated transitional probabilities
# -----------------------------------------------------------------------------------------------------
trans_prob = as.data.frame(trans_prob_mat) # row = transitional state; col = initial state
trans_prob_cummulative = trans_prob_mat # cummlative version of the trans_prob for simulation purpose
for (init in 1:16){
  ind_zero = which(trans_prob_mat[, init] == 0) # record the index where the entry ==0
  for (tran in 1:15){
    trans_prob_cummulative[tran+1, init] = trans_prob_cummulative[tran, init] + trans_prob_cummulative[tran+1, init]
  }
  trans_prob_cummulative[ind_zero, init] = 0 # reset the original entries back to zero
}

TP_Simulation = function(s){
  # Scan through all states and determine the next state based on the calculated transitional probabilities
  #       input: s = each individual state
  #       output: next_state = state after applying transitional probabilities
  tp_one = trans_prob_cummulative[, s]
  for (i in tp_one){
    rd = runif(1) # generate random number
    if (rd < i){
      next_state = which(tp_one == i)
      break # it is important to break the for loop here to terminate scanning which change i 
    }
  }
  return(next_state)
} 

nn = 1000 # determine number of simulation

droso_sim_n_list  = sim_label = list()
for (n in 1:nn) {
  droso_sim_temp = matrix(0, 95, 8)
  droso_sim_temp[, 1] = droso_state_mat[, 1]
  for (t in 1:7){
    droso_sim_temp[, t+1] = map_int(droso_sim_temp[, t], TP_Simulation)
  }
  # renaming and reorginazing
  row.names(droso_sim_temp) = 4:98 
  droso_sim_temp_long = melt(droso_sim_temp)
  # decoding the state back to (gt, kr, kni, hb):
  droso_sim_de = data_frame('AP' = droso_sim_temp_long$Var1, 
                            'time' = droso_sim_temp_long$Var2, 
                            'state' = state_name[droso_sim_temp_long$value])
  state_split = strsplit(droso_sim_de$state, split =' ')
  gt = kr = kni = hb = list()
  for (i in state_split){
    gt = unlist(append(gt, as.integer(i[1])))
    kr = unlist(append(kr, as.integer(i[2])))
    kni = unlist(append(kni, as.integer(i[3])))
    hb = unlist(append(hb, as.integer(i[4])))
  }
  droso_sim_final = data.frame(droso_sim_de) %>%
    mutate(gt, kr, kni, hb) %>% 
    gather(key = 'gene', value = 'Boolean', gt:hb)
  # store each simulation into a large list
  droso_sim_n_list[[n]] = droso_sim_final
  # creating 'sim' vector for naming purpose
  sim_label[[n]] = rep(n, dim(droso_sim_final)[1])
}

# now concatenate all simulations into one huge matrix/df:
droso_sim_n = droso_sim_n_list[[1]]
for (n in 2:nn){ #start with 2, because the 1st one is used for initiaion of the matrix, 'droso_sim_n'
  droso_sim_n = rbind(droso_sim_n, droso_sim_n_list[[n]])
}
droso_sim_n = mutate(droso_sim_n, 'sim' = unlist(sim_label))

# Averaging all the simulation results into probabilities to eliminate the randomness

Avg_Simulation = function(ge, df, n){
  # average over all time and AP position for specific gene, created for map function mainly
  # input: ge = specific gene(string)
  #        df = dataframe/matrix of simluation results
  #        n = number of simulation
  # output: gene_avg_t = averaged results through all time and AP position(list)
  
  gene_avg_x = gene_avg_t = list()
  for (t in 1:8){
    for (x in 4:98){
      sim_subset = filter(df, AP == x, time ==t, gene == ge)
      gene_avg_x[[x]] = sum(sim_subset$Boolean)/n #avg at time t through all AP position
    }
    gene_avg_t[[t]] = unlist(gene_avg_x) # avg through all time and all AP
  }
  return(gene_avg_t)
}

droso_avgsim_all = map(c('gt', 'kr', 'kni', 'hb'), Avg_Simulation, df = droso_sim_n, n = nn)
droso_avgsim = data.frame('AP' = droso_Boolean$AP, 'time' = droso_Boolean$time) %>%
  mutate('gt_avg' = unlist(droso_avgsim_all[[1]]),
         'kr_avg' = unlist(droso_avgsim_all[[2]]),
         'kni_avg' = unlist(droso_avgsim_all[[3]]), 
         'hb_avg' = unlist(droso_avgsim_all[[4]])
  ) 

droso_avgsim_long = gather(droso_avgsim, key = 'gene', value = 'avg_prob', gt_avg:hb_avg)


# visulization the simulation results:
ggplot(droso_avgsim_long, aes(x = AP, y = avg_prob)) +
  geom_line() +
  facet_grid(time ~ gene) +
  ylab('Avg Probabilities') +
  xlab("AP position") +
  ggtitle('Simulation Results for n = 1000') +
  theme(text = element_text(size = 15))
                                              
# compare simulation and experimental Boolean results:
source = c(rep('sim', 3040), rep('WT', 3040))
gene_label = rep(c(rep('gt', 760), rep('kr', 760), rep('kni', 760), rep('hb', 760)), 2)
droso_compare = droso_avgsim %>%
  cbind(droso_Boolean[, c(5, 4, 6, 3)]) %>%
  gather(key = 'gene', value = 'value', gt_avg:hb_Boolean) %>%
  mutate('source' = source, 'gene_label' = gene_label)
  
droso_compare_gt = droso_compare[, c(1,2,3,9)] %>%
  gather(key = 'type', value = 'value', gt_avg:gt_Boolean)
  
ggplot(droso_compare, aes(x = AP, y = value)) +
  geom_line(alpha = 0.8, aes(linetype = source)) + 
  facet_grid(gene_label~time) +
  xlab("AP position") +
  ylab("Probability") + 
  ggtitle("Simulation vs Experimental") + 
  theme(text = element_text(size = 15))

# calculate the WT probabilities: intensity >> Boolean >> probabilities >> bin to 100

Probability_Conversion = function(x, sep,  gene1, gene2){
  # Convert the raw intensity measurements of the gap genes into 
  #     probabilities by counting the Boolean value then average over a strip:
  #   inputs: x = data file;
  #           sep = width of the strip
  #           gen1, gen2 = interested genes
  #   output: a dataframe of converted probabilities of the gap gene
  one_prob = data.frame()
  for (n in seq(0, (100-sep),sep)){
    one = subset(x, (x$AP > n) & (x$AP < (n+sep)))
    gene1_bool = map(one[[gene1]], function(x){case_when(x > gene_mean[gene1]~1, x < gene_mean[gene1] ~ 0)})
    gene2_bool = map(one[[gene2]], function(x){case_when(x > gene_mean[gene2]~1, x < gene_mean[gene2] ~ 0)})
    gene1_prob = sum(unlist(gene1_bool))/length(gene1_bool)
    gene2_prob = sum(unlist(gene2_bool))/length(gene2_bool)
    one_prob = rbind(one_prob, c(gene1_prob, gene2_prob))
  }
  colnames(one_prob) = c(gene1, gene2)
  return(one_prob)
}
gt_kni_prob = map(gt_kni_raw, Probability_Conversion, sep = 1, gene1 = 'gt', gene2 = 'kni')
kr_hb_prob = map(kr_hb_raw, Probability_Conversion, sep = 1, gene1 = 'kr', gene2 = 'hb')

# finally assemble the prob data into one complete dataframe
droso_gap_prob = data.frame()
for (i in 1:8){
  prob_temp = cbind('AP' = 4:98, # only consider data in this range because of 'NAN'
                   gt_kni_prob[[i]][4:98, ], kr_hb_prob[[i]][4:98, ], 
                   'time' = rep(i, dim(gt_kni_prob[[i]][4:98, ])[1]))
  droso_gap_prob = rbind(droso_gap_prob, prob_temp)
}

# compare simulation and experimental Prob. results:
source = c(rep('sim', 3040), rep('WT', 3040))
gene_label = rep(c(rep('gt', 760), rep('kr', 760), rep('kni', 760), rep('hb', 760)), 2)
droso_compare_prob = droso_avgsim %>%
  cbind(droso_gap_prob[, c(2, 4, 3, 5)]) %>%
  gather(key = 'gene', value = 'value', gt_avg:hb) %>%
  mutate('source' = source, 'gene_label' = gene_label)

# visulized the probabilities comparison results:
ggplot(droso_compare_prob, aes(x = AP, y = value)) +
  geom_line(alpha = 0.8, aes(linetype = source)) + 
  facet_grid(time~gene_label) +
  xlab("AP position") +
  ylab("Probability") + 
  ggtitle("Simulation vs Experimental in Prob") + 
  theme(text = element_text(size = 15))

#########################################################################################
# IV.  Spacial Decomposition :
#########################################################################################
# -----------------------------------------------------------------------------------------------------
# Experimental EDA spcaially 
# -----------------------------------------------------------------------------------------------------

# first compare the effect of binning, see if binning change the dynamic of the interaction
ggplot(droso_gap_bin, aes(x = gt, y = kni, color = AP)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~time, nrow =2) +
  ggtitle("gt vs kni in time binned data set ") + 
  theme(text = element_text(size = 15))

ggplot(gt_kni_exp, aes(x = gt, y = kni, color = AP)) +
  geom_point(alpha = 0.8, shape =  '.') +
  facet_wrap(~time, nrow = 2) +
  ggtitle("gt vs kni in time raw data set ") + 
  theme(text = element_text(size = 15))

# binned strong mutual repression pairs interaction 
ggplot(droso_gap_bin, aes(x = gt, y = kr, color = AP)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~time, nrow =2) +
  ggtitle("gt vs kr in time binned data set ") + 
  theme(text = element_text(size = 15))

ggplot(droso_gap_bin, aes(x = kni, y = hb, color = AP)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~time, nrow =2) +
  ggtitle("kni vs hb in time binned data set ") + 
  theme(text = element_text(size = 15))



AP_label = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100")


#########################################################################################
# V.  Time Delay :
#########################################################################################
time_label = c('1-2', '2-3', '3-4', '4-5', '5-6', '6-7', '7-8')
gene_lag  = data.frame()
gene = c('gt', 'kr', 'kni', 'hb')
for (g in gene){
  for (t in 1:7){
    gene_t1 = filter(droso_gap_bin, time == t)[g]
    gene_t2 = filter(droso_gap_bin, time == t+1)[g]
    gene_temp = data.frame(gene_t1, gene_t2, rep(time_label[t], 95), rep(g, 95))
    colnames(gene_temp) = c('time1', 'time2', 'time_label', 'gene') 
    gene_lag = rbind(gene_lag, gene_temp)
  }
}

ggplot(gene_lag, aes(x = time1, y = time2))+
  geom_point(alpha = 0.8) + 
  facet_grid(time_label~gene) +
  xlab("Time (t)") +
  ylab("Time (t + 1)") +
  ggtitle("Time Delayed(lag by 1)") +
  theme(text = element_text(size = 15))

 