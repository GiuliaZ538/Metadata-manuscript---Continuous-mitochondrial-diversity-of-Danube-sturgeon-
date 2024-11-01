#Partition Finder
##partitionfinder v. 2.1.1 was run with python v.2.7, on alignment files saved with extension .phy
python /path/to/partitionfinder-2.1.1/PartitionFinder.py /path/to/folder/alignment_files/

##Settings
## ALIGNMENT FILE ##
alignment = alignment.phy;

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##
models = beast;
##models = GTR, GTR+G, GTR+I+G;

# MODEL SELECCTION: AIC | AICc | BIC #
model_selection = BIC;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]

ND1_1stpos = 1-972\3;
ND1_2ndpos = 2-972\3;
ND1_3rdpos = 3-972\3;
ND2_1stpos = 973-2016\3;
ND2_2ndpos = 974-2016\3;
ND2_3rdpos = 975-2016\3;
COX1_1stpos = 2017-3567\3;
COX1_2ndpos = 2018-3567\3;
COX1_3rdpos = 2019-3567\3;
COX2_1stpos = 3568-4266\3;
COX2_2ndpos = 3569-4266\3;
COX2_3rdpos = 3570-4266\3;
ATP8_1stpos = 4267-4433\3;
ATP8_2ndpos = 4268-4433\3;
ATP8_3rdpos = 4269-4433\3;
ATP6_1stpos = 4434-5114\3;
ATP6_2ndpos = 4435-5114\3;
ATP6_3rdpos = 4436-5114\3;
COX3_1stpos = 5115-5897\3;
COX3_2ndpos = 5116-5897\3;
COX3_3rdpos = 5117-5897\3;
ND3_1stpos = 5898-6245\3;
ND3_2ndpos = 5899-6245\3;
ND3_3rdpos = 5900-6245\3;
ND4L_1stpos = 6246-6539\3; 
ND4L_2ndpos = 6247-6539\3; 
ND4L_3rdpos = 6248-6539\3; 
ND4_1stpos = 6540-7920\3;
ND4_2ndpos = 6541-7920\3;
ND4_3rdpos = 6542-7920\3;
ND5_1stpos = 7921-9759\3;
ND5_2ndpos = 7922-9759\3;
ND5_3rdpos = 7923-9759\3;
ND6_1stpos = 9760-10277\3;
ND6_2ndpos = 9761-10277\3;
ND6_3rdpos = 9762-10277\3;
CYTB_1stpos = 10278-11419\3;
CYTB_2ndpos = 10279-11419\3;
CYTB_3rdpos = 10280-11419\3;

## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##
[schemes]
search = greedy;
