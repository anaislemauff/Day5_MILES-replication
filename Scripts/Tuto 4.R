### load the packages ###
library(phyloseq)  # the microbiome data container and its tools
library(vegan)     # community ecology: diversity indices and distances
library(dplyr)     # data manipulation: grouping, summarizing
library(tidyr)     # reshaping: pivot_longer
library(ggplot2)   # plotting
library(here)      # file access

## Load dataset
# Tutorial 2 saved `ps` here. These tutorials sit in a subfolder alongside the
# earlier ones, so look in the working directory first, then one level up.
# For this new project, here() [1] "/cloud/project"

ps_file <- here("data","tutorial2_objects.RData")
load(ps_file)
ps                     # print the object's summary to confirm it loaded
sample_variables(ps)   # tell us which metadata columns are available to group by?

#######STEP 1#######

#### Check 1: using the phyloseq package --> gives indications about the dataset

nsamples(ps)   # how many samples?
ntaxa(ps)      # how many features?
# What taxonomic ranks are in the table, and how many features at each?
taxonomy <- as.data.frame(as(tax_table(ps), "matrix"))
table(taxonomy$Rank)
# Are these read counts, or proportions? Counts are whole numbers >= 0.
range(otu_table(ps)) ## this show relative abundance, not raw counts (processed data already?)

#### Check 2: Is there any insulin phenotype, recorded on the participants ###
sample_variables(ps)   # every metadata column attached to these samples 
#[Is that specific command to use in the metadata?]

# The grouping variable the whole comparison rests on
table(sample_data(ps)$ir_is_classification, useNA = "ifany")
# IR      IS Unknown 
# 286     377     192  As there are 192 unknown, we will later exclude them

### Check 3: does this taxonomy actually contain the genera MILES named? ###

# The genera MILES reported, written down before we look at anything
miles_genera <- c("Alistipes", "Coprococcus", "Flavonifractor", "Odoribacter",
                  "Oscillibacter", "Pseudoflavonifractor", "Anaerostipes")
# Every genus name this dataset knows about
hmp2_genera <- taxonomy$Taxon[taxonomy$Rank == "Genus"]

data.frame(Genus = miles_genera, present_in_HMP2 = miles_genera %in% hmp2_genera)
# one genera = FALSE (Anaerostipes) --> going to be described as non detected

### Check 4: is there enough data at a single time point? ###
# Samples per visit, biggest first
head(sort(table(sample_data(ps)$VisitID), decreasing = TRUE), 4)


######STEP 2######
##Going to choose one visit (01=69 most of cases), 
#because the visits may not consitent whithin patients, use the highest number of samples during visit 01
visit01 <- subset_samples(ps, VisitID == "01")
nsamples(visit01)   # 855 samples down to 69

#Select the 2 groups for comparison (IS and IR), exclude unknown
visit01_IRIS <- subset_samples(visit01, ir_is_classification %in% c("IR", "IS"))
nsamples(visit01_IRIS)                                 # 69 down to 44
table(sample_data(visit01_IRIS)$ir_is_classification)  # expect IR 28, IS 16

#Keeping on taxonomic rank (here genus)
visit01_IRIS <- subset_taxa(visit01_IRIS, Rank == "Genus")
ntaxa(visit01_IRIS)   # 96 features down to 45 genera

#Drop the genera that were never seen in the data
sum(taxa_sums(visit01_IRIS) == 0)   # how many all-zero genera?
visit01_IRIS <- prune_taxa(taxa_sums(visit01_IRIS) > 0, visit01_IRIS)
ntaxa(visit01_IRIS) #But here nothing changes

#view of the data: 
#We have 44 stool samples, 28 IR and 16 IS, described by genus-level relative abundances.
visit01_IRIS
summary(sample_sums(visit01_IRIS))   # note: these do NOT sum to 1


#####STEP 3######
### 3.1 Normalize then melt
# Divide each sample by its own total, so every sample sums to 1
visit01_rel <- transform_sample_counts(visit01_IRIS, function(x) x / sum(x))

# Check it worked: every sample total should now be 1
all(abs(sample_sums(visit01_rel) - 1) < 1e-9)

# One row per sample x genus, with metadata and the clean Taxon label attached
df <- psmelt(visit01_rel)
dim(df)   # 44 samples x 45 genera = 1980 rows

### 3.2 Which genera dominate each group?
# Mean relative abundance of each genus within each group
genus_means <- df |>
  group_by(ir_is_classification, Taxon) |>
  summarise(mean_abund = mean(Abundance), .groups = "drop")

ggplot(genus_means, aes(x = ir_is_classification, y = mean_abund, fill = Taxon)) +
  geom_col() +                                  # stacked bar, one column per group
  labs(
    title = "Mean genus-level composition by IR / IS group",
    x     = "IR / IS classification",
    y     = "Mean relative abundance"
  ) +
  theme_minimal() +
  theme(legend.key.size = unit(0.3, "cm"),      # 45 genera need a small legend
        legend.text     = element_text(size = 6))
