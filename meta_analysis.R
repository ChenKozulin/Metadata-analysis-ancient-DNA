library(magrittr)
library(ggplot2)
url_to_samples_table <- "https://raw.githubusercontent.com/SPAAM-community/AncientMetagenomeDir/e29eb729e4b5d32b3afb872a7183ff51f6b0dbb5/ancientmetagenome-environmental/samples/ancientmetagenome-environmental_samples.tsv"

samples<-readr::read_tsv(url_to_samples_table)
library (dplyr)
samples_filtered <- samples %>%
  filter(
    grepl("^[0-9.]+$", depth),       # keep only numeric (with possible decimal)
    grepl("^[0-9.]+$", sample_age), # same for sample_age
    site_name!="Unknown"
  ) %>%
  mutate(
    depth = as.numeric(depth),
    sample_age = as.numeric(sample_age)
  )
# Now plot depth against sample_age in a scatterplot to see if there is a potential signal
plot<-ggplot(samples_filtered) +
  geom_point(aes(x = depth, y = sample_age)) 
#Recreate this plot with a log-scaled axis
samples_filtered %>%
  ggplot() +
  geom_point(aes(x = depth, y = sample_age)) +
  scale_y_log10(labels = scales::label_comma()) +
  geom_hline(yintercept = 20000, color = "red")
#5: Filter the dataset to remove all samples that are older than this threshold.
samples_young<-samples_filtered %>%
  filter(sample_age<=20000)
#Recreate the plot from above. The log-scaling can be turned off now
samples_young %>%
  ggplot() +
  geom_point(aes(x = depth, y = sample_age))

cor(samples_young$depth, samples_young$sample_age, method = "pearson")

#Determine the number of sites in the filtered dataset
unique(samples_young$site_name)
length(unique(samples_young$site_name))

# Calculate the number of samples per site with group_by and summarize. 
# Sort the result table by the number of samples with arrange.
samples_per_site <- samples_young %>%
  dplyr::group_by(site_name) %>%
  dplyr::summarise(n = dplyr::n()) %>%
  dplyr::arrange(n)

#Prepare a bar-plot that shows this information, with the sites on the x-axis and 
#the number of samples per site on the y-axis. 
#The bars should be ordered by the number of samples.
samples_per_site_sorted <- samples_per_site %>% arrange(n)
ggplot(samples_per_site_sorted,
       aes(x = reorder(site_name, -n), y = n)) +
  geom_bar(stat = "identity") +
  labs(x = "Site", y = "Number of Samples") +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

#Select only the columns site_name and sample_age and show all rows
samples_site_age <- samples %>%
  select(site_name, sample_age)%>%
  filter(site_name!="Unknown")
print(samples_site_age, n = Inf)

#Further simplify this dataset to only one row per site and add a column 
#(summarize) that shows the distance between min and max age, so the age range per site.
one_row_per_site <- samples_site_age %>%
  group_by(site_name) %>%
  summarise(
    age_range = paste(min(sample_age, na.rm = TRUE),
                      max(sample_age, na.rm = TRUE),
                      sep = "-")
  )%>%
  print(n = Inf)

#Join the sample count per site (as computed above) with the age range per site to get a table with both variables.
sample_count_per_site_orig <- samples_filtered %>%
  group_by(site_name) %>%
  summarise(sample_count = n())
one_row_per_site <- one_row_per_site %>%
  right_join(sample_count_per_site_orig, by = "site_name")

#Calculate the mean sampling interval by dividing the age range by the number of 
#samples and add this information in a new column with mutate.
library(tidyr)
one_row_per_site <- one_row_per_site %>%
  separate(age_range, into = c("min_age", "max_age"), sep = "-", convert = TRUE) %>%
  mutate(
    age_span = max_age - min_age,
    mean_sampling_interval = sample_count / age_span
  )
one_row_per_site <- one_row_per_site %>%
  mutate(age_span = paste(min_age, max_age, sep = "-")) %>%
  select(site_name, age_span, sample_count, mean_sampling_interval)

#Is there a global relationship between sample_age and depth?
#Take the samples_young dataset and recreate the simple scatter plot from above. 
#But now map the site_name to the point colour.
plot2<-ggplot(samples_young) +
  geom_point(aes(x = depth, y = sample_age, color=site_name)) 

#w/o legend
plot2+ 
  guides(color = guide_none())

#Use faceting to split the plot into per-site subplots
plot2+
  facet_wrap(~ site_name)+ 
       guides(color = guide_none())

#adjust the scaling of the subplots
plot2+
  facet_wrap(~ site_name, scales = "free")+ 
  guides(color = guide_none())

#Remove “single-dot” sites and recreate the plot. 
#filter by the standard deviation (sd) along the age or the depth axis.
library(dplyr)

sample_count_per_site <- samples_young %>%
  group_by(site_name) %>%
  summarise(n_samples = n())

samples_no_single <- samples_young %>%
  group_by(site_name) %>%
  filter(n() > 4) %>%   # keep only sites with >1 sample
  ungroup()
ggplot(samples_no_single, aes(x = depth, y = sample_age, color = site_name)) +
  geom_point() +
  facet_wrap(~ site_name, scales = "free")+ 
  guides(color = guide_none())

#this facetted plot visually confirms that depth and sample_age are often correlated 
#for the sites in the environmental samples table of the AncientMetagenomeDir. 
#It also shows a number of notable exceptions that clearly stand out in this plot.


